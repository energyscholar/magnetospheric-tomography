"""Canonical spatial autocorrelation statistics.

Consolidates 15+ ad-hoc Moran's I implementations across the codebase
into a single, tested module. Row-normalized inverse-distance weights
by default (PySAL convention).
"""

import numpy as np


def haversine_km(lat1, lon1, lat2, lon2):
    """Great-circle distance in km between two points (decimal degrees)."""
    R = 6371.0
    la1, lo1, la2, lo2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dlat, dlon = la2 - la1, lo2 - lo1
    a = np.sin(dlat/2)**2 + np.cos(la1)*np.cos(la2)*np.sin(dlon/2)**2
    return R * 2 * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def distance_matrix(lats, lons):
    """(N, N) pairwise distance matrix in km from coordinate arrays."""
    lats = np.asarray(lats)
    lons = np.asarray(lons)
    la = np.radians(lats)
    lo = np.radians(lons)
    dlat = la[:, None] - la[None, :]
    dlon = lo[:, None] - lo[None, :]
    a = np.sin(dlat/2)**2 + np.cos(la[:, None]) * np.cos(la[None, :]) * np.sin(dlon/2)**2
    return 6371.0 * 2 * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def inverse_distance_weights(dist, cutoff_km=None, row_normalize=True):
    """Spatial weight matrix from distance matrix.

    W[i,j] = 1/d[i,j] for i != j (and d > 0).
    Row-normalized by default so each station contributes equally
    regardless of neighbor count.
    """
    with np.errstate(divide='ignore'):
        W = np.where(dist > 0, 1.0 / dist, 0.0)
    if cutoff_km is not None:
        W[dist > cutoff_km] = 0.0
    np.fill_diagonal(W, 0.0)
    if row_normalize:
        row_sums = W.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1.0
        W = W / row_sums
    return W


def queen_contiguity_weights(n_lat, n_lon, lat_mask=None):
    """8-neighbor queen contiguity weight matrix for a regular grid.

    lat_mask: boolean array of shape (n_lat,) -- True for included latitudes.
    Returns (W, cell_indices) where W is (n_cells, n_cells) and
    cell_indices maps flat index to (lat_idx, lon_idx).
    """
    if lat_mask is None:
        lat_mask = np.ones(n_lat, dtype=bool)

    cell_map = {}
    idx = 0
    for i in range(n_lat):
        if not lat_mask[i]:
            continue
        for j in range(n_lon):
            cell_map[(i, j)] = idx
            idx += 1

    n_cells = len(cell_map)
    W = np.zeros((n_cells, n_cells))

    for (i, j), ci in cell_map.items():
        for di in [-1, 0, 1]:
            for dj in [-1, 0, 1]:
                if di == 0 and dj == 0:
                    continue
                ni = i + di
                nj = (j + dj) % n_lon  # wrap longitude
                if (ni, nj) in cell_map:
                    W[ci, cell_map[(ni, nj)]] = 1.0

    row_sums = W.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1.0
    W = W / row_sums

    cell_indices = sorted(cell_map.keys(), key=lambda k: cell_map[k])
    return W, cell_indices


def morans_i(values, W, min_valid=5):
    """Canonical Moran's I with NaN handling.

    Returns (I, expected_I, z_score, n_valid).
    NaN stations are excluded and W is recomputed on the valid subset.
    """
    values = np.asarray(values, dtype=float)
    valid = np.isfinite(values)
    n_valid = int(valid.sum())

    if n_valid < min_valid:
        return np.nan, np.nan, np.nan, n_valid

    if valid.all():
        v = values
        Wv = W
    else:
        idx = np.where(valid)[0]
        v = values[idx]
        Wv = W[np.ix_(idx, idx)]
        row_sums = Wv.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1.0
        Wv = Wv / row_sums

    n = len(v)
    mean_v = v.mean()
    dev = v - mean_v
    ss = np.sum(dev**2)
    if ss < 1e-20:
        return np.nan, -1.0/(n-1), np.nan, n_valid

    S0 = Wv.sum()
    if S0 < 1e-20:
        return np.nan, -1.0/(n-1), np.nan, n_valid

    numer = n * np.sum(Wv * np.outer(dev, dev))
    I = numer / (S0 * ss)

    E_I = -1.0 / (n - 1)

    S1 = 0.5 * np.sum((Wv + Wv.T)**2)
    S2 = np.sum((Wv.sum(axis=1) + Wv.sum(axis=0))**2)
    k = (n * np.sum(dev**4)) / (ss**2)
    var_I = (n * ((n**2 - 3*n + 3)*S1 - n*S2 + 3*S0**2) -
             k * (n * (n - 1)*S1 - 2*n*S2 + 6*S0**2)) / \
            ((n - 1) * (n - 2) * (n - 3) * S0**2) - E_I**2

    z = (I - E_I) / np.sqrt(max(var_I, 1e-20)) if var_I > 0 else np.nan

    return I, E_I, z, n_valid


def morans_i_timeseries(data, W, min_valid=5):
    """Compute Moran's I over time axis.

    data: (n_times, n_stations) array.
    Returns dict with keys: I, E_I, z, n_valid (all shape (n_times,)).

    Uses vectorized path when all values are finite (common for SWMF).
    """
    data = np.asarray(data, dtype=float)
    n_times, n = data.shape

    if n < min_valid:
        return {
            'I': np.full(n_times, np.nan),
            'E_I': np.full(n_times, np.nan),
            'z': np.full(n_times, np.nan),
            'n_valid': np.zeros(n_times, dtype=int),
        }

    has_nan = np.any(~np.isfinite(data))

    if not has_nan:
        means = data.mean(axis=1, keepdims=True)
        dev = data - means
        ss_vec = np.sum(dev**2, axis=1)
        S0 = W.sum()
        numer = n * np.einsum('ti,ij,tj->t', dev, W, dev)
        E_I = -1.0 / (n - 1)

        I_arr = np.where(ss_vec > 1e-20, numer / (S0 * ss_vec), np.nan)

        S1 = 0.5 * np.sum((W + W.T)**2)
        S2 = np.sum((W.sum(axis=1) + W.sum(axis=0))**2)
        k_vec = (n * np.sum(dev**4, axis=1)) / np.where(ss_vec > 1e-20, ss_vec**2, 1.0)
        var_I = (n * ((n**2 - 3*n + 3)*S1 - n*S2 + 3*S0**2) -
                 k_vec * (n * (n - 1)*S1 - 2*n*S2 + 6*S0**2)) / \
                ((n - 1) * (n - 2) * (n - 3) * S0**2) - E_I**2
        z_arr = np.where(var_I > 0, (I_arr - E_I) / np.sqrt(np.maximum(var_I, 1e-20)), np.nan)

        return {
            'I': I_arr,
            'E_I': np.full(n_times, E_I),
            'z': z_arr,
            'n_valid': np.full(n_times, n, dtype=int),
        }

    result = {
        'I': np.full(n_times, np.nan),
        'E_I': np.full(n_times, np.nan),
        'z': np.full(n_times, np.nan),
        'n_valid': np.zeros(n_times, dtype=int),
    }
    for t in range(n_times):
        I, E_I, z, nv = morans_i(data[t], W, min_valid=min_valid)
        result['I'][t] = I
        result['E_I'][t] = E_I
        result['z'][t] = z
        result['n_valid'][t] = nv
    return result


def permutation_test(values, W, n_perms=999):
    """Monte Carlo permutation test for Moran's I significance.

    Returns (I_observed, p_value).
    """
    values = np.asarray(values, dtype=float)
    valid = np.isfinite(values)
    if valid.sum() < 5:
        return np.nan, np.nan

    I_obs = morans_i(values, W)[0]
    if np.isnan(I_obs):
        return np.nan, np.nan

    rng = np.random.default_rng(42)
    count_ge = 0
    v = values[valid].copy()
    idx = np.where(valid)[0]
    Wv = W[np.ix_(idx, idx)]

    for _ in range(n_perms):
        rng.shuffle(v)
        I_perm = morans_i(v, Wv)[0]
        if not np.isnan(I_perm) and abs(I_perm) >= abs(I_obs):
            count_ge += 1

    p = (count_ge + 1) / (n_perms + 1)
    return I_obs, p
