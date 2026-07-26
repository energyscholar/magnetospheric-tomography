"""Relational diagnostics for magnetospheric tomography pipeline.

Computes edge-field metrics, altitude-resolved Moran's I, cross-component
coupling, and variogram shape for SuperMAG + GPS TEC data.

Usage:
    import relational_diagnostics as rd
    metrics = rd.epoch_metrics(sm_data, sm_coords, tec_grid, epoch_minutes=(690, 750))
    rd.print_summary(metrics)
"""

import numpy as np
from math import radians, sin, cos, asin, sqrt

import loader_supermag as sm
import loader_tec as tec
import spatial_stats as ss

RE = 6371.0
CHAPMAN_INT = 4.133
H_SCALE = 60.0


def haversine_km(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = map(radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2
    return 2 * RE * asin(min(1, sqrt(a)))


def chapman_ne(n_max, z_km, h_max=300, H=60):
    zeta = (z_km - h_max) / H
    return n_max * np.exp(0.5 * (1 - zeta - np.exp(-zeta)))


def station_distances(lats, lons):
    """Build great-circle distance matrix (km)."""
    n = len(lats)
    dist = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = haversine_km(lats[i], lons[i], lats[j], lons[j])
            dist[i, j] = dist[j, i] = d
    return dist


def edge_field_metrics(vals_N, vals_E, dist, W):
    """Compute relational metrics from edge fields (A operator output).

    Returns dict with: I_N, I_E, I_Z, gamma_N, gamma_E,
    coupling_pct, vario_ratio.
    """
    eN = vals_N[:, None] - vals_N[None, :]
    eE = vals_E[:, None] - vals_E[None, :]

    gamma_N = np.mean(eN ** 2) / 2
    gamma_E = np.mean(eE ** 2) / 2

    # Cross-component coupling (R operator at rho=0.3)
    rho = 0.3
    aN = (eN - eN.T) / 2
    aE = (eE - eE.T) / 2
    var_coupled = (np.mean((eN + rho * aE) ** 2) + np.mean((eE + rho * aN) ** 2)) / 2
    var_uncoupled = (np.mean(eN ** 2) + np.mean(eE ** 2)) / 2
    coupling_pct = (var_coupled - var_uncoupled) / var_uncoupled * 100 if var_uncoupled > 0 else 0

    I_N, _, _, _ = ss.morans_i(vals_N, W)
    I_E, _, _, _ = ss.morans_i(vals_E, W)

    # Variogram shape: far/near semivariance ratio
    n = len(vals_N)
    upper = np.triu_indices(n, k=1)
    dists_u = dist[upper]
    diffs2_N = (vals_N[upper[0]] - vals_N[upper[1]]) ** 2
    med_d = np.median(dists_u)
    near = diffs2_N[dists_u < med_d]
    far = diffs2_N[dists_u >= med_d]
    near_gamma = np.mean(near) / 2 if len(near) > 0 else 1e-10
    far_gamma = np.mean(far) / 2 if len(far) > 0 else 0
    vario_ratio = far_gamma / near_gamma if near_gamma > 0 else float('inf')

    return {
        'I_N': I_N, 'I_E': I_E,
        'gamma_N': gamma_N, 'gamma_E': gamma_E,
        'coupling_pct': coupling_pct,
        'vario_ratio': vario_ratio,
    }


def altitude_resolved_I(n_max_arr, lats, W, auroral_mask, f_auroral=0.4,
                        H_solar=60.0, H_precip=20.0):
    """Compute Moran's I at multiple altitudes using two-layer Chapman model.

    Solar EUV layer: h_max=300 km, H=H_solar km (global)
    Auroral precipitation layer: h_max=120 km, H=H_precip km (auroral zone only)

    TEC-conserving: f_auroral is fraction of column TEC from precipitation.
    n_precip = f * n_max * (H_solar/H_precip) so that the column integral
    n_solar*H_solar*C + n_precip*H_precip*C = n_max*H_solar*C = TEC_input.

    Returns dict mapping altitude_km -> I_ne.
    """
    n = len(n_max_arr)
    ne_solar = n_max_arr.copy()
    ne_precip = np.zeros(n)

    for i in range(n):
        if auroral_mask[i]:
            ne_solar[i] = n_max_arr[i] * (1 - f_auroral)
            ne_precip[i] = n_max_arr[i] * f_auroral * (H_solar / H_precip)

    altitudes = [100, 120, 150, 200, 300, 400, 500]
    result = {}
    for alt in altitudes:
        ne = np.array([
            chapman_ne(ne_solar[i], alt, 300, H_solar) +
            chapman_ne(ne_precip[i], alt, 120, H_precip)
            for i in range(n)
        ])
        I_ne, _, _, _ = ss.morans_i(ne, W)
        result[alt] = I_ne

    return result


def epoch_metrics(sm_data, sm_coords, tec_grid=None,
                  epoch_minutes=None, epoch_2hr=None):
    """Full relational analysis for one epoch.

    Specify time window via epoch_minutes=(start, end) for minute-level data,
    or epoch_2hr=idx for 2-hour bin index.

    Returns dict with ground metrics + optional altitude-resolved metrics.
    """
    codes = sorted([c for c in sm_data if c in sm_coords])
    lats = np.array([sm_coords[c]['geolat'] for c in codes])
    lons = np.array([sm_coords[c]['geolon'] for c in codes])

    # Extract values
    if epoch_minutes is not None:
        t0, t1 = epoch_minutes
        vals = {}
        for ci, cn in enumerate(['N', 'E', 'Z']):
            v = []
            for c in codes:
                chunk = sm_data[c][t0:t1, ci]
                if np.isfinite(chunk).mean() > 0.3:
                    v.append(float(np.nanmedian(chunk)))
                else:
                    v.append(np.nan)
            vals[cn] = np.array(v)
    elif epoch_2hr is not None:
        vals = {}
        for ci, cn in enumerate(['N', 'E', 'Z']):
            v = []
            for c in codes:
                agg = sm.aggregate_to_2hour(sm_data[c], method='median')
                if epoch_2hr < len(agg) and np.isfinite(agg[epoch_2hr, ci]):
                    v.append(agg[epoch_2hr, ci])
                else:
                    v.append(np.nan)
            vals[cn] = np.array(v)
    else:
        raise ValueError("Specify epoch_minutes or epoch_2hr")

    valid = np.isfinite(vals['N']) & np.isfinite(vals['E']) & np.isfinite(vals['Z'])
    n = int(np.sum(valid))
    if n < 8:
        return None

    dist = station_distances(lats[valid], lons[valid])
    W = ss.inverse_distance_weights(dist)

    # Ground-level edge-field metrics
    metrics = edge_field_metrics(vals['N'][valid], vals['E'][valid], dist, W)
    I_Z, _, _, _ = ss.morans_i(vals['Z'][valid], W)
    metrics['I_Z'] = I_Z
    metrics['n_stations'] = n
    metrics['rms_N'] = float(np.sqrt(np.mean(vals['N'][valid] ** 2)))

    # Altitude-resolved metrics (requires TEC)
    if tec_grid is not None:
        sta_coords_dict = {c: sm_coords[c] for c in [codes[i] for i in range(len(codes)) if valid[i]]}
        tec_at = tec.extract_at_stations(tec_grid, sta_coords_dict)

        n_max_arr = []
        tec_mask = []
        valid_codes = [codes[i] for i in range(len(codes)) if valid[i]]
        for c in valid_codes:
            if c in tec_at and tec_at[c] is not None and tec_at[c] > 0:
                n_max_arr.append(tec_at[c] * 1e16 / (H_SCALE * 1e3 * CHAPMAN_INT))
                tec_mask.append(True)
            else:
                n_max_arr.append(np.nan)
                tec_mask.append(False)

        n_max_arr = np.array(n_max_arr)
        tec_mask = np.array(tec_mask)
        n_tec = np.sum(tec_mask)

        if n_tec >= 8:
            lats_v = lats[valid][tec_mask]
            W_tec = ss.inverse_distance_weights(dist[np.ix_(tec_mask, tec_mask)])
            auroral = np.abs(lats_v) > 60

            # Storm-dependent auroral fraction
            rms = metrics['rms_N']
            f_aur = min(0.7, max(0.1, 0.1 + 0.6 * (rms / 1000)))

            alt_I = altitude_resolved_I(n_max_arr[tec_mask], lats_v, W_tec, auroral, f_aur)
            metrics['altitude_I'] = alt_I
            metrics['I_ne'] = alt_I.get(300, np.nan)
            metrics['delta_120_300'] = alt_I.get(120, 0) - alt_I.get(300, 0)
            metrics['f_auroral'] = f_aur
            metrics['n_tec_stations'] = int(n_tec)

    return metrics


def multi_epoch_table(sm_data, sm_coords, tec_maps=None, resolution='2hr'):
    """Compute metrics across all epochs for an event day.

    resolution: '2hr' for 2-hour bins, '30min' for 30-minute bins
    Returns list of (ut_hours, metrics_dict) tuples.
    """
    results = []
    if resolution == '2hr':
        for epoch in range(12):
            tec_grid = tec_maps[min(epoch, len(tec_maps) - 1)][1] if tec_maps else None
            m = epoch_metrics(sm_data, sm_coords, tec_grid, epoch_2hr=epoch)
            results.append((epoch * 2, m))
    elif resolution == '30min':
        for epoch in range(48):
            t0 = epoch * 30
            t1 = (epoch + 1) * 30
            tec_idx = epoch // 4 if tec_maps else 0
            tec_grid = tec_maps[min(tec_idx, len(tec_maps) - 1)][1] if tec_maps else None
            m = epoch_metrics(sm_data, sm_coords, tec_grid, epoch_minutes=(t0, t1))
            hours = epoch * 0.5
            results.append((hours, m))
    return results


def print_summary(results, show_altitude=True):
    """Print formatted table of multi-epoch results."""
    has_alt = any(m and 'altitude_I' in m for _, m in results)

    if has_alt and show_altitude:
        hdr = f"{'UT':>5} {'RMS':>5} {'I_N':>5} {'I_E':>5} {'I_Z':>5} {'G%':>6} {'VR':>5} {'I120':>5} {'I300':>5} {'D120':>5}"
    else:
        hdr = f"{'UT':>5} {'RMS':>5} {'I_N':>5} {'I_E':>5} {'I_Z':>5} {'G%':>6} {'VR':>5}"
    print(hdr)
    print('-' * len(hdr))

    for ut, m in results:
        if m is None:
            h = int(ut)
            mn = int((ut - h) * 60)
            print(f"{h:02d}:{mn:02d}  -- insufficient data --")
            continue

        h = int(ut)
        mn = int((ut - h) * 60)

        if has_alt and show_altitude and 'altitude_I' in m:
            print(f"{h:02d}:{mn:02d} {m['rms_N']:5.0f} {m['I_N']:+.3f} {m['I_E']:+.3f} {m['I_Z']:+.3f} "
                  f"{m['coupling_pct']:+5.1f}% {m['vario_ratio']:5.2f} "
                  f"{m['altitude_I'].get(120, 0):+.3f} {m['altitude_I'].get(300, 0):+.3f} "
                  f"{m.get('delta_120_300', 0):+.3f}")
        else:
            print(f"{h:02d}:{mn:02d} {m['rms_N']:5.0f} {m['I_N']:+.3f} {m['I_E']:+.3f} {m['I_Z']:+.3f} "
                  f"{m['coupling_pct']:+5.1f}% {m['vario_ratio']:5.2f}")
