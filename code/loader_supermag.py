"""SuperMAG data loader for the unified Moran's I pipeline."""

import json
import numpy as np
from pathlib import Path
from datetime import date

BASE_DIR = Path(__file__).parent.parent
DATA_DIR = BASE_DIR / 'data' / 'supermag'
FILL_VAL = 99999.0


def load_station_coords():
    """Load station coordinates. Returns {CODE: {geolat, geolon, ...}}."""
    with open(DATA_DIR / 'station_coords.json') as f:
        raw = json.load(f)
    out = {}
    for k, v in raw.items():
        code = k.upper()
        out[code] = {
            'geolat': v.get('geolat', v.get('glat', np.nan)),
            'geolon': v.get('geolon', v.get('glon', np.nan)),
            'name': v.get('name', code),
        }
        if 'aacgmlat' in v:
            out[code]['maglat'] = v['aacgmlat']
            out[code]['maglon'] = v['aacgmlon']
    return out


def load_station_ranges():
    """Load station date ranges. Returns {code: {start, stop}}."""
    p = DATA_DIR / 'station_ranges.json'
    if not p.exists():
        return {}
    with open(p) as f:
        return json.load(f)


def available_stations_for_date(year, month):
    """Return list of station codes that have data for (year, month)."""
    ym = f"{year:04d}{month:02d}"
    available = []
    for f in DATA_DIR.glob(f"*_{ym}_full.npz"):
        code = f.stem.split('_')[0].upper()
        available.append(code)
    return sorted(available)


def load_station_month(station, year, month):
    """Load one station-month NPZ. Returns (n_minutes, 3) NEZ array or None."""
    code = station.lower()
    ym = f"{year:04d}{month:02d}"
    fpath = DATA_DIR / f"{code}_{ym}_full.npz"
    if not fpath.exists():
        return None
    try:
        d = np.load(fpath, allow_pickle=True)
        if 'nez' in d:
            nez = d['nez'].astype(np.float64)
        elif 'N' in d and 'E' in d and 'Z' in d:
            nez = np.column_stack([
                d['N'].astype(np.float64),
                d['E'].astype(np.float64),
                d['Z'].astype(np.float64),
            ])
        else:
            return None
    except Exception:
        return None
    nez[nez >= FILL_VAL - 1] = np.nan
    nez[nez <= -(FILL_VAL - 1)] = np.nan
    return nez


def load_event_day(year, month, day, n_days=1, stations=None,
                   min_valid_frac=0.3):
    """Load all available stations for an event period.

    Returns {CODE: (n_minutes, 3) NEZ array}, station_coords dict.
    """
    if stations is None:
        stations = available_stations_for_date(year, month)

    coords = load_station_coords()
    day0 = day - 1  # 0-indexed

    result = {}
    result_coords = {}
    for sta in stations:
        code = sta.upper()
        nez = load_station_month(code, year, month)
        if nez is None:
            continue

        start_row = day0 * 1440
        end_row = (day0 + n_days) * 1440
        if start_row >= len(nez):
            continue
        end_row = min(end_row, len(nez))

        chunk = nez[start_row:end_row].copy()

        valid_frac = np.isfinite(chunk[:, 0]).mean()
        if valid_frac < min_valid_frac:
            continue

        if code not in coords:
            continue

        result[code] = chunk
        result_coords[code] = coords[code]

    return result, result_coords


def aggregate_to_2hour(minute_data, method='median'):
    """Downsample (n_minutes, 3) to (n_2hr_bins, 3).

    Each 2-hour bin = 120 minutes. method: 'median' or 'mean'.
    """
    n_min = minute_data.shape[0]
    n_cols = minute_data.shape[1]
    n_bins = n_min // 120
    if n_bins == 0:
        return np.full((1, n_cols), np.nan)

    reshaped = minute_data[:n_bins * 120].reshape(n_bins, 120, n_cols)
    valid_counts = np.isfinite(reshaped).sum(axis=1)
    agg_fn = np.nanmedian if method == 'median' else np.nanmean
    with np.errstate(all='ignore'):
        out = agg_fn(reshaped, axis=1)
    out[valid_counts < 10] = np.nan
    return out
