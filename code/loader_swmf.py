"""SWMF virtual magnetometer data loader for the unified Moran's I pipeline."""

import numpy as np
from pathlib import Path

BASE_DIR = Path(__file__).parent.parent
CCMC_DIR = BASE_DIR / 'data' / 'ccmc_vmag' / 'SWPCTEST' / 'deltaB' / 'IGRFDIPOLE'

EVENTS = {
    1: {'date': '2003-10-29', 'name': 'Halloween 2003',
        'solar_phase': 'max', 'severity': 'extreme',
        'ym': '200310', 'day0': 28},
    2: {'date': '2006-12-14', 'name': 'Dec 2006',
        'solar_phase': 'declining', 'severity': 'strong',
        'ym': '200612', 'day0': 13},
    3: {'date': '2001-08-31', 'name': 'Aug 2001',
        'solar_phase': 'max', 'severity': 'strong',
        'ym': '200108', 'day0': 30},
    4: {'date': '2005-08-31', 'name': 'Aug 2005',
        'solar_phase': 'declining', 'severity': 'moderate',
        'ym': '200508', 'day0': 30},
    5: {'date': '2010-04-05', 'name': 'Apr 2010',
        'solar_phase': 'min', 'severity': 'moderate',
        'ym': '201004', 'day0': 4},
    6: {'date': '2011-08-05', 'name': 'Aug 2011',
        'solar_phase': 'ascending', 'severity': 'moderate',
        'ym': '201108', 'day0': 4},
}

COMPONENT_MAP = {
    'total_N': 0, 'total_E': 1, 'total_D': 2,
    'mhd_N': 3, 'mhd_E': 4, 'mhd_D': 5,
    'fac_N': 6, 'fac_E': 7, 'fac_D': 8,
    'hall_N': 9, 'hall_E': 10, 'hall_D': 11,
    'ped_N': 12, 'ped_E': 13, 'ped_D': 14,
}


def find_mag_file(event_id):
    """Find the .mag file for a CCMC SWPCTEST event."""
    event_dir = CCMC_DIR / f'Event{event_id}'
    if not event_dir.exists():
        return None
    mag_files = list(event_dir.glob('magnetometers_*.mag'))
    if not mag_files:
        return None
    return mag_files[0]


def parse_mag_file(filepath):
    """Parse SWMF .mag file.

    Returns dict with:
        station_names: list of station codes (excluding DST)
        positions: {station: (lat, lon)} in geographic degrees
        times: array of fractional hours
        data: {station: (n_times, 15) array of dB values}
    """
    stations = {}
    station_positions = {}
    times = []
    idx_to_name = {}

    with open(filepath) as f:
        for line in f:
            parts = line.split()

            if 'magnetometers:' in line:
                colon_idx = parts.index('magnetometers:')
                names = parts[colon_idx + 1:]
                for i, name in enumerate(names):
                    idx_to_name[str(i + 1)] = name
                continue

            if len(parts) < 27:
                continue
            try:
                int(parts[0])
            except ValueError:
                continue

            hr, mn, sc = int(parts[4]), int(parts[5]), int(parts[6])
            sta_idx = parts[8]
            sta = idx_to_name.get(sta_idx, sta_idx)
            x, y, z = float(parts[9]), float(parts[10]), float(parts[11])
            vals = [float(v) for v in parts[12:27]]

            if sta not in stations:
                stations[sta] = []
                r = np.sqrt(x**2 + y**2 + z**2)
                lat = np.degrees(np.arcsin(z / r)) if r > 0 else 0
                lon = np.degrees(np.arctan2(y, x))
                if lon < 0:
                    lon += 360
                station_positions[sta] = (lat, lon)

            stations[sta].append(vals)

            if sta == list(stations.keys())[0]:
                frac_hour = hr + mn/60.0 + sc/3600.0
                if times and frac_hour < times[-1] - 12:
                    frac_hour += 24
                times.append(frac_hour)

    for sta in stations:
        stations[sta] = np.array(stations[sta])

    real_names = [s for s in stations if s != 'DST']

    return {
        'station_names': real_names,
        'positions': {s: station_positions[s] for s in real_names},
        'times': np.array(times),
        'data': {s: stations[s] for s in real_names},
    }


def load_event(event_id):
    """Load and parse a specific SWMF CCMC event. Returns None if not on disk."""
    mag_file = find_mag_file(event_id)
    if mag_file is None:
        return None
    info = EVENTS.get(event_id, {})
    parsed = parse_mag_file(mag_file)
    parsed['event_id'] = event_id
    parsed['event_info'] = info
    return parsed


def extract_component(data, stations, comp_name):
    """Extract a single component timeseries for all stations.

    Returns (n_times, n_stations) array.
    """
    col = COMPONENT_MAP[comp_name]
    return np.column_stack([data[s][:, col] for s in stations])


def extract_horizontal_mag(data, stations, prefix='total'):
    """Extract horizontal magnitude sqrt(N^2 + E^2) for a source.

    prefix: 'total', 'mhd', 'fac', 'hall', 'ped'
    Returns (n_times, n_stations) array.
    """
    n_col = COMPONENT_MAP[f'{prefix}_N']
    e_col = COMPONENT_MAP[f'{prefix}_E']
    vals_n = np.column_stack([data[s][:, n_col] for s in stations])
    vals_e = np.column_stack([data[s][:, e_col] for s in stations])
    return np.sqrt(vals_n**2 + vals_e**2)


def aggregate_to_2hour(times, data_2d):
    """Downsample (n_times, n_stations) to 2-hour bins.

    times: fractional hours array.
    Returns (binned_data, bin_centers).
    """
    hour_bins = np.arange(0, 26, 2)
    n_bins = len(hour_bins) - 1
    n_sta = data_2d.shape[1]
    out = np.full((n_bins, n_sta), np.nan)
    centers = np.full(n_bins, np.nan)

    for b in range(n_bins):
        mask = (times >= hour_bins[b]) & (times < hour_bins[b+1])
        if mask.sum() > 0:
            out[b] = np.nanmean(data_2d[mask], axis=0)
            centers[b] = (hour_bins[b] + hour_bins[b+1]) / 2

    return out, centers


def list_available_events():
    """Return dict of available events with verified file paths."""
    available = {}
    for eid, info in EVENTS.items():
        mag_file = find_mag_file(eid)
        if mag_file is not None:
            available[eid] = {**info, 'file': str(mag_file)}
    return available
