"""GPS TEC (CODE GIM) data loader for the unified Moran's I pipeline.

Parses IONEX files, extracts TEC grids, provides station-level interpolation
and grid-based Moran's I computation.
"""

import subprocess
import tempfile
import numpy as np
from pathlib import Path
from datetime import datetime, date

BASE_DIR = Path(__file__).parent.parent
GIM_DIR = BASE_DIR / 'data' / 'code_gims'

TEC_LATS = np.arange(87.5, -90.0, -2.5)  # 71 values, N to S
TEC_LONS = np.arange(-180.0, 185.0, 5.0)  # 73 values


def _doy(year, month, day):
    return (date(year, month, day) - date(year, 1, 1)).days + 1


def _find_ionex_file(year, doy):
    """Find IONEX file path on disk (compressed or not)."""
    yy = year % 100
    fname_z = f"CODG{doy:03d}0.{yy:02d}I.Z"
    fname_plain = f"CODG{doy:03d}0.{yy:02d}I"

    year_dir = GIM_DIR / str(year)
    if not year_dir.exists():
        return None

    for name in [fname_plain, fname_z]:
        p = year_dir / name
        if p.exists():
            return p
    return None


def _decompress_z(filepath):
    """Decompress a .Z file to a temp location. Returns path to decompressed file."""
    if not str(filepath).endswith('.Z'):
        return str(filepath)

    tmp = tempfile.NamedTemporaryFile(suffix='.ionex', delete=False)
    tmp.close()

    try:
        with open(filepath, 'rb') as f_in:
            result = subprocess.run(['zcat'], stdin=f_in, capture_output=True, timeout=10)
            if result.returncode == 0:
                with open(tmp.name, 'wb') as f_out:
                    f_out.write(result.stdout)
                return tmp.name
    except Exception:
        pass

    try:
        import shutil
        shutil.copy2(filepath, tmp.name + '.Z')
        subprocess.run(['uncompress', '-f', tmp.name + '.Z'],
                       capture_output=True, timeout=10)
        if Path(tmp.name).exists():
            return tmp.name
    except Exception:
        pass

    return None


def parse_ionex(filepath):
    """Parse one IONEX file.

    Returns list of (epoch_datetime, tec_grid) tuples.
    tec_grid: (71, 73) array in TECU. 9999 values → NaN.
    """
    fpath = str(filepath)
    temp_path = None

    if fpath.endswith('.Z'):
        temp_path = _decompress_z(fpath)
        if temp_path is None:
            return []
        fpath = temp_path

    try:
        with open(fpath) as f:
            lines = f.readlines()
    except Exception:
        return []
    finally:
        if temp_path and Path(temp_path).exists():
            try:
                Path(temp_path).unlink()
            except Exception:
                pass

    exponent = -1
    header_end = 0
    n_lat = 71
    n_lon = 73

    for i, line in enumerate(lines):
        if 'EXPONENT' in line:
            try:
                exponent = int(line[:60].strip())
            except ValueError:
                pass
        elif 'END OF HEADER' in line:
            header_end = i + 1
            break

    scale = 10.0 ** exponent
    maps = []
    i = header_end

    while i < len(lines):
        line = lines[i]
        if 'START OF TEC MAP' in line:
            i += 1
            parts = lines[i][:60].split()
            try:
                epoch = datetime(int(parts[0]), int(parts[1]), int(parts[2]),
                                 int(parts[3]), int(parts[4]), int(parts[5]))
            except (ValueError, IndexError):
                i += 1
                continue
            i += 1

            grid = np.full((n_lat, n_lon), np.nan)
            lat_idx = 0

            while i < len(lines) and 'END OF TEC MAP' not in lines[i]:
                if 'LAT/LON1/LON2/DLON/H' in lines[i]:
                    values = []
                    i += 1
                    while (i < len(lines) and
                           'LAT/LON1/LON2/DLON/H' not in lines[i] and
                           'END OF TEC MAP' not in lines[i] and
                           'START OF TEC MAP' not in lines[i]):
                        row_vals = lines[i][:60].split()
                        values.extend([float(v) for v in row_vals])
                        i += 1

                    vals_arr = np.array(values[:n_lon])
                    vals_arr[vals_arr == 9999] = np.nan
                    if lat_idx < n_lat:
                        grid[lat_idx, :len(vals_arr)] = vals_arr * scale
                    lat_idx += 1
                    continue
                i += 1

            maps.append((epoch, grid))
        else:
            i += 1

    return maps


def load_day(year, month, day):
    """Load all TEC maps for a given date. Returns list of (datetime, grid) or None."""
    doy = _doy(year, month, day)
    fpath = _find_ionex_file(year, doy)
    if fpath is None:
        return None
    return parse_ionex(fpath)


def bilinear_interpolate(grid, lat, lon):
    """Interpolate TEC at arbitrary (lat, lon) from a (71, 73) grid."""
    if lon < TEC_LONS[0]:
        lon += 360
    elif lon > TEC_LONS[-1]:
        lon -= 360

    lat_idx = np.searchsorted(-TEC_LATS, -lat) - 1
    lon_idx = np.searchsorted(TEC_LONS, lon) - 1
    lat_idx = max(0, min(lat_idx, len(TEC_LATS) - 2))
    lon_idx = max(0, min(lon_idx, len(TEC_LONS) - 2))

    dlat = TEC_LATS[lat_idx + 1] - TEC_LATS[lat_idx]
    dlon = TEC_LONS[lon_idx + 1] - TEC_LONS[lon_idx]

    lat_frac = (lat - TEC_LATS[lat_idx]) / dlat if abs(dlat) > 1e-10 else 0
    lon_frac = (lon - TEC_LONS[lon_idx]) / dlon if abs(dlon) > 1e-10 else 0
    lat_frac = max(0.0, min(1.0, lat_frac))
    lon_frac = max(0.0, min(1.0, lon_frac))

    v00 = grid[lat_idx, lon_idx]
    v01 = grid[lat_idx, lon_idx + 1]
    v10 = grid[lat_idx + 1, lon_idx]
    v11 = grid[lat_idx + 1, lon_idx + 1]

    if any(np.isnan(v) for v in [v00, v01, v10, v11]):
        return np.nan

    return (v00 * (1-lat_frac) * (1-lon_frac) +
            v01 * (1-lat_frac) * lon_frac +
            v10 * lat_frac * (1-lon_frac) +
            v11 * lat_frac * lon_frac)


def extract_at_stations(grid, station_coords):
    """Extract TEC at station positions via bilinear interpolation.

    station_coords: {CODE: {'geolat': ..., 'geolon': ...}} or {CODE: (lat, lon)}
    Returns {CODE: tec_value} dict.
    """
    result = {}
    for code, coord in station_coords.items():
        if isinstance(coord, dict):
            lat, lon = coord['geolat'], coord['geolon']
        else:
            lat, lon = coord[0], coord[1]
        result[code] = bilinear_interpolate(grid, lat, lon)
    return result


def grid_values_in_band(grid, lat_min=50, lat_max=75):
    """Extract grid cell values within a latitude band (Northern hemisphere).

    Returns (values_1d, lats_1d, lons_1d, lat_mask) for cells within band.
    """
    lat_mask = (TEC_LATS >= lat_min) & (TEC_LATS <= lat_max)
    sub_grid = grid[lat_mask, :]
    sub_lats = TEC_LATS[lat_mask]

    values = []
    lats = []
    lons = []
    for i, lat in enumerate(sub_lats):
        for j, lon in enumerate(TEC_LONS):
            v = sub_grid[i, j]
            if np.isfinite(v):
                values.append(v)
                lats.append(lat)
                lons.append(lon)

    return np.array(values), np.array(lats), np.array(lons), lat_mask


def available_dates():
    """Scan disk for all available IONEX dates. Returns list of (year, doy)."""
    dates = []
    if not GIM_DIR.exists():
        return dates
    for year_dir in sorted(GIM_DIR.iterdir()):
        if not year_dir.is_dir():
            continue
        try:
            year = int(year_dir.name)
        except ValueError:
            continue
        for f in sorted(year_dir.iterdir()):
            name = f.name
            if name.startswith('CODG') and ('I' in name):
                try:
                    doy = int(name[4:7])
                    dates.append((year, doy))
                except ValueError:
                    continue
    return dates


def check_date_available(year, month, day):
    """Check if IONEX data exists for a specific date."""
    doy = _doy(year, month, day)
    return _find_ionex_file(year, doy) is not None
