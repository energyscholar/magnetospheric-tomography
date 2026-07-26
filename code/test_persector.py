#!/usr/bin/env python3
"""Per-sector decomposition: is the altitude switch a Fennoscandia-only artifact,
or does it appear independently within each densely-sampled sector? Decides the reframe.
Run from scripts/."""
import numpy as np
import loader_supermag as sm
import loader_tec as tec
import relational_diagnostics as rd
import event_catalog as catalog

SECTORS = {  # name: (lon_min, lon_max, lat_min, lat_max) in geographic deg (lon 0-360)
    'Fennoscandia/IMAGE': (0, 45, 55, 82),
    'Greenland':          (300, 345, 58, 85),
    'N.America/Canada':   (200, 300, 45, 82),
}

def in_sector(c, box):
    la = c.get('geolat', np.nan); lo = c.get('geolon', np.nan)
    if not (np.isfinite(la) and np.isfinite(lo)): return False
    lo %= 360
    lomin, lomax, lamin, lamax = box
    return (lomin <= lo <= lomax) and (lamin <= la <= lamax)

def deltas(tab):
    return [(ut, m.get('delta_120_300'), m.get('I_N'), m.get('rms_N'))
            for ut, m in tab if m and m.get('delta_120_300') is not None]

for ename in ['Halloween 2003', 'Bastille Day 2000']:
    ev = [e for e in catalog.EVENTS if ename.lower() in e.name.lower()][0]
    sm_data, sm_coords = sm.load_event_day(ev.year, ev.month, ev.day, min_valid_frac=0.3)
    tec_maps = tec.load_day(ev.year, ev.month, ev.day)
    print(f"\n{'='*70}\n{ev.name}  ({len(sm_data)} stations total)\n{'='*70}")

    # full network reference
    full = deltas(rd.multi_epoch_table(sm_data, sm_coords, tec_maps, '2hr'))
    fd = np.array([d for _, d, _, _ in full])
    print(f"FULL NETWORK: mean Δ={np.nanmean(fd):+.3f}  range=[{np.nanmin(fd):+.3f},{np.nanmax(fd):+.3f}]  "
          f"sign-changes={int(np.sum(np.diff(np.sign(fd))!=0))}  n_epochs={len(fd)}")

    for name, box in SECTORS.items():
        keep = [k for k in sm_data if in_sector(sm_coords.get(k, {}), box)]
        if len(keep) < 6:
            print(f"  {name:22s}: {len(keep)} stations — too few, skipped")
            continue
        smd = {k: sm_data[k] for k in keep}; smc = {k: sm_coords[k] for k in keep}
        d = deltas(rd.multi_epoch_table(smd, smc, tec_maps, '2hr'))
        if not d:
            print(f"  {name:22s}: {len(keep)} stations — no altitude epochs (no auroral/TEC coverage)")
            continue
        dd = np.array([x for _, x, _, _ in d])
        IN = np.array([i for _, _, i, _ in d])
        sc = int(np.sum(np.diff(np.sign(dd)) != 0))
        # does it switch main->recovery like the paper claims? (sign of first third vs last third)
        n = len(dd); early = np.nanmean(dd[:max(1, n//3)]); late = np.nanmean(dd[-max(1, n//3):])
        print(f"  {name:22s}: {len(keep):2d} stn | mean Δ={np.nanmean(dd):+.3f} "
              f"range=[{np.nanmin(dd):+.3f},{np.nanmax(dd):+.3f}] sign-changes={sc} "
              f"| early={early:+.3f} late={late:+.3f} {'SWITCH' if np.sign(early)!=np.sign(late) else '—'} "
              f"| meanI_N={np.nanmean(IN):.3f}")
