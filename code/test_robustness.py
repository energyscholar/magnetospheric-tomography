#!/usr/bin/env python3
"""Gjerloev-prep robustness tests for the altitude-coherence paper:
(1) baseline-perturbation stability of I_gnd, (2) IMAGE-drop jackknife.
Validates against cached values first, then perturbs. Run from scripts/."""
import json, numpy as np
from pathlib import Path
import loader_supermag as sm
import loader_tec as tec
import relational_diagnostics as rd
import event_catalog as catalog

rng = np.random.default_rng(0)
BASE = Path(__file__).parent.parent
CACHE = json.load(open(BASE / 'results/paper1_data/halloween_2003_altitude_data.json'))

ev = [e for e in catalog.EVENTS if 'halloween' in e.name.lower()][0]
sm_data, sm_coords = sm.load_event_day(ev.year, ev.month, ev.day, min_valid_frac=0.3)
tec_maps = tec.load_day(ev.year, ev.month, ev.day)
print(f"loaded: {len(sm_data)} stations, {len(tec_maps) if tec_maps else 0} tec epochs")

k0 = list(sm_data.keys())[0]
print(f"sm_data[{k0}] type={type(sm_data[k0]).__name__} shape={getattr(sm_data[k0],'shape',None)}")
print(f"sm_coords[{k0}] = {sm_coords[k0]}")

def run(smd, smc):
    return rd.multi_epoch_table(smd, smc, tec_maps, '2hr')

def IN(tab):    return np.array([m['I_N'] for _, m in tab if m and m.get('I_N') is not None])
def DELTA(tab): return np.array([m.get('delta_120_300') for _, m in tab
                                 if m and m.get('delta_120_300') is not None])

# ---- VALIDATE reproduction ----
tab = run(sm_data, sm_coords)
m0 = tab[0][1]
print("\n== VALIDATION (reproduce cache, epoch 0) ==")
print(f"  repro: I_N={m0.get('I_N'):.4f}  I_120={m0['altitude_I'].get(120):.4f}  I_300={m0['altitude_I'].get(300):.4f}")
c0 = CACHE['epochs'][0]
print(f"  cache: I_N={c0['I_N']:.4f}  I_120={c0['I_120']:.4f}  I_300={c0['I_300']:.4f}")
ok = abs(m0.get('I_N') - c0['I_N']) < 1e-3
print(f"  MATCH: {ok}")

base_IN, base_D = IN(tab), DELTA(tab)
print(f"  base I_N mean={np.nanmean(base_IN):.4f}, delta mean={np.nanmean(base_D):+.4f}, n_epochs={len(base_D)}")

if not ok:
    print("\n!! harness does not reproduce cache — NOT running tests (would be untrustworthy)")
    raise SystemExit

# ---- TEST 1: baseline perturbation ----
def perturb(smd, sigma):
    out = {}
    for k, v in smd.items():
        vv = np.array(v, dtype=float, copy=True)
        vv[:, 0] = vv[:, 0] + rng.normal(0, sigma)   # per-station constant N-baseline offset
        out[k] = vv
    return out

print("\n== TEST 1: baseline-perturbation stability of I_gnd ==")
for sigma in [2.0, 5.0]:
    means, sign_pres, M = [], 0, 100
    for _ in range(M):
        pt = run(perturb(sm_data, sigma), sm_coords)
        means.append(np.nanmean(IN(pt)))
        d = DELTA(pt)
        if len(d) == len(base_D) and np.mean(np.sign(d) == np.sign(base_D)) >= 0.9:
            sign_pres += 1
    print(f"  sigma={sigma:.0f} nT: I_N mean {np.mean(means):.4f} +/- {np.std(means):.4f} "
          f"(base {np.nanmean(base_IN):.4f}); delta-sign preserved {sign_pres}/{M}")

# ---- TEST 2: IMAGE-drop jackknife ----
def is_image(c):
    la = c.get('geolat', np.nan); lo = c.get('geolon', np.nan)
    if not (np.isfinite(la) and np.isfinite(lo)): return False
    lo = lo % 360
    return (55 <= la <= 80) and (0 <= lo <= 40)

keep = [k for k in sm_data if not is_image(sm_coords.get(k, {}))]
dropped = len(sm_data) - len(keep)
print(f"\n== TEST 2: IMAGE-drop jackknife (dropped {dropped}/{len(sm_data)} Fennoscandian stations) ==")
tab2 = run({k: sm_data[k] for k in keep}, {k: sm_coords[k] for k in keep})
print(f"  full   : I_N mean {np.nanmean(base_IN):.4f}, delta mean {np.nanmean(base_D):+.4f}")
print(f"  no-IMAGE: I_N mean {np.nanmean(IN(tab2)):.4f}, delta mean {np.nanmean(DELTA(tab2)):+.4f}, "
      f"n_stations~{tab2[0][1].get('n_stations') if tab2 and tab2[0][1] else 'NA'}")
