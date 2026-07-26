#!/usr/bin/env python3
"""Coordinate-robustness (Gjerloev's remaining data-handling question): does the ground
spatial coherence depend on the N-component projection, or does it hold on the
coordinate-invariant horizontal magnitude |H| = sqrt(N^2+E^2)? Run from scripts/."""
import numpy as np
import loader_supermag as sm
import loader_tec as tec
import relational_diagnostics as rd
import event_catalog as catalog

def IN(tab): return np.array([m['I_N'] for _, m in tab if m and m.get('I_N') is not None])

def to_H(smd):  # replace N-column (0) with horizontal magnitude sqrt(N^2+E^2)
    out = {}
    for k, v in smd.items():
        vv = np.array(v, dtype=float, copy=True)
        vv[:, 0] = np.sqrt(vv[:, 0]**2 + vv[:, 1]**2)
        out[k] = vv
    return out

for ename in ['Halloween 2003', 'Bastille Day 2000']:
    ev = [e for e in catalog.EVENTS if ename.lower() in e.name.lower()][0]
    smd, smc = sm.load_event_day(ev.year, ev.month, ev.day, min_valid_frac=0.3)
    tm = tec.load_day(ev.year, ev.month, ev.day)
    I_N = IN(rd.multi_epoch_table(smd, smc, tm, '2hr'))
    I_H = IN(rd.multi_epoch_table(to_H(smd), smc, tm, '2hr'))
    n = min(len(I_N), len(I_H))
    r = np.corrcoef(I_N[:n], I_H[:n])[0, 1] if n > 2 else np.nan
    print(f"{ev.name:20s}: I_N mean {np.nanmean(I_N):.3f}  I_|H| mean {np.nanmean(I_H):.3f}  "
          f"corr(I_N, I_|H|) = {r:+.3f}  (n={n})")
