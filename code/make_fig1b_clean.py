#!/usr/bin/env python3
"""Figure 1, panel (b), CLEAN native redraw (white bg, no tutorial UI chrome).
Real two-layer coherence from the Halloween 2003 cache, representative storm-phase
epoch: I_120=+0.25, I_300=+0.21, Delta=+0.04. Colors match paper (I_120 blue, I_300 red).
Run from scripts/. Writes ../paper/figures/fig1b_clean.{pdf,png}"""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# pull the real representative epoch from the cache (first valid altitude epoch)
d = json.load(open('../data/derived/paper1_data/halloween_2003_altitude_data.json'))
ep = next(e for e in d['epochs'] if e.get('I_120') is not None and e.get('I_300') is not None)
I120, I300 = ep['I_120'], ep['I_300']
delta = ep.get('delta_120_300', I120 - I300)
C_E, C_F = '#1f5fd0', '#d02020'

fig, (axbar, axps) = plt.subplots(1, 2, figsize=(5.4, 3.3),
                                  gridspec_kw=dict(width_ratios=[1.35, 1.0], wspace=0.55))

# --- left: horizontal bars at two altitudes (300 km on top: altitude increases upward) ---
axbar.barh([1], [I300], color=C_F, height=0.5, alpha=0.85, zorder=3)   # 300 km, upper
axbar.barh([0], [I120], color=C_E, height=0.5, alpha=0.85, zorder=3)   # 120 km, lower
axbar.set_yticks([1, 0]); axbar.set_yticklabels(['300 km', '120 km'], fontsize=10)
axbar.text(I300 + 0.012, 1, f'$I_{{300}}={I300:+.2f}$', va='center', ha='left',
           color=C_F, fontsize=10, fontweight='bold')
axbar.text(I120 + 0.012, 0, f'$I_{{120}}={I120:+.2f}$', va='center', ha='left',
           color=C_E, fontsize=10, fontweight='bold')
axbar.axvline(0, color='0.6', lw=0.8, zorder=1)
axbar.set_xlim(0, max(I120, I300) * 1.6); axbar.set_ylim(-0.7, 1.7)
axbar.set_xlabel("Moran's $I$", fontsize=9.5)
axbar.spines[['top', 'right']].set_visible(False)

# --- right: phase space I_120 vs I_300 ---
lim = 0.4
axps.plot([0, lim], [0, lim], ls=(0, (4, 3)), color='0.6', lw=1.0, zorder=1)
axps.scatter([I120], [I300], s=90, color='0.15', zorder=4)
axps.annotate('this epoch', (I120, I300), textcoords='offset points', xytext=(6, -12),
              fontsize=8, color='0.3')
axps.text(0.30, 0.06, 'E-region\nleads ($\\Delta>0$)', fontsize=8, color=C_E,
          ha='center', va='center')
axps.text(0.08, 0.33, 'F-region\nleads ($\\Delta<0$)', fontsize=8, color=C_F,
          ha='center', va='center')
axps.set_xlim(0, lim); axps.set_ylim(0, lim)
axps.set_xlabel(r'$I_{120}$', fontsize=10, color=C_E)
axps.set_ylabel(r'$I_{300}$', fontsize=10, color=C_F)
axps.set_xticks([0, 0.2, 0.4]); axps.set_yticks([0, 0.2, 0.4])
axps.tick_params(labelsize=8)
axps.set_aspect('equal'); axps.spines[['top', 'right']].set_visible(False)

fig.suptitle('real two-layer coherence  ·  Halloween 2003, storm phase',
             fontsize=9.5, y=1.04, color='0.2')
fig.text(0.5, -0.04, r'$\Delta = I_{120}-I_{300}=' + f'{delta:+.2f}$',
         fontsize=11.5, fontweight='bold', ha='center', va='top')
fig.text(0.005, 1.04, '(b)', fontsize=13, fontweight='bold', ha='left', va='top')
out = '../paper/figures/fig1b_clean'
fig.savefig(out + '.pdf', bbox_inches='tight')
fig.savefig(out + '.png', dpi=200, bbox_inches='tight')
print(f'wrote {out}  (I120={I120:+.3f} I300={I300:+.3f} delta={delta:+.3f})')
