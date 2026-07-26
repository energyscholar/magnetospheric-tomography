#!/usr/bin/env python3
"""Figure 1, panel (a): concept schematic of the altitude-projection idea.

Story, left -> right:
  3-D altitude-stratified current systems  ->  collapse onto the ground (I_gnd)
  ->  [+ GPS-TEC, 2-layer Chapman]  ->  un-projected into I_120 / I_300 -> Delta.

Colors match the paper: I_gnd black, I_120 (E-region) blue, I_300 (F-region) red.
Run from scripts/.  Writes ../paper/figures/fig1a_concept.{pdf,png}
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, FancyArrowPatch
from matplotlib.lines import Line2D

rng = np.random.default_rng(3)
C_E, C_F, C_RC, C_GND = '#1f5fd0', '#d02020', '#7a3ea8', '#333333'
DEPTH = (-1.35, 0.95)      # isometric depth offset (dx, dy) for a horizontal plane
W = 4.2                    # plane width

def plane(ax, x0, y0, color, alpha=0.20, edge=None, lw=1.4, z=1):
    dx, dy = DEPTH
    pts = [(x0, y0), (x0 + W, y0), (x0 + W + dx, y0 + dy), (x0 + dx, y0 + dy)]
    ax.add_patch(Polygon(pts, closed=True, facecolor=color, alpha=alpha,
                         edgecolor=edge or color, lw=lw, zorder=z))
    return pts

def texture(ax, x0, y0, color, n, spread, z=3):
    """scatter of coherence 'cells' on a plane; small spread=patchy, large=smooth."""
    dx, dy = DEPTH
    u = rng.uniform(0.08, 0.92, n); v = rng.uniform(0.12, 0.88, n)
    px = x0 + u * W + v * dx
    py = y0 + v * dy
    s = rng.normal(spread, spread * 0.25, n).clip(8, None)
    ax.scatter(px, py, s=s, c=color, alpha=0.55, edgecolors='none', zorder=z)

fig, ax = plt.subplots(figsize=(9.2, 5.6))
ax.set_xlim(-1.2, 15.4); ax.set_ylim(-0.4, 7.4); ax.axis('off')

# ---------- LEFT: the 3-D stack collapsing to the ground ----------
xL = 0.0
y_gnd, y_E, y_F, y_RC = 0.4, 2.3, 3.9, 5.5

plane(ax, xL, y_RC, C_RC, alpha=0.10, lw=1.0)             # ring current (faint)
ax.annotate('', xy=(xL + 3.4, y_RC + 0.55), xytext=(xL + 0.9, y_RC + 0.15),
            arrowprops=dict(arrowstyle='->', color=C_RC, lw=1.6,
                            connectionstyle='arc3,rad=0.45'), zorder=4)
ax.text(xL + W + DEPTH[0] + 0.15, y_RC + DEPTH[1] + 0.15, 'ring current',
        color=C_RC, fontsize=9, style='italic', ha='left', va='bottom')

plane(ax, xL, y_F, C_F)                                    # F-region 300 km
texture(ax, xL, y_F, C_F, 7, 210)                          # smooth / large-scale (EUV)
ax.text(xL - 0.05, y_F + 0.28, 'F-region  ~300 km', color=C_F, fontsize=10,
        fontweight='bold', ha='left', va='bottom')
ax.text(xL - 0.05, y_F + 0.02, 'solar EUV · Pedersen', color=C_F, fontsize=8,
        ha='left', va='top')

plane(ax, xL, y_E, C_E)                                    # E-region 120 km
texture(ax, xL, y_E, C_E, 34, 42)                          # patchy / small-scale (auroral)
ax.text(xL - 0.05, y_E + 0.28, 'E-region  ~120 km', color=C_E, fontsize=10,
        fontweight='bold', ha='left', va='bottom')
ax.text(xL - 0.05, y_E + 0.02, 'auroral precip. · Hall', color=C_E, fontsize=8,
        ha='left', va='top')

plane(ax, xL, y_gnd, C_GND, alpha=0.14, lw=1.4)            # ground
dx, dy = DEPTH
nst = 11
su = np.linspace(0.1, 0.9, nst); sv = rng.uniform(0.25, 0.75, nst)
stx = xL + su * W + sv * dx; sty = y_gnd + sv * dy
ax.scatter(stx, sty, s=26, marker='^', c='k', zorder=5)
ax.text(xL - 0.05, y_gnd - 0.02, 'ground magnetometers', color='k', fontsize=9,
        ha='left', va='top')

# projection: dashed lines from both layers straight down to the ground plane
for frac in (0.30, 0.62):
    tx = xL + frac * W + 0.55 * dx
    ax.plot([tx, tx], [y_gnd + 0.55 * dy, y_F + 0.55 * dy], ls=(0, (2, 2)),
            color='0.45', lw=1.0, zorder=2)
    ax.annotate('', xy=(tx, y_gnd + 0.55 * dy + 0.05), xytext=(tx, y_E + 0.55 * dy),
                arrowprops=dict(arrowstyle='-|>', color='0.35', lw=1.3), zorder=6)

ax.text(xL + 1.6, y_gnd - 0.62,
        r'projection collapses altitude' + '\n' + r'$\Rightarrow\ I_\mathrm{gnd}$ (altitude-blurred)',
        color='k', fontsize=8.5, ha='center', va='top')

# ---------- MIDDLE: the method arrow ----------
xM = 6.9
ax.add_patch(FancyArrowPatch((xM - 0.3, 3.3), (xM + 1.9, 3.3),
             arrowstyle='-|>', mutation_scale=22, lw=2.2, color='#222', zorder=6))
ax.text(xM + 0.8, 3.75, '+ GPS-TEC\n2-layer Chapman', ha='center', va='bottom',
        fontsize=8.8, fontweight='bold')
ax.text(xM + 0.8, 2.98, 'un-project by altitude', ha='center', va='top',
        fontsize=8, style='italic', color='0.3')

# ---------- RIGHT: recovered altitude-resolved coherence + Delta ----------
xR = 9.6
plane(ax, xR, y_F, C_F); texture(ax, xR, y_F, C_F, 7, 210)
plane(ax, xR, y_E, C_E); texture(ax, xR, y_E, C_E, 34, 42)
ax.text(xR + W + dx + 0.2, y_F + dy + 0.05, r'$I_{300}$', color=C_F,
        fontsize=13, fontweight='bold', ha='left', va='center')
ax.text(xR + W + dx + 0.2, y_E + dy + 0.05, r'$I_{120}$', color=C_E,
        fontsize=13, fontweight='bold', ha='left', va='center')
ax.text(xR + 1.9, y_gnd + 0.15,
        r'$\Delta = I_{120} - I_{300}$', color='k', fontsize=12, ha='center', va='center')
ax.text(xR + 1.9, y_gnd - 0.42, 'altitude coherence switch', color='0.3',
        fontsize=8.5, style='italic', ha='center', va='center')

ax.text(xL - 0.9, 7.15, '(a)', fontsize=13, fontweight='bold', ha='left', va='top')

plt.tight_layout(pad=0.4)
out = '../paper/figures/fig1a_concept'
fig.savefig(out + '.pdf', bbox_inches='tight')
fig.savefig(out + '.png', dpi=200, bbox_inches='tight')
print('wrote', out + '.pdf / .png')
