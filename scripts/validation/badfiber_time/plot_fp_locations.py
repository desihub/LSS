"""
plot_fp_locations.py
Plot DESI focal plane with all fibers shown as gray dots and the 41
correlated-bad-bin fibers highlighted in color by petal.
"""

import numpy as np
import fitsio
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DATADIR = '/global/u1/r/rohlf/bao_jr/dr3-matterhorn-checks/all-fibers-vs-time'
LRGCAT  = '/global/cfs/cdirs/desi/survey/catalogs/DA3/LSS/matterhorn-v2/LSScats/test/LRG_full_noveto.dat.fits'

# load the 41 flagged fibers
flagged = set()
with open(f'{DATADIR}/correlated_bad_bins_lrg_elg.txt') as f:
    for line in f:
        if not line.startswith('#'):
            flagged.add(int(line.split()[0]))
print(f'{len(flagged)} unique flagged fibers', flush=True)

# read focal plane positions from LRG catalog
print('Reading focal plane positions...', flush=True)
cols = ['FIBER', 'FIBERASSIGN_X', 'FIBERASSIGN_Y']
cat = fitsio.read(LRGCAT, columns=cols)

fiber_arr = np.array(cat['FIBER'],         dtype=np.int32)
x_arr     = np.array(cat['FIBERASSIGN_X'], dtype=np.float32)
y_arr     = np.array(cat['FIBERASSIGN_Y'], dtype=np.float32)

# keep only valid (observed) entries
valid = (fiber_arr >= 0) & (fiber_arr < 5000) & (x_arr != 999999)

# median X, Y per fiber
fiber_x = np.full(5000, np.nan)
fiber_y = np.full(5000, np.nan)
for fib in range(5000):
    m = valid & (fiber_arr == fib)
    if m.any():
        fiber_x[fib] = np.median(x_arr[m])
        fiber_y[fib] = np.median(y_arr[m])

print(f'  {np.isfinite(fiber_x).sum()} fibers with positions', flush=True)

# petal color map
PETAL_COLORS = plt.cm.tab10(np.linspace(0, 1, 10))

# ── plot ──────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 8))
ax.set_aspect('equal')
ax.set_facecolor('#111111')
fig.patch.set_facecolor('#111111')

# all fibers — gray
all_x = fiber_x[np.isfinite(fiber_x)]
all_y = fiber_y[np.isfinite(fiber_x)]
ax.scatter(fiber_x, fiber_y, s=1, c='#555555', zorder=1, rasterized=True)

# flagged fibers — colored by petal, larger
for fib in sorted(flagged):
    if np.isfinite(fiber_x[fib]):
        petal = fib // 500
        ax.scatter(fiber_x[fib], fiber_y[fib],
                   s=60, color=PETAL_COLORS[petal], zorder=3,
                   edgecolors='white', linewidths=0.4)
        ax.text(fiber_x[fib] + 4, fiber_y[fib] + 4,
                str(fib), fontsize=4, color='white', zorder=4)

# focal plane boundary circle (~410 mm radius)
theta = np.linspace(0, 2*np.pi, 300)
ax.plot(410*np.cos(theta), 410*np.sin(theta), 'w-', lw=0.5, alpha=0.3)

# petal boundary lines (10 petals, 36 deg each)
for p in range(10):
    angle = np.radians(p * 36 + 18)   # midpoint of each petal sector
    ax.plot([0, 420*np.cos(angle)], [0, 420*np.sin(angle)],
            color='#333333', lw=0.5)

# legend for petals present
petals_present = sorted(set(fib // 500 for fib in flagged))
handles = [plt.scatter([], [], s=40, color=PETAL_COLORS[p],
                       edgecolors='white', linewidths=0.4,
                       label=f'Petal {p}') for p in petals_present]
ax.legend(handles=handles, fontsize=8, loc='lower right',
          facecolor='#222222', labelcolor='white', framealpha=0.8)

ax.set_xlabel('Focal plane X (mm)', fontsize=10, color='white')
ax.set_ylabel('Focal plane Y (mm)', fontsize=10, color='white')
ax.set_title(f'DESI focal plane — {len(flagged)} fibers with correlated LRG+ELG bad bins\n'
             'matterhorn-v2  |  colored by petal  |  fiber number labeled',
             fontsize=10, color='white')
ax.tick_params(colors='white', labelsize=8)
for spine in ax.spines.values():
    spine.set_edgecolor('#444444')

plt.tight_layout()
outpng = f'{DATADIR}/fp_correlated_fibers.png'
fig.savefig(outpng, dpi=150, bbox_inches='tight', facecolor=fig.get_facecolor())
plt.close(fig)
print(f'Saved {outpng}', flush=True)
