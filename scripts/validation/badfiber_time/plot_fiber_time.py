"""
plot_fiber_time.py  <fiber>  <TRACER> [<TRACER> ...]

Plot failure rate vs time (MJD) for one fiber, one panel per tracer.
Loads from pre-built NPZ files in the same directory.

Examples:
  python plot_fiber_time.py 725 BGS LRG ELG QSO
  python plot_fiber_time.py 986 LRG ELG
"""

import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import datetime

if len(sys.argv) < 3:
    sys.exit('Usage: plot_fiber_time.py <fiber> <TRACER> [<TRACER> ...]')

FIBER   = int(sys.argv[1])
TRACERS = [t.upper() for t in sys.argv[2:]]
NBINS   = 32
PETAL   = FIBER // 500

DATADIR = '/global/cfs/cdirs/desi/survey/catalogs/DA3/LSS/matterhorn-v2/redshift_assessment_plots'
SUMFILE = '/global/cfs/cdirs/desi/survey/catalogs/DA3/LSS/matterhorn-v2/redshift_assessment_plots/LRGfibersim_matterhorn-v2.txt'

TRACER_COLOR = {
    'BGS': 'forestgreen',
    'LRG': 'firebrick',
    'ELG': 'steelblue',
    'QSO': 'darkorange',
}

MJD_EPOCH = datetime.date(1858, 11, 17)

def date_to_mjd(yyyymmdd):
    s = str(yyyymmdd)
    d = datetime.date(int(s[:4]), int(s[4:6]), int(s[6:8]))
    return (d - MJD_EPOCH).days

# CCD swaps per petal
CCD_SWAPS_ALL = {
    0: [(date_to_mjd(20260106), 'r', 'sta-dieid164')],
    1: [(date_to_mjd(20220613), 'b', 'sn22805'),
        (date_to_mjd(20230727), 'b', 'STA31230'),
        (date_to_mjd(20231129), 'z', 'M1-27'),
        (date_to_mjd(20241119), 'r', 'sta-dieid159')],
    2: [(date_to_mjd(20201113), 'r', 'M1-28')],
    3: [(date_to_mjd(20250514), 'z', 'M1-02-1')],
    4: [(date_to_mjd(20221108), 'b', 'sn22822'),
        (date_to_mjd(20221108), 'r', 'M1-12')],
    5: [(date_to_mjd(20210622), 'b', 'sn22817'),
        (date_to_mjd(20221108), 'z', 'Mosaic-3')],
    6: [],
    7: [(date_to_mjd(20240724), 'r', 'M1-35'),
        (date_to_mjd(20251007), 'r', 'dieid390'),
        (date_to_mjd(20251007), 'z', 'M1-35')],
    8: [(date_to_mjd(20230507), 'b', 'sn17986')],
    9: [(date_to_mjd(20210615), 'r', 'M1-22-1'),
        (date_to_mjd(20260429), 'r', 'sta-dieid395')],
}
CHAN_COLOR = {'b': 'royalblue', 'r': 'crimson', 'z': 'darkorchid'}

MJD_YEAR = {yr: (datetime.date(yr, 1, 1) - MJD_EPOCH).days for yr in range(2020, 2027)}

# nsig from LRG summary (if available)
nsig_map = {}
try:
    from scipy import special
    data = np.loadtxt(SUMFILE, usecols=range(10))
    for row in data:
        fib  = int(row[0])
        pval = float(row[1])
        nsig_map[fib] = -np.sqrt(2) * special.erfcinv(2 * np.clip(pval, 1e-10, 1-1e-10))
except Exception:
    pass

def load_fiber(tracer, fiber):
    path = f'{DATADIR}/{tracer}_fiber_time.npz'
    d = np.load(path)
    s, e = d['index'][fiber]
    if s < 0:
        return np.array([]), np.array([], dtype=bool)
    return d['mjd'][s:e], d['zsuc'][s:e]

SHARED_BINS = None  # set after loading all panels

def binned_failrate(mjds, sucs, bins):
    bcs, frates, errs = [], [], []
    for b0, b1 in zip(bins[:-1], bins[1:]):
        m = (mjds >= b0) & (mjds < b1)
        n = m.sum()
        bcs.append(0.5 * (b0 + b1))
        if n < 3:
            frates.append(np.nan); errs.append(np.nan)
        else:
            f = (~sucs[m]).sum() / n
            frates.append(f)
            errs.append(np.sqrt(f * (1 - f) / n))
    return np.array(bcs), np.array(frates), np.array(errs)

def draw_panel(ax, tracer, mjds, sucs):
    color = TRACER_COLOR.get(tracer, 'black')
    n_tot  = len(mjds)
    n_fail = (~sucs).sum() if n_tot else 0
    mean_fr = n_fail / n_tot if n_tot else np.nan

    extra = f'  σ={nsig_map[FIBER]:.1f}' if (tracer == 'LRG' and FIBER in nsig_map) else ''
    ax.set_title(f'{tracer}  n={n_tot}  fail={n_fail}{extra}', fontsize=10, pad=3)

    if n_tot < 5:
        ax.text(0.5, 0.5, 'no data', ha='center', va='center',
                transform=ax.transAxes, fontsize=9)
        return

    bcs, fr, er = binned_failrate(mjds, sucs, SHARED_BINS)
    ok = np.isfinite(fr)
    ax.errorbar(bcs[ok], fr[ok], yerr=er[ok],
                fmt='o-', color=color, ms=4, lw=1.5, capsize=3, elinewidth=1,
                label='Fail rate')
    ax.axhline(mean_fr, color='gray', lw=1.0, linestyle='--', alpha=0.7,
               label=f'Mean={mean_fr:.3f}')
    ax.set_ylim(bottom=0)
    ax.set_ylabel('Failure rate', fontsize=10)
    ax.tick_params(labelsize=9)
    ax.legend(fontsize=8, loc='upper left')

    vmin, vmax = mjds.min(), mjds.max()
    for yr, mjd_yr in MJD_YEAR.items():
        if vmin < mjd_yr < vmax:
            ax.axvline(mjd_yr, color='lightgray', lw=0.7, linestyle=':')

    for mjd_swap, ch, label in CCD_SWAPS_ALL.get(PETAL, []):
        if vmin < mjd_swap < vmax:
            col = CHAN_COLOR[ch]
            ax.axvline(mjd_swap, color=col, lw=1.4, linestyle='--', alpha=0.9)
            ax.text(mjd_swap + 5, ax.get_ylim()[1] * 0.5,
                    f'{ch}:{label}', fontsize=7, color=col, rotation=90, va='top')

# ── load data ─────────────────────────────────────────────────────────────────
panels = []
for tracer in TRACERS:
    mjds, sucs = load_fiber(tracer, FIBER)
    panels.append((tracer, mjds, sucs))
    print(f'{tracer}: {len(mjds)} obs, {(~sucs).sum() if len(sucs) else 0} failures', flush=True)

# shared bins from global min/max across all tracers
all_mjds = np.concatenate([m for _, m, _ in panels if len(m)])
SHARED_BINS = np.linspace(all_mjds.min(), all_mjds.max(), NBINS + 1)

# ── figure ────────────────────────────────────────────────────────────────────
n = len(panels)
fig, axes = plt.subplots(n, 1, figsize=(11, 3.5 * n), sharex=True)
if n == 1:
    axes = [axes]

fig.suptitle(f'Fiber {FIBER} (petal {PETAL}) — failure rate vs time  '
             f'[{NBINS} bins, matterhorn-v2]\n'
             'Dashed verticals: CCD swaps  (blue=b, red=r, purple=z)',
             fontsize=12, y=1.01)

for ax, (tracer, mjds, sucs) in zip(axes, panels):
    draw_panel(ax, tracer, mjds, sucs)

axes[-1].set_xlabel('MJD', fontsize=11)

# top axis: monthly ticks, yearly ticks taller with year labels
ax2 = axes[0].twiny()
xmin, xmax = SHARED_BINS[0], SHARED_BINS[-1]
ax2.set_xlim(xmin, xmax)
monthly = [(datetime.date(yr, mo, 1) - MJD_EPOCH).days
           for yr in range(2020, 2027) for mo in range(1, 13)]
yearly  = [(datetime.date(yr,  1, 1) - MJD_EPOCH).days for yr in range(2020, 2027)]
yr_labels = [str(yr) for yr in range(2020, 2027)]
ax2.set_xticks(yearly)
ax2.set_xticklabels(yr_labels, fontsize=8)
ax2.tick_params(axis='x', which='major', direction='in', length=12)
ax2.set_xticks(monthly, minor=True)
ax2.tick_params(axis='x', which='minor', direction='in', length=5)

plt.tight_layout()
tracer_str = '_'.join(TRACERS).lower()
outpng = f'fiber{FIBER}_time_{tracer_str}.png'
fig.savefig(outpng, dpi=150, bbox_inches='tight')
plt.close(fig)
print(f'Saved {outpng}', flush=True)
