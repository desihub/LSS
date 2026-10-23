"""
plot_petal_mosaic.py
Mosaic of LRG failure-rate vs time for every fiber in one petal.
Layout: N_COLS x N_ROWS panels (default 5 x 100 = all 500 fibers in a petal).
Bad fibers (Krolewski nsig <= -4) are highlighted with a red panel title.
Shows fiber fail rate, survey average (light gray), and model prediction
(dark gray dashed).  CCD swap dates overlaid per petal.

Usage:
  python plot_petal_mosaic.py            # petal 0
  python plot_petal_mosaic.py --petal 7  # petal 7
"""

import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import datetime
from scipy import special

parser = argparse.ArgumentParser()
parser.add_argument('--petal', type=int, default=0)
args = parser.parse_args()
PETAL = args.petal

DATADIR = '/global/u1/r/rohlf/bao_jr/dr3-matterhorn-checks/all-fibers-vs-time'
BADFIB  = '/global/u1/r/rohlf/bao_jr/dr3-matterhorn-checks/krolewski-pipeline/badfibers/LRG_badfibers.txt'
SUMFILE = '/global/u1/r/rohlf/bao_jr/dr3-matterhorn-checks/krolewski-pipeline/summaryscale/LRGfibersim_matterhorn-v2.txt'
NBINS   = 32
N_COLS  = 5
N_ROWS  = 100   # 5 x 100 = 500 fibers per petal

MJD_EPOCH = datetime.date(1858, 11, 17)

def date_to_mjd(yyyymmdd):
    s = str(yyyymmdd)
    d = datetime.date(int(s[:4]), int(s[4:6]), int(s[6:8]))
    return (d - MJD_EPOCH).days

CCD_SWAPS_ALL = {
    0: [(date_to_mjd(20260106), 'r')],
    1: [(date_to_mjd(20220613), 'b'), (date_to_mjd(20230727), 'b'),
        (date_to_mjd(20231129), 'z'), (date_to_mjd(20241119), 'r')],
    2: [(date_to_mjd(20201113), 'r')],
    3: [(date_to_mjd(20250514), 'z')],
    4: [(date_to_mjd(20221108), 'b'), (date_to_mjd(20221108), 'r')],
    5: [(date_to_mjd(20210622), 'b'), (date_to_mjd(20221108), 'z')],
    6: [],
    7: [(date_to_mjd(20240724), 'r'), (date_to_mjd(20251007), 'r'),
        (date_to_mjd(20251007), 'z')],
    8: [(date_to_mjd(20230507), 'b')],
    9: [(date_to_mjd(20210615), 'r'), (date_to_mjd(20260429), 'r')],
}
CHAN_COLOR = {'b': 'royalblue', 'r': 'crimson', 'z': 'darkorchid'}
MJD_YEAR = {yr: (datetime.date(yr, 1, 1) - MJD_EPOCH).days for yr in range(2020, 2027)}

# ── bad fiber set and nsig map ─────────────────────────────────────────────────
bad_fibers = set(np.loadtxt(BADFIB, dtype=int).tolist())
nsig_map = {}
try:
    data = np.loadtxt(SUMFILE, usecols=range(10))
    for row in data:
        fib = int(row[0])
        pval = float(row[1])
        nsig_map[fib] = -np.sqrt(2) * special.erfcinv(2 * np.clip(pval, 1e-10, 1-1e-10))
except Exception as e:
    print(f'nsig load failed: {e}')

# ── load NPZ ───────────────────────────────────────────────────────────────────
print('Loading LRG NPZ...', flush=True)
d = np.load(f'{DATADIR}/LRG_fiber_time.npz')
all_mjd  = d['mjd'].astype(np.float64)
all_fail = (~d['zsuc']).astype(np.float64)
all_mod  = d['mod_success_rate'].astype(np.float64)

bins = np.linspace(all_mjd.min(), all_mjd.max(), NBINS + 1)
bcs  = 0.5 * (bins[:-1] + bins[1:])

# survey average across all fibers
n_tot_sv  = np.histogram(all_mjd, bins=bins)[0].astype(np.float64)
n_fail_sv = np.histogram(all_mjd, bins=bins, weights=all_fail)[0]
with np.errstate(invalid='ignore', divide='ignore'):
    sv_fr = np.where(n_tot_sv >= 10, n_fail_sv / n_tot_sv, np.nan)
sv_er = np.where(np.isfinite(sv_fr),
                 np.sqrt(sv_fr * (1 - sv_fr) / np.maximum(n_tot_sv, 1)), np.nan)

scale = (1 - all_fail.mean()) / all_mod.mean()
print(f'  scale={scale:.4f}', flush=True)

def get_fiber(fib):
    s, e = d['index'][fib]
    if s < 0:
        return np.array([]), np.array([], dtype=bool), np.array([])
    return d['mjd'][s:e], d['zsuc'][s:e], d['mod_success_rate'][s:e]

def binned_failrate(mjds, sucs):
    fr, er = [], []
    for b0, b1 in zip(bins[:-1], bins[1:]):
        m = (mjds >= b0) & (mjds < b1)
        n = m.sum()
        if n < 3:
            fr.append(np.nan); er.append(np.nan)
        else:
            f = (~sucs[m]).sum() / n
            fr.append(f)
            er.append(np.sqrt(f * (1 - f) / n))
    return np.array(fr), np.array(er)

def model_pred(mjds, mods):
    pred = np.full(len(bcs), np.nan)
    for i, (b0, b1) in enumerate(zip(bins[:-1], bins[1:])):
        m = (mjds >= b0) & (mjds < b1)
        if m.sum() >= 3:
            pred[i] = 1 - np.clip(mods[m] * scale, 0, 1).mean()
    return pred

# ── fibers in this petal ───────────────────────────────────────────────────────
fibers = list(range(PETAL * 500, (PETAL + 1) * 500))
print(f'Petal {PETAL}: fibers {fibers[0]}–{fibers[-1]}  ({len(fibers)} total)', flush=True)

# ── figure ────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(N_ROWS, N_COLS,
                         figsize=(N_COLS * 2.5, N_ROWS * 2.0),
                         squeeze=False, sharex=True)
fig.suptitle(f'LRG petal {PETAL} — all {len(fibers)} fibers — failure rate vs time  (matterhorn-v2)\n'
             'Red: fiber  |  Light gray: survey avg  |  Dark dashed: model prediction  |  '
             'Red title = Krolewski bad fiber',
             fontsize=10, y=1.002)

for fi, fib in enumerate(fibers):
    ax    = axes[fi // N_COLS, fi % N_COLS]
    mjds, sucs, mods = get_fiber(fib)

    ns     = nsig_map.get(fib, np.nan)
    ns_str = f' σ={ns:.1f}' if np.isfinite(ns) else ''
    is_bad = fib in bad_fibers
    title_color = 'darkorange' if is_bad else 'black'
    ax.set_title(f'fib {fib}{ns_str}', fontsize=5, pad=2, color=title_color,
                 fontweight='bold' if is_bad else 'normal')

    if len(mjds) < 5:
        ax.text(0.5, 0.5, 'no data', ha='center', va='center',
                transform=ax.transAxes, fontsize=5, color='gray')
        ax.tick_params(labelsize=4)
        continue

    fr, er = binned_failrate(mjds, sucs)
    pred   = model_pred(mjds, mods)
    ok     = np.isfinite(fr)
    sv_ok  = np.isfinite(sv_fr) & ok
    pr_ok  = np.isfinite(pred) & ok

    ax.fill_between(bcs[sv_ok],
                    (sv_fr - sv_er)[sv_ok], (sv_fr + sv_er)[sv_ok],
                    color='black', alpha=0.10, zorder=0)
    ax.plot(bcs[sv_ok], sv_fr[sv_ok], color='black', lw=0.5, alpha=0.25)
    ax.plot(bcs[pr_ok], pred[pr_ok], color='black', lw=0.8,
            linestyle='--', alpha=0.55)
    ax.errorbar(bcs[ok], fr[ok], yerr=er[ok],
                fmt='o-', color='firebrick', ms=1.8, lw=0.9,
                capsize=1.0, elinewidth=0.5)

    ax.set_ylim(bottom=0)
    ax.tick_params(labelsize=4)
    ax.set_ylabel('fail', fontsize=4)

    vmin, vmax = mjds.min(), mjds.max()
    for yr, mjd_yr in MJD_YEAR.items():
        if vmin < mjd_yr < vmax:
            ax.axvline(mjd_yr, color='lightgray', lw=0.3, linestyle=':')
    for mjd_swap, ch in CCD_SWAPS_ALL.get(PETAL, []):
        if vmin < mjd_swap < vmax:
            ax.axvline(mjd_swap, color=CHAN_COLOR[ch], lw=0.7,
                       linestyle='--', alpha=0.8)

# top-row year labels
monthly = [(datetime.date(yr, mo, 1) - MJD_EPOCH).days
           for yr in range(2020, 2027) for mo in range(1, 13)]
yearly  = [(datetime.date(yr,  1, 1) - MJD_EPOCH).days for yr in range(2020, 2027)]
for col in range(N_COLS):
    ax2 = axes[0, col].twiny()
    ax2.set_xlim(bins[0], bins[-1])
    ax2.set_xticks(yearly)
    ax2.set_xticklabels([str(y) for y in range(2020, 2027)], fontsize=4)
    ax2.tick_params(axis='x', which='major', direction='in', length=6)
    ax2.set_xticks(monthly, minor=True)
    ax2.tick_params(axis='x', which='minor', direction='in', length=2)

plt.tight_layout()
outpng = f'{DATADIR}/lrg_petal{PETAL}_mosaic.png'
outpdf = f'{DATADIR}/lrg_petal{PETAL}_mosaic.pdf'
fig.savefig(outpng, dpi=100, bbox_inches='tight')
fig.savefig(outpdf, bbox_inches='tight')
plt.close(fig)
print(f'Saved {outpng}', flush=True)
print(f'Saved {outpdf}', flush=True)
