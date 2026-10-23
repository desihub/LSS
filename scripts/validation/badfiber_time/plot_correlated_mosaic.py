"""
plot_correlated_mosaic.py
Mosaic of LRG+ELG failure-rate vs time for the 36 unique fibers found by
find_correlated_bad_bins.py.  Each panel shows both tracers normalized to
their per-fiber mean (fr / mean_fr) so the two very different scales
(LRG ~3%, ELG ~33%) are comparable.  Sorted by peak z_prod.
Output: correlated_mosaic_lrg_elg.png
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import datetime

DATADIR = '/global/cfs/cdirs/desi/survey/catalogs/DA3/LSS/matterhorn-v2/redshift_assessment_plots'
NBINS   = 32

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

# ── load all 41 result rows, sort by fiber number ────────────────────────────
rows = []
with open(f'{DATADIR}/correlated_bad_bins_lrg_elg.txt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        cols = line.split()
        rows.append({
            'fiber':   int(cols[0]),
            'petal':   int(cols[1]),
            'bin_mjd': float(cols[2]),
            'z_lrg':   float(cols[3]),
            'z_elg':   float(cols[4]),
            'z_prod':  float(cols[5]),
        })

rows.sort(key=lambda r: r['fiber'])
print(f'{len(rows)} result rows', flush=True)

# ── load NPZ data ─────────────────────────────────────────────────────────────
print('Loading data...', flush=True)
lrg = np.load(f'{DATADIR}/LRG_fiber_time.npz')
elg = np.load(f'{DATADIR}/ELG_fiber_time.npz')

def get_fiber(d, fib):
    s, e = d['index'][fib]
    if s < 0:
        return np.array([]), np.array([], dtype=bool)
    return d['mjd'][s:e], d['zsuc'][s:e]

# shared bins from global range
all_lrg_mjd = lrg['mjd']; all_elg_mjd = elg['mjd']
mjd_min = min(all_lrg_mjd.min(), all_elg_mjd.min())
mjd_max = max(all_lrg_mjd.max(), all_elg_mjd.max())
bins = np.linspace(mjd_min, mjd_max, NBINS + 1)
bcs  = 0.5 * (bins[:-1] + bins[1:])

def binned_failrate(mjds, sucs):
    frates, errs = [], []
    for b0, b1 in zip(bins[:-1], bins[1:]):
        m = (mjds >= b0) & (mjds < b1)
        n = m.sum()
        if n < 3:
            frates.append(np.nan); errs.append(np.nan)
        else:
            f = (~sucs[m]).sum() / n
            frates.append(f)
            errs.append(np.sqrt(f * (1 - f) / n))
    return np.array(frates), np.array(errs)

# ── figure layout ─────────────────────────────────────────────────────────────
N = len(rows)
NCOLS = 6
NROWS = int(np.ceil(N / NCOLS))
bin_width = bins[1] - bins[0]

fig, axes = plt.subplots(NROWS, NCOLS, figsize=(NCOLS * 3.5, NROWS * 2.8),
                         squeeze=False, sharex=True)
fig.suptitle('Correlated bad bins: LRG (red) & ELG (blue) normalized fail rate vs time\n'
             'fr / mean_fr  |  matterhorn-v2  |  ordered by fiber number  |  flagged bin shaded',
             fontsize=12, y=1.01)

for fi, row in enumerate(rows):
    ax    = axes[fi // NCOLS, fi % NCOLS]
    fib   = row['fiber']
    petal = row['petal']

    lrg_mjd, lrg_suc = get_fiber(lrg, fib)
    elg_mjd, elg_suc = get_fiber(elg, fib)

    lrg_fr, lrg_er = binned_failrate(lrg_mjd, lrg_suc)
    elg_fr, elg_er = binned_failrate(elg_mjd, elg_suc)

    lrg_mean = (~lrg_suc).mean() if len(lrg_suc) else np.nan
    elg_mean = (~elg_suc).mean() if len(elg_suc) else np.nan

    lrg_norm  = lrg_fr / lrg_mean
    elg_norm  = elg_fr / elg_mean
    lrg_enorm = lrg_er / lrg_mean
    elg_enorm = elg_er / elg_mean

    ok_l = np.isfinite(lrg_norm)
    ok_e = np.isfinite(elg_norm)

    ax.errorbar(bcs[ok_l], lrg_norm[ok_l], yerr=lrg_enorm[ok_l],
                fmt='o-', color='firebrick', ms=2, lw=1.0, capsize=1.5,
                elinewidth=0.6, label='LRG')
    ax.errorbar(bcs[ok_e], elg_norm[ok_e], yerr=elg_enorm[ok_e],
                fmt='s-', color='steelblue', ms=2, lw=1.0, capsize=1.5,
                elinewidth=0.6, label='ELG')
    ax.axhline(1.0, color='gray', lw=0.6, linestyle='--', alpha=0.6)

    # highlight flagged bin
    bm = row['bin_mjd']
    ax.axvspan(bm - bin_width/2, bm + bin_width/2, color='gold', alpha=0.3, zorder=0)

    ax.set_ylim(bottom=0)
    ax.tick_params(labelsize=5)

    lrg_n = len(lrg_mjd); elg_n = len(elg_mjd)
    ax.set_title(f'fib {fib} p{petal}  L:{lrg_n} E:{elg_n}\n'
                 f'zL={row["z_lrg"]:.1f} zE={row["z_elg"]:.1f}', fontsize=6, pad=2)

    vmin = min(lrg_mjd.min() if len(lrg_mjd) else mjd_max,
               elg_mjd.min() if len(elg_mjd) else mjd_max)
    vmax = max(lrg_mjd.max() if len(lrg_mjd) else mjd_min,
               elg_mjd.max() if len(elg_mjd) else mjd_min)

    for yr, mjd_yr in MJD_YEAR.items():
        if vmin < mjd_yr < vmax:
            ax.axvline(mjd_yr, color='lightgray', lw=0.4, linestyle=':')

    for mjd_swap, ch in CCD_SWAPS_ALL.get(petal, []):
        if vmin < mjd_swap < vmax:
            ax.axvline(mjd_swap, color=CHAN_COLOR[ch], lw=0.9,
                       linestyle='--', alpha=0.8)

    ax.set_ylabel('fr/mean', fontsize=5)
    ax.set_xlabel('MJD', fontsize=5)

axes[0, 0].legend(fontsize=5, loc='upper left')

for fi in range(N, NROWS * NCOLS):
    axes[fi // NCOLS, fi % NCOLS].set_visible(False)

# top axis on every panel: year labels + monthly ticks
monthly = [(datetime.date(yr, mo, 1) - MJD_EPOCH).days
           for yr in range(2020, 2027) for mo in range(1, 13)]
yearly  = [(datetime.date(yr,  1, 1) - MJD_EPOCH).days for yr in range(2020, 2027)]
yr_labels = [str(yr) for yr in range(2020, 2027)]

for fi in range(N):
    ax2 = axes[fi // NCOLS, fi % NCOLS].twiny()
    ax2.set_xlim(mjd_min, mjd_max)
    ax2.set_xticks(yearly)
    ax2.set_xticklabels(yr_labels, fontsize=5)
    ax2.tick_params(axis='x', which='major', direction='in', length=8)
    ax2.set_xticks(monthly, minor=True)
    ax2.tick_params(axis='x', which='minor', direction='in', length=3)

plt.tight_layout()
outpng = f'{DATADIR}/correlated_mosaic_lrg_elg.png'
fig.savefig(outpng, dpi=150, bbox_inches='tight')
plt.close(fig)
print(f'Saved {outpng}', flush=True)
