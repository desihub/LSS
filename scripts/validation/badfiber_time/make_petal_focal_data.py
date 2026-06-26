"""
make_petal_focal_data.py
Extract all data needed for the interactive focal-plane viewer for one petal.
Outputs: lrg_petal{N}_focal_data.json

JSON structure:
{
  "petal": 0,
  "bins_mjd": [32 bin-centre MJDs],
  "survey_fr": [32 floats],   -- survey-average fail rate
  "survey_er": [32 floats],
  "fibers": {
    "0": { "x": float, "y": float, "is_bad": bool,
           "nsig": float|null, "n_obs": int,
           "fr": [32], "er": [32], "pred": [32] },
    ...
  }
}

Usage:
  python make_petal_focal_data.py             # petal 0
  python make_petal_focal_data.py --petal 7
"""

import argparse, json, os
import numpy as np
import fitsio
from scipy import special

parser = argparse.ArgumentParser()
parser.add_argument('--petal', type=int, default=0)
args = parser.parse_args()
PETAL = args.petal

DATADIR = '/global/u1/r/rohlf/bao_jr/dr3-matterhorn-checks/all-fibers-vs-time'
LRGCAT  = '/global/cfs/cdirs/desi/survey/catalogs/DA3/LSS/matterhorn-v2/LSScats/test/LRG_full_noveto.dat.fits'
BADFIB  = '/global/u1/r/rohlf/bao_jr/dr3-matterhorn-checks/krolewski-pipeline/badfibers/LRG_badfibers.txt'
SUMFILE = '/global/u1/r/rohlf/bao_jr/dr3-matterhorn-checks/krolewski-pipeline/summaryscale/LRGfibersim_matterhorn-v2.txt'
NBINS   = 32

bad_fibers = set(np.loadtxt(BADFIB, dtype=int).tolist())

nsig_map = {}
try:
    data = np.loadtxt(SUMFILE, usecols=range(10))
    for row in data:
        fib = int(row[0])
        pval = float(row[1])
        nsig_map[fib] = float(-np.sqrt(2) * special.erfcinv(2 * np.clip(pval, 1e-10, 1-1e-10)))
except Exception as e:
    print(f'nsig load failed: {e}')

# ── focal plane positions (median per fiber from LRG catalog) ─────────────────
print('Reading focal plane positions...', flush=True)
cat_fp = fitsio.read(LRGCAT, columns=['FIBER', 'FIBERASSIGN_X', 'FIBERASSIGN_Y'])
fp_fiber = np.array(cat_fp['FIBER'],         dtype=np.int32)
fp_x     = np.array(cat_fp['FIBERASSIGN_X'], dtype=np.float32)
fp_y     = np.array(cat_fp['FIBERASSIGN_Y'], dtype=np.float32)
valid    = (fp_fiber >= 0) & (fp_fiber < 5000) & (fp_x != 999999)

fiber_x = np.full(5000, np.nan)
fiber_y = np.full(5000, np.nan)
fib_lo, fib_hi = PETAL * 500, (PETAL + 1) * 500
for fib in range(5000):
    m = valid & (fp_fiber == fib)
    if m.any():
        fiber_x[fib] = float(np.median(fp_x[m]))
        fiber_y[fib] = float(np.median(fp_y[m]))
print(f'  {np.isfinite(fiber_x).sum()} / 5000 total fibers with positions', flush=True)
print(f'  {np.isfinite(fiber_x[fib_lo:fib_hi]).sum()} / 500 petal {PETAL} fibers with positions', flush=True)

# ── load NPZ ──────────────────────────────────────────────────────────────────
print('Loading LRG NPZ...', flush=True)
d = np.load(f'{DATADIR}/LRG_fiber_time.npz')
all_mjd  = d['mjd'].astype(np.float64)
all_fail = (~d['zsuc']).astype(np.float64)
all_mod  = d['mod_success_rate'].astype(np.float64)

bins = np.linspace(all_mjd.min(), all_mjd.max(), NBINS + 1)
bcs  = 0.5 * (bins[:-1] + bins[1:])

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

def binned(mjds, sucs, mods):
    fr, er, pred = [], [], []
    for b0, b1 in zip(bins[:-1], bins[1:]):
        m = (mjds >= b0) & (mjds < b1)
        n = m.sum()
        if n < 3:
            fr.append(None); er.append(None); pred.append(None)
        else:
            f = float((~sucs[m]).sum() / n)
            fr.append(round(f, 6))
            er.append(round(float(np.sqrt(f * (1 - f) / n)), 6))
            pred.append(round(float(1 - np.clip(mods[m] * scale, 0, 1).mean()), 6))
    return fr, er, pred

# ── build output dict ─────────────────────────────────────────────────────────
all_pos = {}
for fib in range(5000):
    if np.isfinite(fiber_x[fib]):
        all_pos[str(fib)] = [round(fiber_x[fib], 1), round(fiber_y[fib], 1)]

out = {
    'petal': PETAL,
    'all_positions': all_pos,
    'bins_mjd': [round(float(v), 2) for v in bcs],
    'survey_fr': [round(float(v), 6) if np.isfinite(v) else None for v in sv_fr],
    'survey_er': [round(float(v), 6) if np.isfinite(v) else None for v in sv_er],
    'fibers': {}
}

print(f'Computing binned rates for petal {PETAL} fibers...', flush=True)
for fib in range(fib_lo, fib_hi):
    mjds, sucs, mods = get_fiber(fib)
    fr, er, pred = binned(mjds, sucs, mods) if len(mjds) >= 5 else ([None]*NBINS, [None]*NBINS, [None]*NBINS)
    out['fibers'][str(fib)] = {
        'x':      round(fiber_x[fib], 2) if np.isfinite(fiber_x[fib]) else None,
        'y':      round(fiber_y[fib], 2) if np.isfinite(fiber_y[fib]) else None,
        'is_bad': fib in bad_fibers,
        'nsig':   round(nsig_map[fib], 2) if fib in nsig_map else None,
        'n_obs':  int(len(mjds)),
        'fr':     fr,
        'er':     er,
        'pred':   pred,
    }

outpath = f'{DATADIR}/lrg_petal{PETAL}_focal_data.json'
with open(outpath, 'w') as f:
    json.dump(out, f, separators=(',', ':'))
print(f'Saved {outpath}  ({os.path.getsize(outpath)//1024} KB)', flush=True)
