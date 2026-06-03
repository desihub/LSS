"""
find_correlated_bad_bins.py

Find fibers where LRG and ELG failure rates are both >3 sigma above their
respective mean failure rates in the same time bin.

Z-score per bin (binomial test vs overall mean):
    z = (k - n * mean_fr) / sqrt(n * mean_fr * (1 - mean_fr))

Fully vectorized with np.bincount — runs in seconds.
"""

import numpy as np
import time

DATADIR = '/global/u1/r/rohlf/bao_jr/dr3-matterhorn-checks/all-fibers-vs-time'
NBINS   = 32
ZSIG    = 3.0
MIN_OBS = 10   # minimum total obs per fiber per tracer

t0 = time.time()
print('Loading data...', flush=True)
lrg = np.load(f'{DATADIR}/LRG_fiber_time.npz')
elg = np.load(f'{DATADIR}/ELG_fiber_time.npz')

lrg_fib  = lrg['fiber'].astype(np.int32)
lrg_mjd  = lrg['mjd'].astype(np.float64)
lrg_fail = (~lrg['zsuc']).astype(np.float64)

elg_fib  = elg['fiber'].astype(np.int32)
elg_mjd  = elg['mjd'].astype(np.float64)
elg_fail = (~elg['zsuc']).astype(np.float64)

print(f'  LRG {len(lrg_fib):,}  ELG {len(elg_fib):,}  ({time.time()-t0:.1f}s)', flush=True)

# per-fiber total obs and total failures (all time)
N = 5000
lrg_ntot  = np.bincount(lrg_fib, minlength=N).astype(np.float64)
lrg_nfail = np.bincount(lrg_fib, weights=lrg_fail, minlength=N)
elg_ntot  = np.bincount(elg_fib, minlength=N).astype(np.float64)
elg_nfail = np.bincount(elg_fib, weights=elg_fail, minlength=N)

lrg_mean = np.where(lrg_ntot > 0, lrg_nfail / lrg_ntot, np.nan)
elg_mean = np.where(elg_ntot > 0, elg_nfail / elg_ntot, np.nan)

# shared bins from global MJD range across both tracers
mjd_min = min(lrg_mjd.min(), elg_mjd.min())
mjd_max = max(lrg_mjd.max(), elg_mjd.max())
bins = np.linspace(mjd_min, mjd_max, NBINS + 1)
bin_centers = 0.5 * (bins[:-1] + bins[1:])

# fibers with enough data in both tracers
valid = (lrg_ntot >= MIN_OBS) & (elg_ntot >= MIN_OBS) & \
        (lrg_mean > 0) & (lrg_mean < 1) & \
        (elg_mean > 0) & (elg_mean < 1)

print(f'  {valid.sum()} fibers with >= {MIN_OBS} obs in both tracers', flush=True)

results = []

for b in range(NBINS):
    b0, b1 = bins[b], bins[b+1]

    # LRG: obs and failures per fiber in this bin
    lrg_mask = (lrg_mjd >= b0) & (lrg_mjd < b1)
    ln = np.bincount(lrg_fib[lrg_mask], minlength=N).astype(np.float64)
    lk = np.bincount(lrg_fib[lrg_mask], weights=lrg_fail[lrg_mask], minlength=N)

    # ELG
    elg_mask = (elg_mjd >= b0) & (elg_mjd < b1)
    en = np.bincount(elg_fib[elg_mask], minlength=N).astype(np.float64)
    ek = np.bincount(elg_fib[elg_mask], weights=elg_fail[elg_mask], minlength=N)

    # z-scores (only where n >= 3)
    lrg_denom = np.sqrt(ln * lrg_mean * (1 - lrg_mean))
    elg_denom = np.sqrt(en * elg_mean * (1 - elg_mean))

    z_lrg = np.where((ln >= 3) & (lrg_denom > 0), (lk - ln * lrg_mean) / lrg_denom, np.nan)
    z_elg = np.where((en >= 3) & (elg_denom > 0), (ek - en * elg_mean) / elg_denom, np.nan)

    # fibers where both exceed threshold
    hit = valid & (z_lrg > ZSIG) & (z_elg > ZSIG)
    for fib in np.where(hit)[0]:
        results.append({
            'fiber': int(fib),
            'petal': int(fib) // 500,
            'bin_mjd': bin_centers[b],
            'z_lrg': z_lrg[fib],
            'z_elg': z_elg[fib],
            'z_prod': z_lrg[fib] * z_elg[fib],
            'n_lrg': int(ln[fib]),
            'n_elg': int(en[fib]),
            'mean_fr_lrg': lrg_mean[fib],
            'mean_fr_elg': elg_mean[fib],
        })

results.sort(key=lambda r: -r['z_prod'])
print(f'\nFound {len(results)} (fiber, bin) pairs  |  '
      f'{len(set(r["fiber"] for r in results))} unique fibers  '
      f'({time.time()-t0:.1f}s)\n', flush=True)

hdr = (f'{"fiber":>6}  {"petal":>5}  {"bin_MJD":>9}  '
       f'{"z_LRG":>7}  {"z_ELG":>7}  {"z_prod":>7}  '
       f'{"n_LRG":>6}  {"n_ELG":>6}  {"fr_LRG":>7}  {"fr_ELG":>7}')
print(hdr)
print('-' * len(hdr))
for r in results:
    print(f'{r["fiber"]:>6}  {r["petal"]:>5}  {r["bin_mjd"]:>9.1f}  '
          f'{r["z_lrg"]:>7.2f}  {r["z_elg"]:>7.2f}  {r["z_prod"]:>7.1f}  '
          f'{r["n_lrg"]:>6}  {r["n_elg"]:>6}  '
          f'{r["mean_fr_lrg"]:>7.3f}  {r["mean_fr_elg"]:>7.3f}')

outpath = f'{DATADIR}/correlated_bad_bins_lrg_elg.txt'
with open(outpath, 'w') as f:
    f.write(f'# Fibers with LRG and ELG both >{ZSIG}sigma above mean fail rate in same time bin\n')
    f.write(f'# NBINS={NBINS}, MIN_OBS={MIN_OBS}, MJD range [{mjd_min:.1f}, {mjd_max:.1f}]\n')
    f.write(f'# {hdr}\n')
    for r in results:
        f.write(f'{r["fiber"]:>6}  {r["petal"]:>5}  {r["bin_mjd"]:>9.1f}  '
                f'{r["z_lrg"]:>7.2f}  {r["z_elg"]:>7.2f}  {r["z_prod"]:>7.1f}  '
                f'{r["n_lrg"]:>6}  {r["n_elg"]:>6}  '
                f'{r["mean_fr_lrg"]:>7.3f}  {r["mean_fr_elg"]:>7.3f}\n')
print(f'\nSaved {outpath}')
