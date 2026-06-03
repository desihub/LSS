"""
make_fiber_time_data.py
Build per-tracer data files for failure-rate-vs-time plots.

For each tracer (BGS, LRG, ELG, QSO) writes:
  <tracer>_fiber_time.npz  containing:
    fiber  int16  (N,)   fiber number
    mjd    float32(N,)   mean tile MJD
    zsuc   bool   (N,)   True = good redshift
    index  int32  (5000,2) start/end row for each fiber (end exclusive)
                           index[f] = [-1,-1] if fiber f has no entries

Output dir: same directory as this script.
"""

import numpy as np
import fitsio
import datetime
import time
from collections import defaultdict

OUTDIR  = '/global/u1/r/rohlf/bao_jr/dr3-matterhorn-checks/all-fibers-vs-time'
LSSDIR  = '/global/cfs/cdirs/desi/survey/catalogs/DA3/LSS/matterhorn-v2/LSScats/test'
EXPCAT  = '/global/cfs/cdirs/desi/spectro/redux/matterhorn/exposures-matterhorn.fits'

TRACERS = {
    'LRG': {
        'cat':   f'{LSSDIR}/LRG_full_noveto.dat.fits',
        'tsnr':  'TSNR2_ELG',
        'tsnr_cut': 80,
        'dchi2_cut': 15,
        'zmin': 0.0, 'zmax': 1.5,
    },
    'ELG': {
        'cat':   f'{LSSDIR}/ELG_LOPnotqso_full_noveto.dat.fits',
        'tsnr':  'TSNR2_ELG',
        'tsnr_cut': 80,
        'dchi2_cut': 9,
        'zmin': 0.6, 'zmax': 1.6,
    },
    'BGS': {
        'cat':   f'{LSSDIR}/BGS_BRIGHT_full_noveto.dat.fits',
        'tsnr':  'TSNR2_BGS',
        'tsnr_cut': 1000,
        'dchi2_cut': 25,
        'zmin': 0.001, 'zmax': 0.6,
    },
    'QSO': {
        'cat':   f'{LSSDIR}/QSO_full_noveto.dat.fits',
        'tsnr':  'TSNR2_ELG',
        'tsnr_cut': 80,
        'dchi2_cut': 25,
        'zmin': 0.8, 'zmax': 3.5,
    },
}

# ── exposure catalog: TILEID -> mean MJD ─────────────────────────────────────
print('Loading exposure catalog...', flush=True)
t0 = time.time()
exp = fitsio.read(EXPCAT, columns=['TILEID', 'MJD'])
tile_mjd_dict = defaultdict(list)
for r in exp:
    tile_mjd_dict[int(r['TILEID'])].append(float(r['MJD']))
# vectorizable lookup arrays
max_tileid = max(tile_mjd_dict.keys()) + 1
tile_mjd_arr = np.full(max_tileid, np.nan, dtype=np.float32)
for tid, mjds in tile_mjd_dict.items():
    tile_mjd_arr[tid] = np.mean(mjds)
print(f'  {len(tile_mjd_dict):,} tiles  ({time.time()-t0:.1f}s)', flush=True)

# ── process each tracer ───────────────────────────────────────────────────────
for tracer, cfg in TRACERS.items():
    print(f'\n── {tracer} ──', flush=True)
    t0 = time.time()

    cols = ['FIBER', 'TILEID', 'ZWARN', 'DELTACHI2', 'Z',
            'COADD_FIBERSTATUS', cfg['tsnr']]
    cat = fitsio.read(cfg['cat'], columns=cols)
    print(f'  read {len(cat):,} rows  ({time.time()-t0:.1f}s)', flush=True)

    fiber_arr = np.array(cat['FIBER'],             dtype=np.int16)
    tileid_arr= np.array(cat['TILEID'],            dtype=np.int32)
    zwarn_arr = np.array(cat['ZWARN'],             dtype=np.int32)
    dchi2_arr = np.array(cat['DELTACHI2'],         dtype=np.float32)
    z_arr     = np.array(cat['Z'],                 dtype=np.float32)
    fibst_arr = np.array(cat['COADD_FIBERSTATUS'], dtype=np.int32)
    tsnr_arr  = np.array(cat[cfg['tsnr']],         dtype=np.float32)
    del cat

    # z_tot mask
    z_tot = ((zwarn_arr != 999999) &
             ((fibst_arr == 0) | (fibst_arr == 8)) &
             (tsnr_arr > cfg['tsnr_cut']))

    # z_suc mask
    z_suc = ((zwarn_arr == 0) &
             (dchi2_arr > cfg['dchi2_cut']) &
             (z_arr > cfg['zmin']) &
             (z_arr < cfg['zmax']))

    # apply z_tot filter
    fiber_f  = fiber_arr[z_tot]
    tileid_f = tileid_arr[z_tot]
    zsuc_f   = z_suc[z_tot]

    # vectorized MJD lookup
    tid_clipped = np.clip(tileid_f, 0, len(tile_mjd_arr) - 1)
    mjd_f = tile_mjd_arr[tid_clipped]

    # drop rows with no MJD
    good = np.isfinite(mjd_f)
    fiber_f = fiber_f[good]
    mjd_f   = mjd_f[good]
    zsuc_f  = zsuc_f[good]

    print(f'  {len(fiber_f):,} rows after z_tot + MJD filter  ({time.time()-t0:.1f}s)', flush=True)

    # sort by fiber
    order   = np.argsort(fiber_f, kind='stable')
    fiber_f = fiber_f[order]
    mjd_f   = mjd_f[order]
    zsuc_f  = zsuc_f[order]

    # build index: index[f] = [start, end) row for fiber f; [-1,-1] if absent
    index = np.full((5000, 2), -1, dtype=np.int32)
    unique_fibers, starts, counts = np.unique(fiber_f, return_index=True, return_counts=True)
    for uf, s, c in zip(unique_fibers, starts, counts):
        if 0 <= uf < 5000:
            index[uf] = [s, s + c]

    n_fibers = (index[:, 0] >= 0).sum()
    print(f'  {n_fibers} unique fibers covered', flush=True)

    outpath = f'{OUTDIR}/{tracer}_fiber_time.npz'
    np.savez_compressed(outpath,
                        fiber=fiber_f,
                        mjd=mjd_f,
                        zsuc=zsuc_f,
                        index=index)
    import os
    size_mb = os.path.getsize(outpath) / 1e6
    print(f'  Saved {outpath}  ({size_mb:.1f} MB)  total {time.time()-t0:.1f}s', flush=True)

print('\nDone.', flush=True)
