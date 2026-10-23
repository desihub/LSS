"""
Compare column statistics (mean, std, min, max) between the DR3 daily catalog
and matterhorn-v2, per tracer. Compares full.dat.fits (pre-clustering-cut)
from both test/ directories. Prints a formatted table.

Authors: James Rohlf and Claude Sonnet 4.6
Last revised: 2026-05-21
"""
import os
import numpy as np
from astropy.io import fits

DR3_DAILY      = '/global/cfs/cdirs/desi/survey/catalogs/DA3/LSS/daily/LSScats/test'
DR3_MATTERHORN = '/global/cfs/cdirs/desi/survey/catalogs/DA3/LSS/matterhorn-v2/LSScats/test'

TRACERS = [
    ('LRG',              0.4, 1.1),
    ('ELG_LOPnotqso',    0.8, 1.6),
    ('ELGnotqso',        0.8, 1.6),
    ('QSO',              0.8, 2.1),
    ('BGS_BRIGHT',       0.1, 0.4),
    ('BGS_FAINT',        0.1, 0.4),
    ('BGS_ANY',          0.1, 0.4),
    ('ELG_LOP',          0.8, 1.6),
]

COLS = ['Z_not4clus', 'WEIGHT_ZFAIL', 'FRAC_TLOBS_TILES', 'COMP_TILE',
        'FRACZ_TILELOCID', 'mod_success_rate', 'TSNR2_LRG', 'TSNR2_ELG',
        'TSNR2_QSO']


def get_stats(catdir, tracer, zmin, zmax):
    fpath = os.path.join(catdir, f'{tracer}_full.dat.fits')
    if not os.path.exists(fpath):
        return None
    with fits.open(fpath) as h:
        d = h[1].data
    mask = (d['Z_not4clus'] > zmin) & (d['Z_not4clus'] < zmax)
    nall = len(d)
    ncut = mask.sum()
    stats = {'N_all': nall, 'N_zcut': ncut}
    for col in COLS:
        if col not in d.columns.names:
            stats[col] = None
            continue
        arr = d[col][mask].astype('f8')
        stats[col] = (arr.mean(), arr.std(), arr.min(), arr.max())
    return stats


for tracer, zmin, zmax in TRACERS:
    label = f'{tracer}  z={zmin}-{zmax}'
    print(f'\n{"="*72}')
    print(f'  {label}')
    print(f'{"="*72}')

    sd = get_stats(DR3_DAILY,      tracer, zmin, zmax)
    sm = get_stats(DR3_MATTERHORN, tracer, zmin, zmax)

    if sd is None: print('  DR3 daily: FILE MISSING'); continue
    if sm is None: print('  matterhorn-v2: FILE MISSING'); continue

    W = 22
    print(f'  {"":^{W}}  {"DR3 daily":^36}  {"matterhorn-v2":^36}')
    print(f'  {"":^{W}}  {"mean":>9} {"sigma":>9} {"range":>16}  '
          f'{"mean":>9} {"sigma":>9} {"range":>16}')
    print(f'  {"-"*W}  {"-"*36}  {"-"*36}')

    for col in COLS:
        if sd[col] is None and sm[col] is None:
            print(f'  {col:{W}}  {"(missing in both)":^36}')
            continue
        if sd[col] is None:
            mm, smv, lo_m, hi_m = sm[col]
            rng_m = f'[{lo_m:.4g}, {hi_m:.4g}]'
            print(f'  {col:{W}}  {"(missing in daily)":^36}  '
                  f'{mm:9.4f} {smv:9.4f} {rng_m:>16}')
            continue
        if sm[col] is None:
            md, sdv, lo_d, hi_d = sd[col]
            rng_d = f'[{lo_d:.4g}, {hi_d:.4g}]'
            print(f'  {col:{W}}  {md:9.4f} {sdv:9.4f} {rng_d:>16}  '
                  f'{"(missing in matterhorn)":^36}')
            continue
        md, sdv, lo_d, hi_d = sd[col]
        mm, smv, lo_m, hi_m = sm[col]
        rng_d = f'[{lo_d:.4g}, {hi_d:.4g}]'
        rng_m = f'[{lo_m:.4g}, {hi_m:.4g}]'
        print(f'  {col:{W}}  {md:9.4f} {sdv:9.4f} {rng_d:>16}  '
              f'{mm:9.4f} {smv:9.4f} {rng_m:>16}')

    print(f'\n  N (all rows):  daily={sd["N_all"]:>10,}   matterhorn={sm["N_all"]:>10,}'
          f'   ratio={sm["N_all"]/sd["N_all"]:.4f}x')
    print(f'  N (z-cut):     daily={sd["N_zcut"]:>10,}   matterhorn={sm["N_zcut"]:>10,}'
          f'   ratio={sm["N_zcut"]/sd["N_zcut"]:.4f}x')

print(f'\n{"="*72}')
print('Done.')
