from astropy.table import Table,vstack
import LSS.common_tools as cm
from astropy.io import fits
import numpy as np

val = 0
for seed in range(451,500):
    if seed == 477:
        continue

    seednum = f'seed0{seed}'
    elg_path = f'/global/cfs/cdirs/desi/mocks/cai/holi/v5.0/{seednum}/ELG/forFA0_withcontaminants.fits'
    elg1 = Table.read(elg_path)

    lrg_path = f'/global/cfs/cdirs/desi/mocks/cai/holi/v4.00/{seednum}/LRG/forFA0.fits'
    lrg1 = Table.read(lrg_path)

    qso_path = f'/global/cfs/cdirs/desi/mocks/cai/holi/v4.00/{seednum}/QSO/forFA0_withcontaminants.fits'
    qso1 = Table.read(qso_path)

    
    
    qsofile = f'/global/cfs/cdirs/desi/mocks/cai/holi/altMTL/qsos/qso{val}.txt'

    np.savetxt(qsofile, np.array([qso1['TARGETID'], qso1['RSDZ']]).T, fmt='%d %.3f')
    print(f'saving qsos to {qsofile}')


    targets = vstack([elg1, lrg1, qso1])

    del elg1
    del lrg1
    del qso1

    seen = set()
    has_duplicates = False
    for valo in targets['TARGETID']:
        if valo in seen:
            has_duplicates = True
            break
        seen.add(valo)

    print("All unique?", not has_duplicates)

    for name, col in targets.columns.items():
        print(f"{name}: {col.dtype}")

    cm.write_LSS_scratchcp(targets, f'/global/cfs/cdirs/desi/mocks/cai/holi/altMTL/forFA{val}.fits', extname='TARGETS')
    fits.setval(f'/global/cfs/cdirs/desi/mocks/cai/holi/altMTL/forFA{val}.fits', 'OBSCON', value='DARK', ext=1)

    print('done', seed)
    val += 1
    del targets
