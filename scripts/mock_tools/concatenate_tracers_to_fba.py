from astropy.table import Table,vstack
import LSS.common_tools as cm
from astropy.io import fits
import numpy as np

val = 0
list_of_realizations = [500,601,602,603,604,605,606,607,608,609,610,611,612,613,614,615,616,618,619,620,621,622,623,624,625,626,627,628,629,630,631,632,633,634,635,636,637,638,639,640,641,642,643,644,645,646,648,649,650] 
for seed in list_of_realizations[4:]:
    realization = str(seed).zfill(4)
    seednum = f'seed{realization}'
    elg_path = f'/global/cfs/cdirs/desi/mocks/cai/holi/v5.0/{seednum}/ELG/forFA0_withcontaminants.fits'
    elg1 = Table.read(elg_path)

    lrg_path = f'/global/cfs/cdirs/desi/mocks/cai/holi/v4.00/{seednum}/LRG/forFA0.fits'
    lrg1 = Table.read(lrg_path)

    qso_path = f'/global/cfs/cdirs/desi/mocks/cai/holi/v4.00/{seednum}/QSO_v2/forFA0_withcontaminants.fits'
    qso1 = Table.read(qso_path)

    
    qsofile = f'/pscratch/sd/d/desica/DA2/mocks/holi_v1/qsos/qso{seed}.txt'

    np.savetxt(qsofile, np.array([qso1['TARGETID'], qso1['RSDZ']]).T, fmt='%d %.3f')
    print(f'saving qsos to {qsofile}')

    targets = vstack([elg1, lrg1, qso1])

    del elg1
    del lrg1
    del qso1
    '''
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
    '''
    cm.write_LSS_scratchcp(targets, f'/pscratch/sd/d/desica/DA2/mocks/holi_v1/forFA{seed}.fits', extname='TARGETS')
    fits.setval(f'/pscratch/sd/d/desica/DA2/mocks/holi_v1/forFA{seed}.fits', 'OBSCON', value='DARK', ext=1)

    print('done', seed)
    val += 1
    del targets
