from astropy.table import Table,vstack
import LSS.common_tools as cm
from astropy.io import fits
import numpy as np

val = 0
list_of_realizations = [201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 477, 617, 647]

for seed in list_of_realizations:
    realization = str(seed).zfill(4)
    seednum = f'seed{realization}'
    elg_path = f'/global/cfs/cdirs/desi/mocks/cai/holi/v5.0_Y5/{seednum}/ELG/forFA0_withcontaminants.fits'
    elg1 = Table.read(elg_path)

    lrg_path = f'/global/cfs/cdirs/desi/mocks/cai/holi/v4.00/{seednum}/LRG/forFA0.fits'
    lrg1 = Table.read(lrg_path)

    qso_path = f'/global/cfs/cdirs/desi/mocks/cai/holi/v4.00/{seednum}/QSO/forFA0_withcontaminants.fits'
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
