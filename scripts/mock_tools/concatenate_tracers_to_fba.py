from astropy.table import Table,vstack
import LSS.common_tools as cm
from astropy.io import fits
import numpy as np

elg1 = Table.read('/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v1.0/AbacusSummit_base_c000_ph000/CutSky/ELG_v5/z0.950/forFA0_Y3_noimagingmask_applied_testfine_withcontaminants.fits')
#elg1 = Table.read('/global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/Holi/seed0202/elgv3.00_nomask.fits')
#mask1 = elg1['RSDZ']<=1.1
#elg2 = Table.read('/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v1.0/AbacusSummit_base_c000_ph000/CutSky/ELG_v5/z1.175/forFA0_Y3.fits')
#mask2 = elg2['RSDZ']>1.1


#E1 = elg1[mask1]
#E2 = elg2[mask2]
lrg1 = Table.read('/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v1.0/AbacusSummit_base_c000_ph000/CutSky/LRG_v4/z0.500/forFA0_base_Y3_noimagingmask_applied_testfine.fits')
#lrg1 = Table.read('/global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/Holi/seed0202/lrgv3.13_nomask.fits')
#lrg2 = Table.read('/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v1.0/AbacusSummit_base_c000_ph000/CutSky/LRG/z0.725/forFA0_Y3.fits')

#mask1 = lrg1['RSDZ']<=0.6
#mask2 = lrg2['RSDZ']>0.6
#L1 = lrg1[mask1]
#L2 = lrg2[mask2]
Q = Table.read('/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v1.0/AbacusSummit_base_c000_ph000/CutSky/QSO_v5/z1.400/forFA0_Y3_noimagingmask_applied_testfine_withcontaminants.fits')
#Q = Table.read('/global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/Holi/seed0202/qsov3.00_nomask.fits')
targets = vstack([elg1, lrg1, Q])
#targets = vstack([E1,E2,L1,L2,Q])

del elg1
#del elg2
del lrg1
#del lrg2
#del E1
#del E2
#del L1
#del L2
del Q


##n=len(targets)  ##A Ashley le falta estoo!
##targets['TARGETID'] = (np.arange(1,n+1)+1e8).astype(int) #different tracer types need to have different targetids
#print('start writing')



seen = set()
has_duplicates = False
for val in targets['TARGETID']:
    if val in seen:
        has_duplicates = True
        break
    seen.add(val)

print("All unique?", not has_duplicates)

for name, col in targets.columns.items():
    print(f"{name}: {col.dtype}")

cm.write_LSS_scratchcp(targets, '/global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/AbacusHighFidelity/forFA0_noimagingmask_applied_testfine_withcontaminants.fits', extname='TARGETS')
#cm.write_LSS_scratchcp(targets, '/global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/Holi/seed0202/forFA202_noimagingmask.fits', extname='TARGETS')
fits.setval('/global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/AbacusHighFidelity/forFA0_noimagingmask_applied_testfine_withcontaminants.fits', 'OBSCON', value='DARK', ext=1)
#fits.setval('/global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/Holi/seed0202/forFA202_noimagingmask.fits', 'OBSCON', value='DARK', ext=1)

#targets.write('/global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/AbacusHighFidelity/forFA0.fits')

