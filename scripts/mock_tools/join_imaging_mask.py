from astropy.table import Table, join
import LSS.common_tools as cm
from LSS.globals import main
from astropy.io import fits

mainp = main(tp = 'LRG', specver = 'loa-v1')

d=Table.read('/global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/AbacusHighFidelity/forFA0_noimagingmask_applied_testfine_withcontaminants_maskbitsnobs.fits')
targets = cm.cutphotmask(d, bits=mainp.imbits)

cat = Table.read('/global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/AbacusHighFidelity/forFA0_noimagingmask_applied_testfine_withcontaminants.fits')

matched = join(cat, targets, keys='TARGETID', join_type='inner')
cm.write_LSS(matched, '/global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/AbacusHighFidelity/forFA0.fits', extname='TARGETS')
fits.setval('/global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/AbacusHighFidelity/forFA0.fits', 'OBSCON', value='DARK', ext=1)


