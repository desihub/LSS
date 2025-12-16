from astropy.table import Table, join
import LSS.common_tools as cm
from LSS.globals import main
#from astropy.io import fits
import os

##pa_elg = '/global/cfs/cdirs/desi/mocks/cai/holi/v5.0/seed0{id_}/ELG/'

mainp = main(tp = 'LRG', specver = 'loa-v1')

pa = '/global/cfs/cdirs/desi/mocks/cai/holi/v4.00/seed0{id_}/{ty}/'

for i in range(450,500):
    for tye in ['QSO', 'LRG']:
        ff = os.path.join(pa, 'imagingmask.fits').format(id_=i, ty=tye)
        if not os.path.isfile(ff):
            print(ff, 'dont exist')
            continue

        d=Table.read(ff)

        targets = cm.cutphotmask(d, bits=mainp.imbits)

        cm.write_LSS(targets, os.path.join(pa, 'forFA0.fits').format(id_=i, ty=tye), extname='TARGETS')
        print(os.path.join(pa, 'forFA0.fits').format(id_=i, ty=tye), 'done')
#fits.setval('/global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/AbacusHighFidelity/forFA0.fits', 'OBSCON', value='DARK', ext=1)


