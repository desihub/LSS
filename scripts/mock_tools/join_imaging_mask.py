from astropy.table import Table, join
import LSS.common_tools as cm
from LSS.globals import main
#from astropy.io import fits
import os
import numpy as np

##pa_elg = '/global/cfs/cdirs/desi/mocks/cai/holi/v5.0/seed0{id_}/ELG/'

mainp = main(tp = 'LRG', specver = 'loa-v1')

path_parent = '/pscratch/sd/a/acarnero/codes/brickmask'
files = ['out_list_elg.txt']
#files = ['out_list_lrg.txt', 'out_list_ac.txt', 'out_list_ab.txt', 'out_list_aa.txt']



#pa = '/global/cfs/cdirs/desi/mocks/cai/holi/v4.00/seed0{id_}/{ty}/'

#for i in range(450,500):
#    for tye in ['QSO', 'LRG']:
#        ff = os.path.join(pa, 'imagingmask.fits').format(id_=i, ty=tye)
#        if not os.path.isfile(ff):
#            print(ff, 'dont exist')
#            continue
for ff in files:
    ifil_list = np.loadtxt(os.path.join(path_parent, ff), unpack=True, dtype=str)
    for ifil in ifil_list:
        d=Table.read(ifil)

        targets = cm.cutphotmask(d, bits=mainp.imbits)
        directory = os.path.dirname(ifil)
        cm.write_LSS(targets, os.path.join(directory, 'forFA0.fits'), extname='TARGETS')
        print(os.path.join(directory, 'forFA0.fits'), 'done')
#fits.setval('/global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/AbacusHighFidelity/forFA0.fits', 'OBSCON', value='DARK', ext=1)


