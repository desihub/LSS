from astropy.table import Table, join
import LSS.common_tools as cm
from LSS.globals import main
#from astropy.io import fits
import os
import numpy as np


mainp = main(tp = 'LRG', specver = 'loa-v1')

path_parent = '/global/homes/d/desica/BRICKMASKcode/brickmask'
files = ['out_list_elg.txt']




for ff in list_of_files:
    ifil_list = np.loadtxt(os.path.join(path_parent, ff), unpack=True, dtype=str)
    for ifil in ifil_list:
        d=Table.read(ifil)

        targets = cm.cutphotmask(d, bits=mainp.imbits)
        directory = os.path.dirname(ifil)
        cm.write_LSS(targets, os.path.join(directory, 'forFA0.fits'), extname='TARGETS')
        print(os.path.join(directory, 'forFA0.fits'), 'done')


