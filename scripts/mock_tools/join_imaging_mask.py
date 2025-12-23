from astropy.table import Table, join
import LSS.common_tools as cm
from LSS.globals import main
#from astropy.io import fits
import os
import numpy as np
from glob import glob

mainp = main(tp = 'LRG', specver = 'loa-v1')

path_parent = '/global/homes/d/desica/LSScode/LSS/scripts/mock_tools'
#files = ['imaging_lrg1.out', 'imaging_lrg2.out', 'imaging_qso1.out', 'imaging_qso2.out']
files = ['out_qso_list_aa.txt', 'out_qso_list_ab.txt', 'out_qso_list_ac.txt', 'out_qso_list_ad.txt', 'out_qso_list_ae.txt', 'out_qso_list_af.txt']



for ff in files:
    ifil_list = np.loadtxt(os.path.join(path_parent, ff), unpack=True, dtype=str)
    for ifil in ifil_list:

        output_path = os.path.join(os.path.dirname(ifil), 'forFA0.fits')
        if os.path.isfile(output_path):
            continue
        else:
            pass
        
        d=Table.read(glob(ifil)[0])

        targets = cm.cutphotmask(d, bits=mainp.imbits)
        directory = os.path.dirname(ifil)
        cm.write_LSS(targets, os.path.join(directory, 'forFA0.fits'), extname='TARGETS')
        print(os.path.join(directory, 'forFA0.fits'), 'done')


