from astropy.io import fits # Access to FITS (Flexible Image Transport System) files.
from astropy.table import Table, hstack, vstack, Column # A class to represent tables of heterogeneous data.
import fitsio
import numpy as np
import glob
import os
import h5py
import argparse
import sys

from desitarget.targetmask import obsconditions
from desimodel.footprint import is_point_in_desi

import LSS.common_tools as common
from LSS.imaging import get_pixel_bitmasknobs as bitmask #get_nobsandmask
from LSS.main.cattools import count_tiles_better
from LSS.globals import main


if os.environ['NERSC_HOST'] == 'cori':
    scratch = 'CSCRATCH'
elif os.environ['NERSC_HOST'] == 'perlmutter':
    scratch = 'PSCRATCH'
else:
    print('NERSC_HOST is not cori or permutter but is '+os.environ['NERSC_HOST'])
    sys.exit('NERSC_HOST not known (code only works on NERSC), not proceeding') 

parser = argparse.ArgumentParser()
parser.add_argument("--realization", help="number for the realization",default=1,type=int)
parser.add_argument("--prog", help="dark or bright",default='dark')
parser.add_argument("--prep", help="prepare file for fiberassign?",default='y')
args = parser.parse_args()


if args.prog == 'dark':
    types = ['ELG', 'LRG', 'QSO']
    desitar = {'ELG':34,'LRG':1,'QSO':4}
    priority = {'ELG':3000,'LRG':3200,'QSO':3400}
    mainp = main(tp='QSO',specver='iron')

indir = '/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/EZmock/'

status = 3
datat = []
for type_ in types:
    fname = 'EZmock_'+type_+'_complete_AbacusSummit_base_c000_ph000_NScomb_'+args.realization.zfill(4)+'.fits.gz'
    data = fitsio.read(indir+fname,columns=['RA','DEC','Z','Z_COSMO','STATUS','NOBS_G','NOBS_R','NOBS_Z','MASKBITS'])#f[1].data
    sel = data['STATUS'] == 3
    data = data[sel]
    data = data.cutphotmask(targets,bits=mainp.imbits) #already done?
    data = Table(data)
    data['DESI_TARGET'] = desitar[type_]
    data['PRIORITY_INIT'] = priority[type_]
    data['PRIORITY'] = priority[type_]
    datat.append(data)
targets = vstack(datat)
del datat
n=len(targets)
targets.rename_column('Z_COSMO', 'TRUEZ') 
targets.rename_column('Z', 'RSDZ') 
targets['BGS_TARGET'] = np.zeros(n, dtype='i8')
targets['MWS_TARGET'] = np.zeros(n, dtype='i8')
targets['SUBPRIORITY'] = np.random.uniform(0, 1, n)
targets['BRICKNAME'] = np.full(n, '000p0000')    #- required !?!
targets['OBSCONDITIONS'] = obsconditions.mask(args.prog.upper()) #np.zeros(n, dtype='i8')+int(3) 
targets['NUMOBS_MORE'] = np.zeros(n, dtype='i8')+int(1) 
targets['NUMOBS_INIT'] = np.zeros(n, dtype='i8')+int(1)
targets['SCND_TARGET'] = np.zeros(n, dtype='i8')+int(0)
targets['ZWARN'] = np.zeros(n, dtype='i8')+int(0)
targets['TARGETID'] = np.arange(1,n+1)

targets.write(out_file_name, overwrite = True)

fits.setval(out_file_name, 'EXTNAME', value='TARGETS', ext=1)
fits.setval(out_file_name, 'OBSCON', value=args.prog.upper(), ext=1)




sys.exit()

