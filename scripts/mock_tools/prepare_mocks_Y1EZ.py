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
parser.add_argument("--realization", help="number for the realization",default='0')
parser.add_argument("--prog", help="dark or bright",default='dark')
parser.add_argument("--prep", help="prepare file for fiberassign?",default='y')
args = parser.parse_args()


if args.prog == 'dark':
    types = ['ELG', 'LRG', 'QSO']
    desitar = {'ELG':34,'LRG':1,'QSO':4}
    priority = {'ELG':3000,'LRG':3200,'QSO':3400}
    mainp = main(tp='QSO',specver='iron')
    numobs = {'ELG':2, 'LRG':2, 'QSO':4}

inroot = '/global/cfs/cdirs/desi/survey/catalogs/'
inmock = '/Y1/mocks/SecondGenMocks/EZmock/'
indir = inroot+inmock

outroot = os.getenv('SCRATCH')
outdir = outroot+inmock+'/forFA/'
if not os.path.exists(outdir):
	os.makedirs(outdir)
out_file_name = outdir + 'forFA'+args.realization+'.fits'

status = 3
datat = []
percentage_elg_hip = 0.1

for type_ in types:
    fname = 'EZmock_'+type_+'_complete_AbacusSummit_base_c000_ph000_NScomb_'+args.realization.zfill(4)+'.fits.gz'
    data = fitsio.read(indir+type_+'/'+fname,columns=['RA','DEC','Z','Z_COSMO','STATUS','NOBS_G','NOBS_R','NOBS_Z','MASKBITS'])#f[1].data
    sel = data['STATUS'] == 3
    data = data[sel]
    data = common.cutphotmask(data,bits=mainp.imbits) #already done?
    data = Table(data)
    if type_ == 'ELG':
        import sample_elg_ezmock as se
        lop, vlo = se.create_subsample(data)

        data_lop = Table.from_pandas(lop)
        data_vlo = Table.from_pandas(vlo)

        df_lop=data_lop.to_pandas()
        df_vlo=data_vlo.to_pandas()
        num_HIP_LOP = int(len(df_lop) * percentage_elg_hip)
        df_HIP_LOP = df_lop.sample(n=num_HIP_LOP)
        remaining_LOP = df_lop.drop(df_HIP_LOP.index)
        df_HIP_LOP.reset_index(drop=True, inplace=True)
        remaining_LOP.reset_index(drop=True, inplace=True)

        num_HIP_VLO = int(len(df_vlo) * percentage_elg_hip)
        df_HIP_VLO = df_vlo.sample(n=num_HIP_VLO)
        remaining_VLO = df_vlo.drop(df_HIP_VLO.index)
        df_HIP_VLO.reset_index(drop=True, inplace=True)
        remaining_VLO.reset_index(drop=True, inplace=True)

        remaining_LOP['PRIORITY_INIT'] = 3100
        remaining_LOP['PRIORITY'] = 3100
        remaining_LOP['DESI_TARGET'] = 2**5 + 2**1
        remaining_VLO['PRIORITY_INIT'] = 3000
        remaining_VLO['PRIORITY'] = 3000
        remaining_VLO['DESI_TARGET'] = 2**7 + 2**1

        df_HIP_LOP['PRIORITY_INIT'] = 3200
        df_HIP_LOP['PRIORITY'] = 3200
        df_HIP_LOP['DESI_TARGET'] = 2**6 + 2**1 + 2**5

        df_HIP_VLO['PRIORITY_INIT'] = 3200
        df_HIP_VLO['PRIORITY'] = 3200
        df_HIP_VLO['DESI_TARGET'] = 2**6 + 2**1 + 2**5

        remaining_LOP['NUMOBS_MORE'] = numobs[type_]
        remaining_LOP['NUMOBS_INIT'] = numobs[type_]
        remaining_VLO['NUMOBS_MORE'] = numobs[type_]
        remaining_VLO['NUMOBS_INIT'] = numobs[type_]
        df_HIP_LOP['NUMOBS_MORE'] = numobs[type_]
        df_HIP_LOP['NUMOBS_INIT'] = numobs[type_]
        df_HIP_VLO['NUMOBS_MORE'] = numobs[type_]
        df_HIP_VLO['NUMOBS_INIT'] = numobs[type_]

        datat.append(Table.from_pandas(remaining_LOP))
        datat.append(Table.from_pandas(remaining_VLO))
        datat.append(Table.from_pandas(df_HIP_LOP))
        datat.append(Table.from_pandas(df_HIP_VLO))
    else:

        data['DESI_TARGET'] = desitar[type_]
        data['PRIORITY_INIT'] = priority[type_]
        data['PRIORITY'] = priority[type_]
        data['NUMOBS_MORE'] = numobs[type_] 
        data['NUMOBS_INIT'] = numobs[type_]
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
targets['SCND_TARGET'] = np.zeros(n, dtype='i8')+int(0)
targets['ZWARN'] = np.zeros(n, dtype='i8')+int(0)
targets['TARGETID'] = np.arange(1,n+1)

targets.write(out_file_name, overwrite = True)

fits.setval(out_file_name, 'EXTNAME', value='TARGETS', ext=1)
fits.setval(out_file_name, 'OBSCON', value=args.prog.upper(), ext=1)




sys.exit()

