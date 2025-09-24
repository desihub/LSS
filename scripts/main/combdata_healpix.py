#standard python
import sys
import os
import shutil
import unittest
from datetime import datetime
import json
import numpy as np
import healpy as hp
import fitsio
import glob
import argparse
from astropy.table import Table,join,unique,vstack
from matplotlib import pyplot as plt
from desitarget.io import read_targets_in_tiles
from desitarget.mtl import inflate_ledger
from desimodel.footprint import is_point_in_desi
import desimodel.footprint as foot
from desitarget import targetmask

#import logging
#logging.getLogger().setLevel(logging.ERROR)


#sys.path.append('../py') #this requires running from LSS/bin, *something* must allow linking without this but is not present in code yet

#from this package
#try:
import LSS.main.cattools as ct
import LSS.common_tools as common
from LSS.globals import main
from LSS.qso_cat_utils import qso_catalog_maker,build_qso_catalog_from_healpix

if os.environ['NERSC_HOST'] == 'cori':
    scratch = 'CSCRATCH'
elif os.environ['NERSC_HOST'] == 'perlmutter':
    scratch = 'PSCRATCH'
else:
    print('NERSC_HOST is not cori or permutter but is '+os.environ['NERSC_HOST'])
    sys.exit('NERSC_HOST not known (code only works on NERSC), not proceeding') 

parser = argparse.ArgumentParser()
parser.add_argument("--basedir", help="base directory for output, default is SCRATCH",default=scratch)
parser.add_argument("--version", help="catalog version; use 'test' unless you know what you are doing!",default='test')
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA",default='main')
parser.add_argument("--prog", help="dark or bright is supported",default='dark')
parser.add_argument("--verspec",help="version for redshifts",default='himalayas')
parser.add_argument("--doqso",help="whether or not to combine qso data",default='y')
parser.add_argument("--mkemlin",help="whether or not to make emission line files",default='y')




args = parser.parse_args()
print(args)

basedir = args.basedir
version = args.version
specrel = args.verspec
prog = args.prog
progu = prog.upper()


#share basedir location '/global/cfs/cdirs/desi/survey/catalogs'
maindir = basedir +'/'+args.survey+'/LSS/'

if not os.path.exists(maindir+'/logs'):
    os.mkdir(maindir+'/logs')
    print('made '+maindir+'/logs')

if not os.path.exists(maindir+'/LSScats'):
    os.mkdir(maindir+'/LSScats')
    print('made '+maindir+'/LSScats')

dirout = maindir+'LSScats/'+version+'/'
if not os.path.exists(dirout):
    os.mkdir(dirout)
    print('made '+dirout)

ldirspec = maindir+specrel+'/'
if not os.path.exists(ldirspec):
    os.mkdir(ldirspec)
    print('made '+ldirspec)
if not os.path.exists(ldirspec+'healpix'):
    os.mkdir(ldirspec+'healpix')
    print('made '+ldirspec+'healpix')


if  args.doqso == 'y':
    #outf = ldirspec+'QSO_catalog_'+prog+'_healpix_v'+args.version+'.fits'
    surpipe = 'main'
    if args.survey == 'SV3':
        surpipe = 'sv3'
    build_qso_catalog_from_healpix( release=args.verspec, survey=surpipe, program=args.prog, dir_output=ldirspec, npool=20, keep_qso_targets=False, keep_all=False,qsoversion='test')
#     dirspec = '/global/cfs/cdirs/desi/spectro/redux/'+args.verspec+'/healpix/'+surpipe+'/'+args.prog+'/'
#     
#     subdirs = os.listdir(dirspec)
#     qsocats = []
#     kl = ['TARGET_RA','TARGET_DEC','DESI_TARGET','TARGETID', 'Z',  'TSNR2_LYA', 'TSNR2_QSO', 'DELTA_CHI2_MGII', 'A_MGII', 'SIGMA_MGII', 'B_MGII', 'VAR_A_MGII', 'VAR_SIGMA_MGII', 'VAR_B_MGII', 'Z_RR', 'Z_QN', 'C_LYA', 'C_CIV', 'C_CIII', 'C_MgII', 'C_Hbeta', 'C_Halpha', 'Z_LYA', 'Z_CIV', 'Z_CIII', 'Z_MgII', 'Z_Hbeta', 'Z_Halpha', 'QSO_MASKBITS']
#     n = 0
#     for sd in subdirs:
#         fd = dirspec+sd+'/'
#         ssdir = os.listdir(fd)
#         print(sd)
#         for ssd in ssdir:
#             du = fd+ssd+'/'
#             app = surpipe+'-'+prog+'-'+ssd+'.fits'
#             rr = du+'redrock-'+app
#             mgii = du+'qso_mgii-'+app
#             qn = du+'qso_qn-'+app
#             old_extname_redrock = False
#             old_extname_for_qn = False #if int(tdate) >= 20220118 else True
#             try:
#                 qso_cati = Table.from_pandas(qso_catalog_maker(rr, mgii, qn, old_extname_redrock, old_extname_for_qn))
#                 names = list(qso_cati.dtype.names)
#                 kll = kl#[]
#                 #for i in range(0,len(kl)):
#                 #    if kl[i] in names:
#                 #        kll.append(kl[i])
#                 #    else:
#                 #        print(kl[i])
#                 qso_cati.keep_columns(kll)
#                 ‘HPXPIXEL’
#                 qsocats.append(qso_cati)
#                 n += 1
#             except:
#                 print('healpix '+ssd +' failed')
#         #if n > 3:
#         #    break
#     qso_cat = vstack(qsocats,metadata_conflicts='silent')
#     common.write_LSS(qso_cat, outf)   



