'''
ENVIRONMENT
=============
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
PYTHONPATH=$PYTHONPATH:$HOME/LSS/py 


EXAMPLE USE
===========

srun -n 64 -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/apply_fnl_blinding.py --type LRG --basedir_out /global/cfs/cdirs/desi/survey/catalogs/DA2/LSS/loa-v1 --version v1.1 
(using -n 64 slightly slows down fNL blinding but greatly speeds up writing of updated randoms)

'''

import sys
import os
import logging
import shutil
import unittest
from datetime import datetime
import json
import numpy as np
#from numpy.random import MT19937
#from numpy.random import RandomState, SeedSequence
from numpy.random import random
import fitsio
import glob
import argparse
from astropy.io import fits
from astropy.table import Table,join,unique,vstack
from matplotlib import pyplot as plt

import LSS.main.cattools as ct
import LSS.blinding_tools as blind
from LSS.tabulated_cosmo import TabulatedDESI

import LSS.recon_tools as rectools
from LSS.cosmodesi_io_tools import catalog_fn
import LSS.common_tools as common

from cosmoprimo.fiducial import DESI
from cosmoprimo.utils import DistanceToRedshift
from cosmoprimo import Cosmology


if os.environ['NERSC_HOST'] == 'cori':
    scratch = 'CSCRATCH'
elif os.environ['NERSC_HOST'] == 'perlmutter':
    scratch = 'PSCRATCH'
else:
    print('NERSC_HOST is not cori or permutter but is '+os.environ['NERSC_HOST'])
    sys.exit('NERSC_HOST not known (code only works on NERSC), not proceeding')




# to remove jax warning (from cosmoprimo)
logging.getLogger("jax._src.lib.xla_bridge").addFilter(logging.Filter("No GPU/TPU found, falling back to CPU."))

logname = 'LSS_blinding'
logger = logging.getLogger(logname)
logger.setLevel(logging.INFO)

# create console handler and set level to debug
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)

# create formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# add formatter to ch
ch.setFormatter(formatter)

# add ch to logger
logger.addHandler(ch)


parser = argparse.ArgumentParser()
parser.add_argument("--type", help="tracer type to be selected")
parser.add_argument("--basedir_in", help="base directory for input, default is location for official catalogs", default='/dvs_ro/cfs/cdirs/desi/survey/catalogs/')
parser.add_argument("--basedir_out", help="base directory for output, default is C(P)SCRATCH", default=os.environ[scratch])
parser.add_argument("--KP", help="extra directory for KP catalogs", default='fNL')
parser.add_argument("--version", help="catalog version", default='test')
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA", default='DA2')
parser.add_argument("--verspec", help="version for redshifts", default='loa-v1')
parser.add_argument("--notqso", help="if y, do not include any qso targets", default='n')
parser.add_argument("--minr", help="minimum number for random files", default=0, type=int)# use 1 for abacus mocks
parser.add_argument("--maxr", help="maximum for random files, default is 1", default=18, type=int) # use 2 for abacus mocks

parser.add_argument("--get_par_mode", help="how to get the row of the file with w0/wa values", choices=['random', 'from_file','specified'],default='from_file')

parser.add_argument("--fnlblind", help="if y, do the fnl blinding", default='n')

parser.add_argument("--fiducial_f", help="fiducial value for f", default=0.8)

#relevant if args.get_par_mode is specified
parser.add_argument("--specified_w0",
                    help="Specify a blind w0 value",
                    default=None)
parser.add_argument("--specified_wa",
                    help="Specify a blind wa value ",
                    default=None)
parser.add_argument("--specified_fnl",
                    help="Specify a blind fnl value ",
                    default=None)


parser.add_argument("--visnz",help="whether to look at the original, blinded, and weighted n(z)",default='n')
parser.add_argument("--useMPI",help="whether to try to use MPI or not",default='y')


args = parser.parse_args()

mpicomm = None
if args.useMPI == 'y':
    try:
        from pypower import mpi
        mpicomm = mpi.COMM_WORLD
        #mpicomm = pyrecon.mpi.COMM_WORLD  # MPI version
    except AttributeError:
        mpicomm = None  # non-MPI version
        sys.exit('Not in MPI mode. The fNL blinding requires MPI, the script will exit before attempting fNL blinding')
        #sys.exit('The following script need to be run with the MPI version of pyrecon. Please use module swap pyrecon:mpi')

if mpicomm is None:
    sys.exit('NOT using MPI. Required, exiting')

root = mpicomm is None or mpicomm.rank == 0


if root: common.printlog(args,logger)

type = args.type
version = args.version
specrel = args.verspec

notqso = 'notqso' if (args.notqso == 'y') else ''
if root: common.printlog('blinding catalogs for tracer type ' + type + notqso,logger)

prog = 'BRIGHT' if (type[:3] == 'BGS' or type == 'bright' or type == 'MWS_ANY') else 'DARK'
progl = prog.lower()

maindir = args.basedir_in +'/'+args.survey+'/LSS/'
if 'mock' not in args.verspec:
    

    ldirspec = maindir+specrel+'/'

    dirin = ldirspec+'LSScats/'+version+'/'
    LSSdir = ldirspec+'LSScats/'
    
    nzmd = 'data'
    dirout = args.basedir_out + '/LSScats/' + version + '/'+args.KP+'/blinded/'
    dirfid = args.basedir_out + '/LSScats/' + version + '/'+args.KP+'/' #where any unblinded catalogs would be

elif 'AbacusSummit_v4_1' in args.verspec: #e.g., use 'mocks/FirstGenMocks/AbacusSummit/Y1/mock1' to get the 1st mock with fiberassign
                                          # VERSPEC="mocks/SecondGenMocks/AbacusSummit_v4_1/mock<number>"
    # dirin = args.basedir_in +'/'+args.survey+'/'+args.verspec+'/LSScats/'+version+'/'
    dirin = args.basedir_in +'/'+args.survey+'/'+args.verspec+'/'

    # LSSdir = args.basedir_in +'/'+args.survey+'/'+args.verspec+'/LSScats/'
    dchi2=None
    tsnrcut=0
    if "FirstGen" in args.verspec:
        randens = 10460.
    nzmd = 'mock'
    dirout = args.basedir_out + '/blinded_mock/'

else:
    sys.exit('verspec '+args.verspec+' not supported')

# dirout = args.basedir_out + '/LSScats/' + version + '/BAO/blinded/'

def mkdir(dirname):
    """Try to create ``dirname`` and catch :class:`OSError`."""
    try:
        os.makedirs(dirname)  # MPI...
    except OSError:
        return

if root:    
    mkdir(dirout)


tp2z = {'LRG': 0.8, 'ELG': 1.1, 'QSO': 1.6,'BGS':0.25}
tp2bias = {'LRG': 2., 'ELG': 1.3, 'QSO': 2.3,'BGS':1.8}

regl = ['_S', '_N']
gcl = ['_SGC', '_NGC']

if root:
    ztp = tp2z[args.type[:3]]
    bias = tp2bias[args.type[:3]]

    w0wa = np.loadtxt('/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/w0wa_initvalues_zeffcombined_1000realisations.txt')

    if args.get_par_mode == 'specified':
        [w0_blind, wa_blind] = [float(args.specified_w0),float(args.specified_wa)]
        if w0_blind is None or wa_blind is None:
            sys.exit('you must provide arguments for --specified_w0 and --specified_wa in the specified get_par_mode')
    
    if args.get_par_mode == 'random':
        #if args.type != 'LRG':
        #    sys.exit('Only do LRG in random mode, read from LRG file for other tracers')
        ind = int(random() * 1000)
        [w0_blind, wa_blind] = w0wa[ind]

    if args.get_par_mode == 'from_file':
        #fn = LSSdir + 'filerow.txt'
        fn = maindir + 'filerow.txt' #location switched for DR2 so that same blinding can be applied to different spectroscopic reductions
        if not os.path.isfile(fn):
            common.printlog('will write '+fn,logger)
            ind_samp = int(random()*1000)
            fo = open(fn.replace('dvs_ro','global'),'w')
            fo.write(str(ind_samp)+'\n')
            fo.close()
        common.printlog('reading from '+fn,logger)
        ind = int(np.loadtxt(fn))    
        [w0_blind,wa_blind] = w0wa[ind]

    #choose f_shift to compensate shift in monopole amplitude
    cosmo_fid = DESI()
    cosmo_shift = cosmo_fid.clone(w0_fld=w0_blind, wa_fld=wa_blind)

    DM_fid = cosmo_fid.comoving_angular_distance(ztp)
    DH_fid = 1. / cosmo_fid.hubble_function(ztp)

    DM_shift = cosmo_shift.comoving_angular_distance(ztp)
    DH_shift = 1. / cosmo_shift.hubble_function(ztp)

    vol_fac =  (DM_shift**2 * DH_shift) / (DM_fid**2 * DH_fid)

    #a, b, c for quadratic formula
    a = 0.2 / bias**2
    b = 2 / (3 * bias)
    c = 1 - (1 + 0.2 * (args.fiducial_f / bias)**2. + 2/3 * args.fiducial_f / bias) / vol_fac

    f_shift = (-b + np.sqrt(b**2. - 4.*a*c))/(2*a)
    dfper = (f_shift - args.fiducial_f)/args.fiducial_f
    maxfper = 0.1
    if abs(dfper) > maxfper:
        dfper = maxfper*dfper/abs(dfper)
        f_shift = (1+dfper)*args.fiducial_f
    fgrowth_blind = f_shift

 
    common.printlog('doing fNL blinding',logger)
from mockfactory.blinding import get_cosmo_blind, CutskyCatalogBlinding
#logger = logging.getLogger('recon')
if root:
    f_blind = fgrowth_blind
    if args.get_par_mode == 'specified':
        fnl_blind = args.specified_fnl
        if fnl_blind is None:
            sys.exit('you must provide arguments for --specified_fnl  in the specified get_par_mode')
        fnl_blind = float(fnl_blind )
        print('fnl value is '+str(fnl_blind))
    else:
        # generate blinding value from the choosen index above
        np.random.seed(ind)
        fnl_blind = np.random.uniform(low=-15, high=15, size=1)[0]
    
if not root:
    w0_blind, wa_blind, f_blind, fnl_blind = None, None, None, None
w0_blind = mpicomm.bcast(w0_blind, root=0)
wa_blind = mpicomm.bcast(wa_blind, root=0)
f_blind = mpicomm.bcast(f_blind, root=0)
fnl_blind = mpicomm.bcast(fnl_blind, root=0)

# collect effective redshift and bias for the considered tracer
zeff = tp2z[args.type[:3]]
bias = tp2bias[args.type[:3]]

# build blinding cosmology
cosmo_blind = get_cosmo_blind('DESI', z=zeff)
#cosmo_blind.params['w0_fld'] = w0_blind
#cosmo_blind.params['wa_fld'] = wa_blind
#cosmo_blind._derived['f'] = f_blind
cosmo_blind._derived['fnl'] = fnl_blind   # on fixe la valeur pour de bon
blinding = CutskyCatalogBlinding(cosmo_fid='DESI', cosmo_blind=cosmo_blind, bias=bias, z=zeff, position_type='rdz', mpicomm=mpicomm, mpiroot=0)

# loop over the different region of the sky
#regions = ['N', 'S'] if args.reg_md == 'NS' else ['NGC', 'SGC']
regions = ['NGC', 'SGC']
for region in regions:
    # path of data and randoms:
    cat_dir = dirout
    if args.KP == 'fNL':
        cat_dir = dirfid
    catalog_kwargs = dict(tracer=args.type, region=region, ctype='clustering', nrandoms=(args.maxr - args.minr))
    data_fn = catalog_fn(**catalog_kwargs, cat_dir=cat_dir, name='data')
    data_outfn = catalog_fn(**catalog_kwargs, cat_dir=dirout, name='data')
    randoms_fn = catalog_fn(**catalog_kwargs, cat_dir=cat_dir, name='randoms')
    if np.ndim(randoms_fn) == 0: randoms_fn = [randoms_fn]

    data_positions, data_weights = None, None
    randoms_positions, randoms_weights = None, None
    if root:
        logger.info('Loading {}.'.format(data_fn))
        data = Table.read(data_fn)
        data_positions, data_weights = [np.array(data['RA'], dtype='float64'), np.array(data['DEC'], dtype='float64'), np.array(data['Z'], dtype='float64')], data['WEIGHT']

        logger.info('Loading {}'.format(randoms_fn))
        randoms = vstack([Table(fitsio.read(fn.replace('global','dvs_ro'))) for fn in randoms_fn])
        
        randoms_positions, randoms_weights = [np.array(randoms['RA'], dtype='float64'), np.array(randoms['DEC'], dtype='float64'), np.array(randoms['Z'], dtype='float64')], randoms['WEIGHT']
        del randoms
    # add fnl blinding weight to the data weight
    new_data_weights = blinding.png(data_positions, data_weights=data_weights,
                                    randoms_positions=randoms_positions, randoms_weights=randoms_weights,
                                    method='data_weights', shotnoise_correction=True)

    
    
    # overwrite the data!
    if root:        
        del randoms_positions
        del randoms_weights
        fnl_blind_weights = new_data_weights / data['WEIGHT']
        data['WEIGHT'] = new_data_weights
        data['WEIGHT_BLIND'] = fnl_blind_weights
        data['WEIGHT_COMP'] *= fnl_blind_weights
        common.write_LSS_scratchcp(data, data_outfn,logger=logger)

if root:
    dngc = fitsio.read(dirout+args.type+'_NGC_clustering.dat.fits')
    dsgc = fitsio.read(dirout+args.type+'_SGC_clustering.dat.fits')
    data = Table(np.concatenate([dngc,dsgc]))
    data.keep_columns(['TARGETID','WEIGHT_BLIND'])
    data.rename_column('TARGETID', 'TARGETID_DATA')

    del dngc
    del dsgc
    #now, adjust the random weight
    def _parfun(rannum):
        for region in regions:
            common.printlog('doing '+region+' random '+str(rannum),logger)
            ranf_in = dirfid + args.type + notqso +'_'+ region + '_' + str(rannum) + '_clustering.ran.fits'
            ranf_out = dirout + args.type + notqso +'_'+ region + '_' + str(rannum) + '_clustering.ran.fits'
            rans = fitsio.read(ranf_in)
            rans = join(rans,data,keys=['TARGETID_DATA'])
            rans['WEIGHT'] *= rans['WEIGHT_BLIND']
            common.write_LSS_scratchcp(rans,ranf_out,logger=logger)

    
    nran = args.maxr-args.minr
    inds = np.arange(args.minr,args.maxr)

    from multiprocessing import Pool
    with Pool() as pool:
        res = pool.map(_parfun, inds)


    common.printlog('done with fNL blinding, done with script',logger)
