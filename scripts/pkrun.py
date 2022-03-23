# To run: srun -n 64 python pkrun.py --type ELG...

import os
import argparse
import logging
import numpy as np
from astropy.table import Table, vstack
from matplotlib import pyplot as plt

# To install mockfactory: python -m pip install git+https://github.com/adematti/mockfactory
from pypower import CatalogFFTPower, PowerSpectrumStatistics, utils, setup_logging
from mockfactory import Catalog
from LSS.tabulated_cosmo import TabulatedDESI
cosmo = TabulatedDESI()
distance = cosmo.comoving_radial_distance

os.environ['NUMEXPR_MAX_THREADS'] = '8'
parser = argparse.ArgumentParser()
parser.add_argument("--type", help="tracer type to be selected")
parser.add_argument("--basedir", help="where to find catalogs",default='/global/cfs/cdirs/desi/survey/catalogs')
parser.add_argument("--version", help="catalog version",default='test')
parser.add_argument("--verspec",help="version for redshifts",default='everest')
parser.add_argument("--survey",help="e.g., SV3 or main",default='SV3')
parser.add_argument("--nran",help="number of random files to combine together (1-18 available)",default=10,type=int)
parser.add_argument("--weight_type",help="types of weights to use; use angular_bitwise for PIP; default just uses WEIGHT column",default='default')

#only relevant for reconstruction
parser.add_argument("--rectype",help="IFT or MG supported so far",default='IFT')
parser.add_argument("--convention",help="recsym or reciso supported so far",default='reciso')

setup_logging()
args = parser.parse_args()

ttype = args.type
basedir = args.basedir
version = args.version
specrel = args.verspec
survey = args.survey
nran = int(args.nran)
weight_type = args.weight_type

edges = {'min':0., 'step':0.001}

dirpk = os.environ['CSCRATCH']+'/'+survey+'pk/'

utils.mkdir(dirpk)
print('made '+dirpk) 

lssdir = basedir+'/'+survey+'/LSS/'+specrel+'/LSScats/'
dirname = lssdir + version
#dirout = svdir+'LSScats/'+version+'/'

zmask = ['']
minn = 0

subt = None

if ttype[:3] == 'LRG':
    zl = [0.4,0.6,0.8,1.1]


if ttype[:3] == 'ELG':# or type == 'ELG_HIP':
    #minn = 5
    zl = [0.8,1.1,1.5]
    #zmask = ['','_zmask']
    
    #zmin = 0.8
    #zmax = 1.6

if ttype == 'QSO':
    zl = [0.8,1.1,1.5,2.1]
    #zmin = 1.
    #zmax = 2.1

if ttype == 'QSOh':
    zl = [2.1,3.5]
    ttype = 'QSO'
    #zmin = 1.
    #zmax = 2.1

if ttype[:3] == 'BGS':
    #minn = 2
    zl = [0.1,0.3,0.5]
    #zmin = 0.1
    #zmax = 0.5 

if ttype[:3] == 'BGS' and ttype[-1] == 'l':
    #minn = 2
    zl = [0.1,0.3]
    ttype = ttype[:-1]
    #zmin = 0.1
    #zmax = 0.5 

if ttype[:3] == 'BGS' and ttype[-1] == 'h':
    #minn = 2
    zl = [0.3,0.5]
    ttype = ttype[:-1]
    #zmin = 0.1
    #zmax = 0.5 


wa = ''
if survey in ['main', 'DA02']:
    wa = 'zdone'

def compute_power_spectrum(edges, tracer='LRG', region='_N', nrandoms=4, zlim=(0., np.inf), weight_type=None, ells=(0, 2, 4), boxsize=5000., nmesh=1024, dtype='f4'):
    if ttype == 'ELGrec' or ttype == 'LRGrec':
        data_fn = os.path.join(dirname, tracer+wa+ region+'_clustering_'+args.rectype+args.convention+'.dat.fits')
        data = Catalog.load_fits(data_fn)

        shifted_fn = os.path.join(dirname, tracer+wa+ region+'_clustering_'+args.rectype+args.convention+'.ran.fits')
        shifted = Catalog.load_fits(shifted_fn)

        randoms_fn = [os.path.join(dirname, '{}{}_{:d}_clustering.ran.fits'.format(tracer+wa, region, iran)) for iran in range(nrandoms)]
        randoms = Catalog.concatenate(*(Catalog.load_fits(fn) for fn in randoms_fn), keep_order=False)
    else:
        data_fn = os.path.join(dirname, '{}{}_clustering.dat.fits'.format(tracer+wa, region))
        data = Catalog.load_fits(data_fn)

        shifted = None

        randoms_fn = [os.path.join(dirname, '{}{}_{:d}_clustering.ran.fits'.format(tracer+wa, region, iran)) for iran in range(nrandoms)]
        randoms = Catalog.concatenate(*(Catalog.load_fits(fn) for fn in randoms_fn), keep_order=False)
    
    def get_positions_weights(catalog, name='data'):
        mask = (catalog['Z'] >= zlim[0]) & (catalog['Z'] < zlim[1])
        positions = [catalog['RA'][mask], catalog['DEC'][mask], distance(catalog['Z'][mask])]
        nmask = catalog.mpicomm.allreduce(mask.sum())
        if catalog.mpicomm.rank == 0:
            catalog.log_info('Using {} rows for {}'.format(nmask, name))
        weights = np.ones_like(positions[0])
        if name == 'data':
            if 'zfail' in weight_type:
                weights *= catalog['WEIGHT_ZFAIL'][mask]
            if 'default' in weight_type and 'bitwise' not in weight_type:
                weights *= catalog['WEIGHT'][mask]
            if 'RF' in weight_type:
                weights *= catalog['WEIGHT_RF'][mask]*catalog['WEIGHT_COMP'][mask]
            if 'completeness_only' in weight_type:
                weights = catalog['WEIGHT_COMP'][mask]
            if 'FKP' in weight_type:
                weights *= catalog['WEIGHT_FKP'][mask]
            if 'bitwise' in weight_type:
                weights = list(catalog['BITWEIGHTS'][mask].T) + [weights]
        if name in ['randoms', 'shifted']:
            if 'default' in weight_type:
                weights *= catalog['WEIGHT'][mask]
            if 'RF' in weight_type:
                weights *= catalog['WEIGHT_RF'][mask]*catalog['WEIGHT_COMP'][mask]
            if 'zfail' in weight_type:
                weights *= catalog['WEIGHT_ZFAIL'][mask]
            if 'completeness_only' in weight_type:
                weights = catalog['WEIGHT_COMP'][mask]
            if 'FKP' in weight_type:
                weights *= catalog['WEIGHT_FKP'][mask]
        return positions, weights
    
    data_positions, data_weights = get_positions_weights(data, name='data')
    randoms_positions, randoms_weights = get_positions_weights(randoms, name='randoms')
    shifted_positions, shifted_weights = None, None
    if shifted is not None:
        shifted_positions, shifted_weights = get_positions_weights(shifted, name='shifted')

    result = CatalogFFTPower(data_positions1=data_positions, data_weights1=data_weights,
                             randoms_positions1=randoms_positions, randoms_weights1=randoms_weights,
                             shifted_positions1=shifted_positions, shifted_weights1=shifted_weights,
                             edges=edges, ells=ells, boxsize=boxsize, nmesh=nmesh, resampler='tsc', interlacing=2,
                             position_type='rdd', dtype=dtype)
    return result


regl = ['_N','_S']

tcorr = ttype
tw = ttype
if survey == 'main':
    regl = ['_DN','_DS','_N','_S']
    if ttype == 'LRGrec':
        regl = ['_DN','_N']
        tcorr = 'LRG'
        tw = 'LRG'+args.rectype+args.convention
    if ttype == 'ELGrec':
        regl = ['_DN','_N']
        tcorr = 'ELG'
        tw = 'ELG'+args.rectype+args.convention

bsl = [1, 5, 10]

nzr = len(zl)
if len(zl) == 2:
    nzr = len(zl)-1
for i in range(0,nzr):
    if i == len(zl)-1:
        zmin=zl[0]
        zmax=zl[-1]
    else:
        zmin = zl[i]
        zmax = zl[i+1]
    print(zmin,zmax)
    for reg in regl:
        print(reg)
        fnroot = tw+survey+reg+'_'+str(zmin)+str(zmax)+version+'_'+weight_type+'lin'
        pfn = os.path.join(dirpk, 'power_{}.npy'.format(fnroot))
        result = compute_power_spectrum(edges, tracer=tcorr, region=reg, nrandoms=args.nran, zlim=(zmin,zmax), weight_type=weight_type).poles
        result.save(pfn)
        for bs in bsl:
            #result = PowerSpectrumStatistics.load(pfn)
            fn_txt = os.path.join(dirpk, 'pk024{}{}.dat'.format(fnroot, bs))
            imax = (result.shape[0]//bs)*bs # cutting last bins to allow rebinning factor bs to divide shape
            result[:imax:bs].save_txt(fn_txt, complex=False)
