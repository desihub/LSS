import os
import argparse
import logging
import numpy as np
from astropy.table import Table, vstack
from matplotlib import pyplot as plt

# To install pycorr: python -m pip install git+https://github.com/adematti/pycorr
from pycorr import TwoPointCorrelationFunction, TwoPointEstimator, utils, project_to_multipoles, project_to_wp, setup_logging,KMeansSubsampler
from LSS.tabulated_cosmo import TabulatedDESI
cosmo = TabulatedDESI()
distance = cosmo.comoving_radial_distance

parser = argparse.ArgumentParser()
parser.add_argument("--type", help="tracer type to be selected")
parser.add_argument("--basedir", help="where to find catalogs",default='/global/cfs/cdirs/desi/survey/catalogs')
parser.add_argument("--version", help="catalog version",default='test')
parser.add_argument("--verspec",help="version for redshifts",default='everest')
parser.add_argument("--survey",help="e.g., SV3 or main",default='SV3')
parser.add_argument("--nran",help="number of random files to combine together (1-18 available)",default=4,type=int)
parser.add_argument("--weight_type",help="types of weights to use; use default_angular_bitwise for PIP with angular upweighting; default just uses WEIGHT column",default='default')
parser.add_argument("--bintype",help="log or lin",default='lin')
parser.add_argument("--njack",help="number of jack-knife subsamples",default=20)
parser.add_argument("--nthreads",help="number of threads for parallel comp",default=64,type=int)
parser.add_argument("--vis",help="set to y to plot each xi ",default='n')
parser.add_argument("--onlyfull",help="only do min max z",default='n')
parser.add_argument("--outdir",help='base directory for output',default=None)

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

if args.bintype == 'log':
    bine = np.logspace(-1.5, 2.2, 80)
if args.bintype == 'lin':
    bine = np.linspace(1e-4, 200, 201)

dirxi = os.environ['CSCRATCH']+'/'+survey+'xi/'
if args.outdir is not None:
    dirxi = args.outdir

if not os.path.exists(dirxi):
    os.mkdir(dirxi)
    print('made '+dirxi) 

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
    if 'safez' in ttype:
        zl = [0.9,1.48]
        ttype = ttype.strip('safez')
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

if args.onlyfull == 'y':
    zl = [zl[0],zl[-1]]
    print(zl)
wa = ''
if survey in ['main', 'DA02']:
    wa = 'zdone'

if 'completeness_only' in weight_type and 'bitwise' in weight_type:
    sys.exit('inconsistent choices were put into weight_type, not proceeding!')

def sel_reg(ra,dec,reg):
    wra = (ra > 100-dec)
    wra &= (ra < 280 +dec)
    if reg == 'DN':
        w = dec < 32.375
        w &= wra
    if reg == 'DS':
        w = ~wra
        w &= dec > -25
    return w        


def compute_correlation_function(mode, tracer='LRG', region='_N', nrandoms=4, zlim=(0., np.inf), weight_type=None, nthreads=8, dtype='f8', wang=None,fnroot=''):
    if ttype == 'ELGrec' or ttype == 'LRGrec':
        data_fn = os.path.join(dirname, tracer+wa+ region+'_clustering_'+args.rectype+args.convention+'.dat.fits')
        data = Table.read(data_fn)

        shifted_fn = os.path.join(dirname, tracer+wa+ region+'_clustering_'+args.rectype+args.convention+'.ran.fits')
        shifted = Table.read(shifted_fn)

        randoms_fn = [os.path.join(dirname, '{}{}_{:d}_clustering.ran.fits'.format(tracer+wa, region, iran)) for iran in range(nrandoms)]
        randoms = vstack([Table.read(fn) for fn in randoms_fn])
    else:
        data_fn = os.path.join(dirname, '{}{}_clustering.dat.fits'.format(tracer+wa, region))
        data = Table.read(data_fn)

        shifted = None

        randoms_fn = [os.path.join(dirname, '{}{}_{:d}_clustering.ran.fits'.format(tracer+wa, region, iran)) for iran in range(nrandoms)]
        randoms = vstack([Table.read(fn) for fn in randoms_fn])
  
    corrmode = mode
    if mode == 'wp':
        corrmode = 'rppi'
    if mode == 'multi':
        corrmode = 'smu'
    if corrmode == 'smu':
        edges = (bine, np.linspace(-1., 1., 201)) #s is input edges and mu evenly spaced between 0 and 1
    if corrmode == 'rppi':
        edges = (bine, bine) #transverse and radial separations are  coded to be the same here
        if mode == 'wp':
            edges = (bins,np.linspace(0,40.,41)) #if you want wp, only go out to pi = 40; consider setting pi_max as argument
    def get_positions_weights(catalog, name='data'):
        mask = (catalog['Z'] >= zlim[0]) & (catalog['Z'] < zlim[1])
        print('Using {:d} rows for {}'.format(mask.sum(),name))
        positions = [catalog['RA'][mask].astype(float), catalog['DEC'][mask].astype(float), distance(catalog['Z'][mask]).astype(float)]
        #if weight_type is None:
        #    weights = None
        #else:
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
        if name == 'randoms':
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
    subsampler = KMeansSubsampler(mode='angular', positions=data_positions, nsamples=int(args.njack), nside=512, random_state=42, position_type='rdd')
    labels = subsampler.label(randoms_positions)
    data_samples = subsampler.label(data_positions)
    randoms_samples = subsampler.label(randoms_positions)

    if shifted is not None:
        shifted_positions, shifted_weights = get_positions_weights(shifted, name='shifted')

    kwargs = {}
    if 'angular' in weight_type and wang is None:
        
        if tracer != 'LRG_main':
            data_fn = os.path.join(dirname, '{}_full.dat.fits'.format(tracer))
            randoms_fn = [os.path.join(dirname, '{}_{:d}_full.ran.fits'.format(tracer, iran)) for iran in range(nrandoms)]
        else:
            data_fn = os.path.join(dirname, '{}_full.dat.fits'.format(tracer[:3]))
            randoms_fn = [os.path.join(dirname, '{}_{:d}_full.ran.fits'.format(tracer[:3], iran)) for iran in range(nrandoms)]
        
        parent_data = Table.read(data_fn)
        parent_randoms = vstack([Table.read(fn) for fn in randoms_fn])
        
        def get_positions_weights(catalog, fibered=False):
            mask = np.ones(len(catalog),dtype='bool')
            if reg != '' and reg != '_DS' and reg != '_DN':
                mask &= catalog['PHOTSYS'] == region.strip('_')
            if reg == '_DS' or reg == '_DN':
                mask &= sel_reg(catalog['RA'],catalog['DEC'],region.strip('_'))
            
            if fibered: mask &= catalog['LOCATION_ASSIGNED']
            positions = [catalog['RA'][mask], catalog['DEC'][mask], catalog['DEC'][mask]]
            if fibered and 'bitwise' not in weight_type: weights = list(catalog['BITWEIGHTS'][mask].T)
            else: weights = np.ones_like(positions[0])
            return positions, weights
    
        fibered_data_positions, fibered_data_weights = get_positions_weights(parent_data, fibered=True)
        parent_data_positions, parent_data_weights = get_positions_weights(parent_data)
        parent_randoms_positions, parent_randoms_weights = get_positions_weights(parent_randoms)
        
        tedges = np.logspace(-3.5, 0.5, 31)
        # First D1D2_parent/D1D2_PIP angular weight
        wangD1D2 = TwoPointCorrelationFunction('theta', tedges, data_positions1=fibered_data_positions, data_weights1=fibered_data_weights,
                                               randoms_positions1=parent_data_positions, randoms_weights1=parent_data_weights,
                                               estimator='weight', engine='corrfunc', position_type='rdd', nthreads=nthreads, dtype=dtype)

        # First D1R2_parent/D1R2_IIP angular weight
        # Input bitwise weights are automatically turned into IIP
        wangD1R2 = TwoPointCorrelationFunction('theta', tedges, data_positions1=fibered_data_positions, data_weights1=fibered_data_weights,
                                               data_positions2=parent_randoms_positions, data_weights2=parent_randoms_weights,
                                               randoms_positions1=parent_data_positions, randoms_weights1=parent_data_weights,
                                               randoms_positions2=parent_randoms_positions, randoms_weights2=parent_randoms_weights,
                                               estimator='weight', engine='corrfunc', position_type='rdd', nthreads=nthreads, dtype=dtype)
        wang = {}
        wang['D1D2_twopoint_weights'] = wangD1D2
        wang['D1R2_twopoint_weights'] = wangD1R2
    
    kwargs.update(wang or {})

    result = TwoPointCorrelationFunction(corrmode, edges, data_positions1=data_positions, data_weights1=data_weights,
                                         randoms_positions1=randoms_positions, randoms_weights1=randoms_weights,
                                         shifted_positions1=shifted_positions, shifted_weights1=shifted_weights,
                                         data_samples1=data_samples, randoms_samples1=randoms_samples,
                                         engine='corrfunc', position_type='rdd', nthreads=nthreads, dtype=dtype, **kwargs)
    #save paircounts
    fn = dirxi+'paircounts_'+fnroot+'.npy'
    result.save(fn)
    return result, wang

ranwt1=False

regl = ['_N','_S','']

tcorr = ttype
tw = ttype
if survey in ['main','DA02']:
    regl = ['_DN','_DS','_N','_S','']
    if ttype == 'LRGrec':
        regl = ['_DN','_N']
        tcorr = 'LRG'
        tw = 'LRG'+args.rectype+args.convention
    if ttype == 'ELGrec':
        regl = ['_DN','_N']
        tcorr = 'ELG'
        tw = 'ELG'+args.rectype+args.convention
        

nzr = len(zl)

bsl = [1,4,5,10]
ells = (0, 2, 4)
nells = len(ells)

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
        #(sep, xiell), wang = compute_correlation_function('multi', bs, tracer=tcorr, region=reg, nrandoms=args.nran, zlim=(zmin,zmax), weight_type=weight_type,nthreads=args.nthreads)
        fnroot = tw+survey+reg+'_'+str(zmin)+str(zmax)+version+'_'+weight_type+args.bintype
        pfn = dirxi+'paircounts_'+fnroot+'.npy'
        result,wang = compute_correlation_function('multi', tracer=tcorr, region=reg, nrandoms=args.nran, zlim=(zmin,zmax), weight_type=weight_type,nthreads=args.nthreads,fnroot=fnroot)
        for bs in bsl:
            result = TwoPointEstimator.load(pfn)
            result.rebin((bs, 1))
            sep,xiell,cov = project_to_multipoles(result,ells=ells)#, wang
            std = np.array_split(np.diag(cov)**0.5, nells)
            fo = open(dirxi+'xi024'+fnroot+str(bs)+'.dat','w')
            for i in range(0,len(sep)):
                fo.write(str(sep[i])+' '+str(xiell[0][i])+' '+str(xiell[1][i])+' '+str(xiell[2][i])+' '+str(std[0][i])+' '+str(std[1][i])+' '+str(std[2][i])+'\n')
            fo.close()
            if args.vis == 'y':
                if args.bintype == 'log':
                    plt.loglog(sep,xiell[0])
                if args.bintype == 'lin':
                    plt.plot(sep,sep**2.*xiell[0])
                plt.title(ttype+' '+str(zmin)+'<z<'+str(zmax)+' in '+reg)
                plt.show()    
