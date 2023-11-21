'''
Documentation needs to be updated
EXAMPLE USE
===========



GENERAL NOTES
=============


NOTES FOR TESTING AND VALIDATION
================================

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
from LSS.globals import main
import LSS.blinding_tools as blind
from LSS.tabulated_cosmo import TabulatedDESI

import LSS.recon_tools as rectools
from LSS.cosmodesi_io_tools import catalog_fn
import LSS.common_tools as common

import pyrecon
from pyrecon import MultiGridReconstruction, IterativeFFTReconstruction, IterativeFFTParticleReconstruction, utils, setup_logging


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


parser = argparse.ArgumentParser()
parser.add_argument("--type", help="tracer type to be selected")
parser.add_argument("--basedir_in", help="base directory for input, default is location for official catalogs", default='/dvs_ro/cfs/cdirs/desi/survey/catalogs/')
parser.add_argument("--basedir_out", help="base directory for output, default is C(P)SCRATCH", default=os.environ[scratch])
parser.add_argument("--version", help="catalog version", default='test')
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA", default='Y1')
parser.add_argument("--verspec", help="version for redshifts", default='iron')
parser.add_argument("--notqso", help="if y, do not include any qso targets", default='n')
parser.add_argument("--use_map_veto",help="string to add on the end of full file reflecting if hp maps were used to cut",default='_HPmapcut')
#parser.add_argument("--reg_md", help="whether to run on split N/S or NGC/SGC", default='GC')

#parser.add_argument("--split_GC", help="whether to make the split NGC/SGC", default='y')

parser.add_argument("--get_par_mode", help="how to get the row of the file with w0/wa values", choices=['random', 'from_file','specified'],default='from_file')

parser.add_argument("--baoblind", help="if y, do the bao blinding shift", default='n')
parser.add_argument("--compmd", help="whether the extra completeness gets added to data or random", choices=['dat','ran'],default='ran')

parser.add_argument("--mkclusdat", help="if y, make the clustering data files after the BAO blinding (needed for RSD blinding)", default='n')
parser.add_argument("--wsyscol", help="column name to use for WEIGHT_SYS", default='WEIGHT_SN')
parser.add_argument("--mkclusran", help="if y, make the clustering random files after the BAO blinding (needed for RSD blinding)", default='n')
parser.add_argument("--minr", help="minimum number for random files", default=0, type=int)# use 1 for abacus mocks
parser.add_argument("--maxr", help="maximum for random files, default is 1", default=1, type=int) # use 2 for abacus mocks
parser.add_argument("--dorecon", help="if y, run the recon needed for RSD blinding", default='n')
parser.add_argument("--rsdblind", help="if y, do the RSD blinding shift", default='n')
parser.add_argument("--fnlblind", help="if y, do the fnl blinding", default='n')
parser.add_argument("--resamp", help="resample the randoms to make sure all is consistent with how weights changed", default='n')
parser.add_argument("--getFKP", help="calculate n(z) and FKP weights on final clustering catalogs", default='n')

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
        mpicomm = pyrecon.mpi.COMM_WORLD  # MPI version
    except AttributeError:
        mpicomm = None  # non-MPI version
        print('Not in MPI mode. The fNL blinding requires MPI, the script will exit before attempting fNL blinding')
        #sys.exit('The following script need to be run with the MPI version of pyrecon. Please use module swap pyrecon:mpi')

if mpicomm is None:
    print('NOT using MPI. If you specified a number of processes, e.g. "srun ... -n 32", greater than 1, things will not work well')

root = mpicomm is None or mpicomm.rank == 0


if root: print(args)

type = args.type
version = args.version
specrel = args.verspec

notqso = 'notqso' if (args.notqso == 'y') else ''
if root: print('blinding catalogs for tracer type ' + type + notqso)

prog = 'BRIGHT' if (type[:3] == 'BGS' or type == 'bright' or type == 'MWS_ANY') else 'DARK'
progl = prog.lower()

mainp = main(args.type,survey='Y1')
zmin = mainp.zmin
zmax = mainp.zmax
tsnrcol = mainp.tsnrcol

#share basedir location '/global/cfs/cdirs/desi/survey/catalogs'
if 'mock' not in args.verspec:
    maindir = args.basedir_in +'/'+args.survey+'/LSS/'

    ldirspec = maindir+specrel+'/'

    dirin = ldirspec+'LSScats/'+version+'/'
    dirin_blind = ldirspec+'LSScats/'+version+'/blinded/'
    LSSdir = ldirspec+'LSScats/'
    tsnrcut = mainp.tsnrcut
    dchi2 = mainp.dchi2
    randens = 2500.
    nzmd = 'data'
elif 'Y1/mock' in args.verspec: #e.g., use 'mocks/FirstGenMocks/AbacusSummit/Y1/mock1' to get the 1st mock with fiberassign
    dirin = args.basedir_in +'/'+args.survey+'/'+args.verspec+'/LSScats/'+version+'/'
    LSSdir = args.basedir_in +'/'+args.survey+'/'+args.verspec+'/LSScats/'
    dchi2=None
    tsnrcut=0
    randens = 10460.
    nzmd = 'mock'

else:
    sys.exit('verspec '+args.verspec+' not supported')

dirout = args.basedir_out + '/LSScats/' + version + '/doubleblinded/'

def mkdir(dirname):
    """Try to create ``dirname`` and catch :class:`OSError`."""
    try:
        os.makedirs(dirname)  # MPI...
    except OSError:
        return
    
mkdir(dirout)


#if root and (not os.path.exists(dirout)):
#    os.makedirs(dirout)
#    print('made ' + dirout)


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
        fn = LSSdir + 'filerow.txt'
        if not os.path.isfile(fn):
            ind_samp = int(random()*1000)
            fo = open(fn,'w')
            fo.write(str(ind_samp)+'\n')
            fo.close()
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

    fb_in = dirin + type + notqso
    fbr_in = fb_in
    if type == 'BGS_BRIGHT-21.5':
        fbr_in = dirin +'BGS_BRIGHT'
    #fcr_in = fbr_in +'{}_0_clustering.ran.fits'
    fcd_in = dirin_blind+ type + notqso+ '{}_clustering.dat.fits'
    print('input file is '+fcd_in)
    nzf_in = dirin + type + notqso + '_full_nz.txt'
    wo = 'y'
    if os.path.isfile(nzf_in):
        wo = 'n'
    if type[:3] == 'QSO':
        dz = 0.02
        #zmin = 0.8
        #zmax = 3.5
        P0 = 6000
    else:
        dz = 0.01
        #zmin = 0.01
        #zmax = 1.6

    if type[:3] == 'LRG':
        P0 = 10000
        #zmin = 0.4
        #zmax = 1.1
    if type[:3] == 'ELG':
        P0 = 4000
        #zmin = 0.6
        #zmax = 1.6
    if type[:3] == 'BGS':
        P0 = 7000
        #zmin = 0.1
        #zmax = 0.5



    if args.baoblind == 'y':
        
        data_clus_fn = dirin_blind + type + notqso+'_clustering.dat.fits'
        if os.path.isfile(data_clus_fn) == False:
            dl = []
            regl = ['_NGC','_SGC']
            for reg in regl:
                dl.append(fitsio.read(dirin_blind + type + notqso+reg+'_clustering.dat.fits'))  
                data = Table(np.concatenate(dl))
                if 'PHOTSYS' not in list(data.dtype.names):
                    data = common.addNS(data)
        else:
            data = Table(fitsio.read(data_clus_fn))
        
        #redo 'WEIGHT' column because of completeness factorization (applied back later if FKP option is chosen)
        cols = list(data.dtype.names)
        if 'WEIGHT_SYS' not in cols:
            if args.wsyscol is not None:
                fd['WEIGHT_SYS'] = np.copy(fd[args.wsyscol])
            else:
                print('did not find WEIGHT_SYS, putting it in as all 1')
                fd['WEIGHT_SYS'] = np.ones(len(fd))
        data['WEIGHT'] = data['WEIGHT_COMP']*data['WEIGHT_ZFAIL']*data['WEIGHT_SYS']
        #cl = gcl
        #for reg in cl:
        #    fnd = fcd_in.format(reg)
            # fndr = dirout + type + notqso + fcr_in.format(reg)
        #    data = Table(fitsio.read(fnd))
            # data_real = Table(fitsio.read(fndr))
        # fin = fitsio.read(fnd)
        # cols = list(fin.dtype.names)
        # nz_in = common.mknz_full(fcd_in, fcr_in, type[:3], bs=dz, zmin=zmin, zmax=zmax, write=wo, randens=randens, md=nzmd)

        # if 'WEIGHT_FKP' not in cols:
        #     print('adding FKP weights')
        #     common.addFKPfull(fcd_in, nz_in, type[:3], bs=dz, zmin=zmin, zmax=zmax, P0=P0, md=nzmd)

        
        # data = Table(fitsio.read(fcd_in))
        data['Z'] = np.clip(data['Z'],0.01,3.6)
        #outf = dirout + fnd.split('/')[-1]
        outf = dirout +  type + notqso+'_clustering.dat.fits'
        print('output going to '+outf)
        blind.apply_zshift_DE(data, outf, w0=w0_blind, wa=wa_blind, zcol='Z')

        #fb_out = dirout + type + notqso
        #fcd_out = fb_out + '_full.dat.fits'
        # nz_out = common.mknz_full(outf, fcr_in, type[:3], bs=dz, zmin=zmin, zmax=zmax, randens=randens, md=nzmd, zcol='Z')

        # ratio_nz = nz_in / nz_out

        fd = Table(fitsio.read(outf))
        zl = fd['Z']
        zr = zl > zmin
        zr &= zl < zmax
        fd = fd[zr]
        common.write_LSS(fd, outf)


    if args.visnz == 'y':
        print('min/max of weights for nz:')
        print(np.min(wl),np.max(wl))
        fdin = fitsio.read(fcd_in)
        a = plt.hist(fdin['Z'][gz],bins=100,range=(zmin,zmax),histtype='step',label='input')
        b = plt.hist(fd['Z'][gz],bins=100,range=(zmin,zmax),histtype='step',label='blinded')
        c = plt.hist(fd['Z'][gz],bins=100,range=(zmin,zmax),histtype='step',weights=fd['WEIGHT_SYS'][gz],label='blinded+reweight')
        plt.legend()
        plt.show()



    #if args.type == 'LRG':
    #    hdul = fits.open(fcd_out,mode='update')
    #    hdul['LSS'].header['FILEROW'] = ind
    #    hdul.close()
    #    hdtest = fitsio.read_header(dirout + 'LRG_full.dat.fits', ext='LSS')['FILEROW']
    #    if hdtest != ind:
    #        sys.exit('ERROR writing/reading row from blind file')


    if args.mkclusdat == 'y':
        ct.mkclusdat(dirout + type + notqso, tp=type, dchi2=dchi2, tsnrcut=tsnrcut, zmin=zmin, zmax=zmax,compmd=args.compmd)

    if args.mkclusran == 'y':
        rcols = ['Z', 'WEIGHT', 'WEIGHT_SYS', 'WEIGHT_COMP', 'WEIGHT_ZFAIL','WEIGHT_FKP','TARGETID_DATA','WEIGHT_SN']
        tsnrcol = 'TSNR2_ELG'
        if args.type[:3] == 'BGS':
            tsnrcol = 'TSNR2_BGS'
        #for rannum in range(args.minr, args.maxr):
        ranin = dirin + args.type + notqso + '_'
        if args.type == 'BGS_BRIGHT-21.5':
            ranin = dirin + 'BGS_BRIGHT' + notqso + '_'
        data_clus_fn = dirout + type + notqso+'_clustering.dat.fits'
        if os.path.isfile(data_clus_fn) == False:
            dl = []
            regl = ['_NGC','_SGC']
            for reg in regl:
                dl.append(fitsio.read(dirout + type + notqso+reg+'_clustering.dat.fits'))  
                dtot = np.concatenate(dl)
                if 'PHOTSYS' not in list(dtot.dtype.names):
                    dtot = common.addNS(Table(dtot))
                
                clus_arrays = [dtot]  
        else:
            clus_arrays = [fitsio.read(data_clus_fn)]
        #for reg in ['N','S']:
        #    clus_arrays.append(fitsio.read(dirout + type + notqso+'_'+reg+'_clustering.dat.fits'))
        
        def _parfun(rannum):
            ct.mkclusran(ranin, dirout + args.type + notqso + '_', rannum, rcols=rcols, tsnrcut=tsnrcut, tsnrcol=tsnrcol,clus_arrays=clus_arrays,use_map_veto=args.use_map_veto)#, ntilecut=ntile, ccut=ccut)
            #for clustering, make rannum start from 0
            if 'Y1/mock' in args.verspec:
                for reg in regl:
                    ranf = dirout + args.type + notqso + reg + '_' + str(rannum) + '_clustering.ran.fits'
                    ranfm = dirout + args.type + notqso + reg + '_' + str(rannum - 1) + '_clustering.ran.fits'
                    os.system('mv ' + ranf + ' ' + ranfm)
        nran = args.maxr-args.minr
        inds = np.arange(args.minr,args.maxr)
        if args.useMPI == 'y':
            from multiprocessing import Pool
            nproc = 9
            #nproc = nran*2
            with Pool(processes=nproc) as pool:
                res = pool.map(_parfun, inds)
        else:
            for ii in inds:
                _parfun(ii)
                #ct.mkclusran(ranin, dirout + args.type + notqso + '_', ii, rcols=rcols, tsnrcut=tsnrcut, tsnrcol=tsnrcol,clus_arrays=clus_arrays)
                print(ii,clus_arrays[0].dtype.names)
        #if args.split_GC == 'y':
        fb = dirout + args.type + notqso + '_'
        #ct.clusNStoGC(fb, args.maxr - args.minr)
        ct.splitclusGC(fb, args.maxr - args.minr)



#     if args.mkclusran == 'y':
#         rcols = ['Z', 'WEIGHT', 'WEIGHT_SYS', 'WEIGHT_COMP', 'WEIGHT_ZFAIL','WEIGHT_FKP','TARGETID_DATA','WEIGHT_SN']
#         tsnrcol = 'TSNR2_ELG'
#         if args.type[:3] == 'BGS':
#             tsnrcol = 'TSNR2_BGS'
#         #for rannum in range(args.minr, args.maxr):
#         ranin = dirin + args.type + notqso + '_'
#         if args.type == 'BGS_BRIGHT-21.5':
#             ranin = dirin + 'BGS_BRIGHT' + notqso + '_'
#         clus_arrays = []
#         for reg in ['N','S']:
#             clus_arrays.append(fitsio.read(dirout + type + notqso+'_'+reg+'_clustering.dat.fits'))
#         def _parfun(rannum):
#             ct.mkclusran(ranin, dirout + args.type + notqso + '_', rannum, rcols=rcols, tsnrcut=tsnrcut, tsnrcol=tsnrcol,clus_arrays=clus_arrays,use_map_veto=args.use_map_veto)#, ntilecut=ntile, ccut=ccut)
#             #for clustering, make rannum start from 0
#             if 'Y1/mock' in args.verspec:
#                 for reg in regl:
#                     ranf = dirout + args.type + notqso + reg + '_' + str(rannum) + '_clustering.ran.fits'
#                     ranfm = dirout + args.type + notqso + reg + '_' + str(rannum - 1) + '_clustering.ran.fits'
#                     os.system('mv ' + ranf + ' ' + ranfm)
#         nran = args.maxr-args.minr
#         inds = np.arange(args.minr,args.maxr)
#         if args.useMPI == 'y':
#             from multiprocessing import Pool
#             nproc = 9
#             #nproc = nran*2
#             with Pool(processes=nproc) as pool:
#                 res = pool.map(_parfun, inds)
#         else:
#             for ii in inds:
#                 _parfun(ii)
#                 #ct.mkclusran(ranin, dirout + args.type + notqso + '_', ii, rcols=rcols, tsnrcut=tsnrcut, tsnrcol=tsnrcol,clus_arrays=clus_arrays)
#                 print(ii,clus_arrays[0].dtype.names)
#         #if args.split_GC == 'y':
#         fb = dirout + args.type + notqso + '_'
#         ct.clusNStoGC(fb, args.maxr - args.minr)

    sys.stdout.flush()

if args.dorecon == 'y':
    distance = TabulatedDESI().comoving_radial_distance

    f, bias = rectools.get_f_bias(args.type)

    Reconstruction = IterativeFFTReconstruction#MultiGridReconstruction
    
    setup_logging() 

    #regions = ['N', 'S'] if args.reg_md == 'NS' else ['NGC', 'SGC']
    regions = ['NGC', 'SGC']
    for region in regions:
        catalog_kwargs = dict(tracer=args.type, region=region, ctype='clustering', nrandoms=(int(args.maxr) - int(args.minr)))
        data_fn = catalog_fn(**catalog_kwargs, cat_dir=dirout, name='data')
        randoms_fn = catalog_fn(**catalog_kwargs, cat_dir=dirout, name='randoms')
        #print(randoms_fn)
        data_rec_fn = catalog_fn(**catalog_kwargs, cat_dir=dirout, rec_type='IFFTrsd', name='data')
        randoms_rec_fn = catalog_fn(**catalog_kwargs, cat_dir=dirout, rec_type='IFFTrsd', name='randoms')
        rectools.run_reconstruction(Reconstruction, distance, data_fn, randoms_fn, data_rec_fn, randoms_rec_fn, f=f, bias=bias, convention='rsd', dtype='f8', zlim=(zmin, zmax), mpicomm=mpicomm)

if root:
    #re-sample redshift dependent columns from data
    nran = args.maxr-args.minr
    if args.resamp == 'y':
        regions = ['NGC', 'SGC']
        rcols = ['Z', 'WEIGHT', 'WEIGHT_SYS', 'WEIGHT_COMP', 'WEIGHT_ZFAIL','WEIGHT_FKP','TARGETID_DATA','WEIGHT_SN']
        for reg in regions:
            flin = dirout + args.type + notqso + '_'+reg    
            def _parfun(rannum):
                ct.clusran_resamp(flin,rannum,rcols=rcols,compmd=args.compmd)#, ntilecut=ntile, ccut=ccut)
            
            inds = np.arange(nran)
            from multiprocessing import Pool
            with Pool(processes=nran*2) as pool:
                res = pool.map(_parfun, inds)
if root and (args.rsdblind == 'y'):
    #if args.reg_md == 'NS':
    #    cl = regl
    #if args.reg_md == 'GC':
    cl = gcl
    for reg in cl:
        fnd = dirout + type + notqso + reg + '_clustering.dat.fits'
        fndr = dirout + type + notqso + reg + '_clustering.IFFTrsd.dat.fits'
        data = Table(fitsio.read(fnd))
        data_real = Table(fitsio.read(fndr))

        out_file = fnd
        blind.apply_zshift_RSD(data, data_real, out_file,
                               fgrowth_fid=args.fiducial_f,
                               fgrowth_blind=fgrowth_blind)#,
                               #comments=f"f_blind: {fgrowth_blind}, w0_blind: {w0_blind}, wa_blind: {wa_blind}")

if args.fnlblind == 'y':
    if mpicomm is None:
        sys.exit('fNL blinding requires MPI, exiting')
    from mockfactory.blinding import get_cosmo_blind, CutskyCatalogBlinding
    logger = logging.getLogger('recon')
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
    cosmo_blind = DESI()
    cosmo_blind._derived['fnl'] = fnl_blind  # only parameter used for PNG blinding
    blinding = CutskyCatalogBlinding(cosmo_fid='DESI', cosmo_blind=cosmo_blind, bias=bias, z=zeff, position_type='rdz', mpicomm=mpicomm, mpiroot=0)

    # loop over the different region of the sky
    #regions = ['N', 'S'] if args.reg_md == 'NS' else ['NGC', 'SGC']
    regions = ['NGC', 'SGC']
    for region in regions:
        # path of data and randoms:
        catalog_kwargs = dict(tracer=args.type, region=region, ctype='clustering', nrandoms=(args.maxr - args.minr))
        data_fn = catalog_fn(**catalog_kwargs, cat_dir=dirout, name='data')
        randoms_fn = catalog_fn(**catalog_kwargs, cat_dir=dirout, name='randoms')
        if np.ndim(randoms_fn) == 0: randoms_fn = [randoms_fn]

        data_positions, data_weights = None, None
        randoms_positions, randoms_weights = None, None
        if root:
            logger.info('Loading {}.'.format(data_fn))
            data = Table.read(data_fn)
            data_positions, data_weights = [np.array(data['RA'], dtype='float64'), np.array(data['DEC'], dtype='float64'), np.array(data['Z'], dtype='float64')], data['WEIGHT']

            logger.info('Loading {}'.format(randoms_fn))
            randoms = vstack([Table(fitsio.read(fn)) for fn in randoms_fn])
            randoms_positions, randoms_weights = [np.array(randoms['RA'], dtype='float64'), np.array(randoms['DEC'], dtype='float64'), np.array(randoms['Z'], dtype='float64')], randoms['WEIGHT']

        # add fnl blinding weight to the data weight
        new_data_weights = blinding.png(data_positions, data_weights=data_weights,
                                        randoms_positions=randoms_positions, randoms_weights=randoms_weights,
                                        method='data_weights', shotnoise_correction=True)

        # overwrite the data!
        if root:
            fnl_blind_weights = new_data_weights / data['WEIGHT']
            data['WEIGHT'] = new_data_weights
            data['WEIGHT_COMP'] = data['WEIGHT_COMP'] * fnl_blind_weights
            common.write_LSS(data, data_fn)

if root:
    #re-sample redshift dependent columns from data
    nran = args.maxr-args.minr
    if args.resamp == 'y':
        regions = ['NGC', 'SGC']
        rcols = ['Z', 'WEIGHT', 'WEIGHT_SYS', 'WEIGHT_COMP', 'WEIGHT_ZFAIL','WEIGHT_FKP','TARGETID_DATA','WEIGHT_SN']
        for reg in regions:
            flin = dirout + args.type + notqso + '_'+reg    
            def _parfun(rannum):
                ct.clusran_resamp(flin,rannum,rcols=rcols,compmd=args.compmd)#, ntilecut=ntile, ccut=ccut)
            
            inds = np.arange(nran)
            from multiprocessing import Pool
            with Pool(processes=nran*2) as pool:
                res = pool.map(_parfun, inds)

    
    if type[:3] == 'QSO':
        dz = 0.02
        zmin = 0.8
        zmax = 3.5
        P0 = 6000

    if type[:3] == 'LRG':
        P0 = 10000
        zmin = 0.4
        zmax = 1.1
    if type[:3] == 'ELG':
        P0 = 4000
        zmin = 0.8
        zmax = 1.6
    if type[:3] == 'BGS':
        P0 = 7000
        zmin = 0.1
        zmax = 0.4

    if args.getFKP == 'y':
        for reg in gcl:
            fb = dirout+args.type+reg
            fcr = fb+'_0_clustering.ran.fits'
            fcd = fb+'_clustering.dat.fits'
            fout = fb+'_nz.txt'
            common.mknz(fcd,fcr,fout,bs=dz,zmin=zmin,zmax=zmax,randens=randens)
            common.addnbar(fb,bs=dz,zmin=zmin,zmax=zmax,P0=P0,nran=nran)




    os.system('rm '+dirout+args.type+'*_S_*')
    os.system('rm '+dirout+args.type+'*_N_*')
    os.system('rm '+dirout+args.type+'*IFFT*')
    os.system('rm '+dirout+args.type+'*full*')
