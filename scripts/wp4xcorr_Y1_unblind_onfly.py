#!/usr/bin/env python
# coding: utf-8

"""
Note
----
This is a modifed version of wp_Y1_unblind_only.py for using two tracers.
```
The script can be called with multiple processes as (e.g. on 2 nodes, 64 threads for each):
```
srun -n 2 python wp_Y1_unblind_onfly.py --nthreads 64 ...
```
Privilege nthreads over MPI processes.
"""

import os
import argparse
import logging

import numpy as np

from astropy.table import Table, vstack
from matplotlib import pyplot as plt

from pycorr import TwoPointCorrelationFunction, TwoPointEstimator, KMeansSubsampler, utils, setup_logging

from LSS.tabulated_cosmo import TabulatedDESI
import LSS.cosmodesi_io_tools as io
import LSS.main.cattools as ct
import LSS.common_tools as common



logger = logging.getLogger('xirunpc')




def compute_angular_weights(nthreads=8, dtype='f8', tracer='ELG', tracer2=None, mpicomm=None, mpiroot=None, **kwargs):

    autocorr = tracer2 is None
    catalog_kwargs = kwargs

    fibered_data_positions1, fibered_data_weights1, fibered_data_positions2, fibered_data_weights2 = None, None, None, None
    parent_data_positions1, parent_data_weights1, parent_data_positions2, parent_data_weights2 = None, None, None, None
    parent_randoms_positions1, parent_randoms_weights1, parent_randoms_positions2, parent_randoms_weights2 = None, None, None, None

    if mpicomm is None or mpicomm.rank == mpiroot:

        fibered_data = io.read_full_positions_weights(name='data', fibered=True, tracer=tracer, **catalog_kwargs)
        parent_data, parent_randoms = io.read_full_positions_weights(name=['data', 'randoms'], fibered=False, tracer=tracer, **catalog_kwargs)
        fibered_data_positions1, fibered_data_weights1 = io.concatenate_data_randoms(fibered_data)
        (parent_data_positions1, parent_data_weights1), (parent_randoms_positions1, parent_randoms_weights1) = io.concatenate_data_randoms(parent_data, parent_randoms, **catalog_kwargs)
        if not autocorr:
            fibered_data = io.read_full_positions_weights(name='data', fibered=True, tracer=tracer2, **catalog_kwargs)
            parent_data, parent_randoms = io.read_full_positions_weights(name=['data', 'randoms'], fibered=False, tracer=tracer2, **catalog_kwargs)
            fibered_data_positions2, fibered_data_weights2 = io.concatenate_data_randoms(fibered_data)
            (parent_data_positions2, parent_data_weights2), (parent_randoms_positions2, parent_randoms_weights2) = io.concatenate_data_randoms(parent_data, parent_randoms, **catalog_kwargs)

    tedges = np.logspace(-4., 0.5, 41)
    # First D1D2_parent/D1D2_PIP angular weight
    wangD1D2 = TwoPointCorrelationFunction('theta', tedges, data_positions1=fibered_data_positions1, data_weights1=fibered_data_weights1,
                                            data_positions2=fibered_data_positions2, data_weights2=fibered_data_weights2,
                                            randoms_positions1=parent_data_positions1, randoms_weights1=parent_data_weights1,
                                            randoms_positions2=parent_data_positions2, randoms_weights2=parent_data_weights2,
                                            estimator='weight', engine='corrfunc', position_type='rdd', nthreads=nthreads,
                                            dtype=dtype, mpicomm=mpicomm, mpiroot=mpiroot)

    # First D1R2_parent/D1R2_IIP angular weight
    # Input bitwise weights are automatically turned into IIP
    if autocorr:
         parent_randoms_positions2, parent_randoms_weights2 = parent_randoms_positions1, parent_randoms_weights1
    wangD1R2 = TwoPointCorrelationFunction('theta', tedges, data_positions1=fibered_data_positions1, data_weights1=fibered_data_weights1,
                                            data_positions2=parent_randoms_positions2, data_weights2=parent_randoms_weights2,
                                            randoms_positions1=parent_data_positions1, randoms_weights1=parent_data_weights1,
                                            randoms_positions2=parent_randoms_positions2, randoms_weights2=parent_randoms_weights2,
                                            estimator='weight', engine='corrfunc', position_type='rdd', nthreads=nthreads,
                                            dtype=dtype, mpicomm=mpicomm, mpiroot=mpiroot)
    wangR1D2 = None
    if not autocorr:
        wangR1D2 = TwoPointCorrelationFunction('theta', tedges, data_positions1=parent_randoms_positions1, data_weights1=parent_randoms_weights1,
                                               data_positions2=fibered_data_positions2, data_weights2=fibered_data_weights2,
                                               randoms_positions1=parent_randoms_positions1, randoms_weights1=parent_randoms_weights1,
                                               randoms_positions2=parent_data_positions2, randoms_weights2=parent_data_weights2,
                                               estimator='weight', engine='corrfunc', position_type='rdd', nthreads=nthreads,
                                               dtype=dtype, mpicomm=mpicomm, mpiroot=mpiroot)

    wang = {}
    wang['D1D2_twopoint_weights'] = wangD1D2
    wang['D1R2_twopoint_weights'] = wangD1R2
    wang['R1D2_twopoint_weights'] = wangR1D2

    return wang


def compute_correlation_function(corr_type, edges, distance, nthreads=8, dtype='f8', wang=None, split_randoms_above=30., weight_type='default', tracer='ELG', tracer2=None, recon_dir=None,rec_type=None, njack=120, option=None, mpicomm=None, mpiroot=None, cat_read=None, dat_cat=None, ran_cat=None, dat_cat2=None, ran_cat2=None, rpcut=None, **kwargs):

    autocorr = tracer2 is None
    catalog_kwargs = kwargs.copy()
    catalog_kwargs['weight_type'] = weight_type
    #catalog_kwargs['recon_dir'] = recon_dir
    with_shifted = rec_type is not None

    if 'angular' in weight_type and wang is None:
        wang = compute_angular_weights(nthreads=nthreads, dtype=dtype, weight_type=weight_type, tracer=tracer, tracer2=tracer2, mpicomm=mpicomm, mpiroot=mpiroot, **kwargs)

    data_positions1, data_weights1, data_samples1, data_positions2, data_weights2, data_samples2 = None, None, None, None, None, None
    randoms_positions1, randoms_weights1, randoms_samples1, randoms_positions2, randoms_weights2, randoms_samples2 = None, None, None, None, None, None
    shifted_positions1, shifted_weights1, shifted_samples1, shifted_positions2, shifted_weights2, shifted_samples2 = None, None, None, None, None, None
    jack_positions = None

    if mpicomm is None or mpicomm.rank == mpiroot:
        data, randoms = io.read_clustering_positions_weights(distance, name=['data', 'randoms'], recon_dir=recon_dir,rec_type=rec_type, tracer=tracer, option=option, cat_read=cat_read, dat_cat=dat_cat, ran_cat=ran_cat, **catalog_kwargs)

        if (with_shifted) & (cat_read == None):
            shifted = randoms  # above returned shifted randoms
            randoms = io.read_clustering_positions_weights(distance, name='randoms', rec_type=False, tracer=tracer, option=option, **catalog_kwargs)
        (data_positions1, data_weights1), (randoms_positions1, randoms_weights1) = io.concatenate_data_randoms(data, randoms, **catalog_kwargs)
        if with_shifted:
            shifted_positions1, shifted_weights1 = io.concatenate_data_randoms(data, shifted, **catalog_kwargs)[1]
        jack_positions = data_positions1

        if not autocorr:
            data, randoms = io.read_clustering_positions_weights(distance, name=['data', 'randoms'], recon_dir=recon_dir,rec_type=rec_type, tracer=tracer, option=option, cat_read=cat_read, dat_cat=dat_cat2, ran_cat=ran_cat2, **catalog_kwargs)
            if with_shifted & (cat_read == None):
                shifted = randoms
                randoms = io.read_clustering_positions_weights(distance, name='randoms', rec_type=False, tracer=tracer2, option=option, **catalog_kwargs)
            (data_positions2, data_weights2), (randoms_positions2, randoms_weights2) = io.concatenate_data_randoms(data, randoms, **catalog_kwargs)
            if with_shifted:
                shifted_positions2, shifted_weights2 = io.concatenate_data_randoms(data, shifted, **catalog_kwargs)[1]
            jack_positions = [np.concatenate([p1, p2], axis=0) for p1, p2 in zip(jack_positions, data_positions2)]

    if njack >= 2:
        subsampler = KMeansSubsampler('angular', positions=jack_positions, nsamples=njack, nside=512, random_state=42, position_type='rdd',
                                      dtype=dtype, mpicomm=mpicomm, mpiroot=mpiroot)

        if mpicomm is None or mpicomm.rank == mpiroot:
 
            data_samples1 = subsampler.label(data_positions1)
            randoms_samples1 = [subsampler.label(p) for p in randoms_positions1]
            if with_shifted:
                shifted_samples1 = [subsampler.label(p) for p in shifted_positions1]
            if not autocorr:
                data_samples2 = subsampler.label(data_positions2)
                randoms_samples2 = [subsampler.label(p) for p in randoms_positions2]
                if with_shifted:
                    shifted_samples2 = [subsampler.label(p) for p in shifted_positions2]

    kwargs = {}
    kwargs.update(wang or {})
    selection_attrs = None
    if rpcut is not None: selection_attrs = {'rp': (rpcut, np.inf)}
    randoms_kwargs = dict(randoms_positions1=randoms_positions1, randoms_weights1=randoms_weights1, randoms_samples1=randoms_samples1,
                          randoms_positions2=randoms_positions2, randoms_weights2=randoms_weights2, randoms_samples2=randoms_samples2,
                          shifted_positions1=shifted_positions1, shifted_weights1=shifted_weights1, shifted_samples1=shifted_samples1,
                          shifted_positions2=shifted_positions2, shifted_weights2=shifted_weights2, shifted_samples2=shifted_samples2)

    zedges = np.array(list(zip(edges[0][:-1], edges[0][1:])))
    mask = zedges[:,0] >= split_randoms_above
    zedges = [zedges[~mask], zedges[mask]]
    split_edges, split_randoms = [], []
    for ii, zedge in enumerate(zedges):
        if zedge.size:
            split_edges.append([np.append(zedge[:,0], zedge[-1,-1])] + list(edges[1:]))
            split_randoms.append(ii > 0)

    results = []
    if mpicomm is None:
        nran = len(randoms_positions1)
    else:
        nran = mpicomm.bcast(len(randoms_positions1) if mpicomm.rank == mpiroot else None, root=mpiroot)
    for i_split_randoms, edges in zip(split_randoms, split_edges):
        result = 0
        D1D2 = None
        for iran in range(nran if i_split_randoms else 1):
            tmp_randoms_kwargs = {}
            if i_split_randoms:
                # On scales above split_randoms_above, sum correlation function over multiple randoms
                for name, arrays in randoms_kwargs.items():
                    if arrays is None:
                        continue
                    else:
                        tmp_randoms_kwargs[name] = arrays[iran]
            else:
                # On scales below split_randoms_above, concatenate randoms
                for name, arrays in randoms_kwargs.items():
                    if arrays is None:
                        continue
                    elif isinstance(arrays[0], (tuple, list)):  # e.g., list of bitwise weights
                        array = [np.concatenate([arr[iarr] for arr in arrays], axis=0) for iarr in range(len(arrays[0]))]
                    else:
                        array = np.concatenate(arrays, axis=0)
                    tmp_randoms_kwargs[name] = array
            tmp = TwoPointCorrelationFunction(corr_type, edges, data_positions1=data_positions1, data_weights1=data_weights1, data_samples1=data_samples1,
                                              data_positions2=data_positions2, data_weights2=data_weights2, data_samples2=data_samples2,
                                              engine='corrfunc', position_type='rdd', nthreads=nthreads, dtype=dtype, **tmp_randoms_kwargs, **kwargs,
                                              D1D2=D1D2, mpicomm=mpicomm, mpiroot=mpiroot, selection_attrs=selection_attrs)
            D1D2 = tmp.D1D2
            result += tmp
        results.append(result)
        
        if tracer2 is not None:
            print('skipping random split')
            break              # currently, the random split (done in the second iteration of this loop), doesn't work for multiple tracers
    return results[0].concatenate_x(*results), wang


def get_edges(corr_type='smu', bin_type='lin'):

    if bin_type == 'log':
        sedges = np.geomspace(0.01, 100., 49)
    elif bin_type == 'lin':
        sedges = np.linspace(0., 200, 201)
    else:
        raise ValueError('bin_type must be one of ["log", "lin"]')
    if corr_type == 'smu':
        edges = (sedges, np.linspace(-1., 1., 201)) #s is input edges and mu evenly spaced between -1 and 1
    elif corr_type == 'rppi':
        rpmax = 80
        if bin_type == 'log':
            sedges = np.geomspace(0.01, rpmax, rpmax+1)
        elif bin_type == 'lin':
            sedges = np.linspace(0., rpmax, rpmax+1)

        if bin_type == 'lin':
            edges = (sedges, np.linspace(-40., 40, 101)) #transverse and radial separations are coded to be the same here
        else:
            edges = (sedges, np.linspace(0., 40., 41))
    elif corr_type == 'theta':
        edges = (np.linspace(0., 4., 101),)
    else:
        raise ValueError('corr_type must be one of ["smu", "rppi", "theta"]')
    return edges


def corr_fn(file_type='npy', region='', tracer='ELG', tracer2=None, zmin=0, zmax=np.inf, recon_dir='n',rec_type=False, weight_type='default', bin_type='lin', njack=0, nrandoms=8, split_randoms_above=10, out_dir='.', option=None, wang=None, rpcut=None, pimax=40):
    if tracer2: tracer += '_' + tracer2
    if rec_type: tracer += '_' + rec_type
    if region: tracer += '_' + region
    if option:
        zmax = str(zmax) + option
    #if recon_dir != 'n':
    #    out_dir += recon_dir+'/'
    split = '_split{:.0f}'.format(split_randoms_above) if split_randoms_above < np.inf else ''
    wang = '{}_'.format(wang) if wang is not None else ''
    root = '{}{}_{}_{}_{}_{}_njack{:d}_nran{:d}{}'.format(wang, tracer, zmin, zmax, weight_type, bin_type, njack, nrandoms, split)
    if rpcut is not None:
        root += '_rpcut{}'.format(rpcut)
    if file_type == 'npy':
        return os.path.join(out_dir, 'allcounts_{}.npy'.format(root))
    return os.path.join(out_dir, '{}_{}.txt'.format(file_type, root))


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--tracer', help='tracer(s) to be selected - 2 for cross-correlation', type=str, nargs='+', default=['LRG'])
    # Set basedir to the actual location of catalogs when starting with full (not clustering) catalogs for mocks. No need to set survey, verspec, version then.
    parser.add_argument('--basedir', help='where to find catalogs', type=str, default='/global/cfs/cdirs/desi/survey/catalogs/')
    parser.add_argument('--survey', help='e.g., SV3, DA02, etc.', type=str, default='Y1')
    parser.add_argument('--verspec', help='version for redshifts', type=str, default='iron')
    parser.add_argument('--version', help='catalog version', type=str, default='test')
    parser.add_argument("--use_map_veto",help="string to add on the end of full file reflecting if hp maps were used to cut",default='_HPmapcut')
    parser.add_argument('--region', help='regions; by default, run on N, S; pass NS to run on concatenated N + S', type=str, nargs='*', choices=['NGC','SGC'], default=None)
    parser.add_argument('--zlim', help='z-limits, or options for z-limits, e.g. "highz", "lowz", "fullonly"', type=str, nargs='*', default=None)
    parser.add_argument('--maglim', help='absolute r-band magnitude limits', type=str, nargs='*', default=None)
    parser.add_argument('--corr_type', help='correlation type', type=str, nargs='*', choices=[ 'rppi'], default=[ 'rppi'])
    parser.add_argument('--weight_type', help='types of weights to use; use "default_angular_bitwise" for PIP with angular upweighting; "default" just uses WEIGHT column', type=str, default='default')
    # Need to add support for fkp weights for use_arrays option
    parser.add_argument('--bin_type', help='binning type', type=str, choices=['log', 'lin'], default='lin')
    parser.add_argument('--nran', help='number of random files to combine together (1-18 available)', type=int, default=4)
    parser.add_argument('--split_ran_above', help='separation scale above which RR are summed over each random file;\
                                                   typically, most efficient for xi < 1, i.e. sep > 10 Mpc/h;\
                                                   see https://arxiv.org/pdf/1905.01133.pdf', type=float, default=80)#, default=20) # because currently the split doesn't work for multiple tracers
    parser.add_argument('--njack', help='number of jack-knife subsamples; 0 for no jack-knife error estimates', type=int, default=60)
    parser.add_argument('--nthreads', help='number of threads', type=int, default=64)
    parser.add_argument('--outdir', help='base directory for output (default: SCRATCH)', type=str, default=None)
    #parser.add_argument('--mpi', help='whether to use MPI', action='store_true', default=False)
    parser.add_argument('--vis', help='show plot of each xi?', action='store_true', default=False)
    parser.add_argument('--rebinning', help='whether to rebin the xi or just keep the original .npy file', default='y')
    # arguments relevant for when running directly from full catalogs.
    #parser.add_argument('--use_arrays', help = 'use pre-stored arrays rather than reading from memory again', default = 'y') 

    #only relevant for reconstruction
    parser.add_argument('--rec_type', help='reconstruction algorithm + reconstruction convention', choices=['IFTPrecsym', 'IFTPreciso','IFTrecsym', 'IFTreciso', 'MGrecsym', 'MGreciso'], type=str, default=None)
    parser.add_argument('--recon_dir', help='if recon catalogs are in a subdirectory, put that here', type=str, default='n')

    parser.add_argument('--rpcut', help='apply this rp-cut', type=float, default=None)
    
    parser.add_argument('--pimax', help='distance along the LOS to integrate, in Mpc/h', type=float, default=40)

    setup_logging()
    args = parser.parse_args()
    
    
    args.use_arrays = 'y'
    if args.rebinning == 'n':
        args.rebinning = False
    if args.rebinning == 'y':
        args.rebinning = True

    mpicomm, mpiroot = None, None
    if True:#args.mpi:
        from pycorr import mpi
        mpicomm = mpi.COMM_WORLD
        mpiroot = 0

    if os.path.normpath(args.basedir) == os.path.normpath('/global/cfs/cdirs/desi/survey/catalogs/'):
        cat_dir = io.catalog_dir(base_dir=args.basedir, survey=args.survey, verspec=args.verspec, version=args.version)
    elif os.path.normpath(args.basedir) == os.path.normpath('/global/project/projectdirs/desi/users/acarnero/mtl_mock000_univ1/'):
        cat_dir = args.basedir
        args.region = ['']
    else:
        cat_dir = args.basedir
    if mpicomm is None or mpicomm.rank == mpiroot:
        logger.info('Catalog directory is {}.'.format(cat_dir))

    if args.outdir is None:
        out_dir = os.path.join(io.get_scratch_dir(), args.survey)
    else:
        out_dir = args.outdir
    if mpicomm is None or mpicomm.rank == mpiroot:
        logger.info('Output directory is {}.'.format(out_dir))
        
    if args.pimax is None:
        pimax = 40
    else:
        pimax = args.pimax


    if args.use_arrays == 'y':
        print("Using arrays")
        from LSS.globals import main
        
        tracer, tracer2 = args.tracer[0], None
        if len(args.tracer) > 1:
            tracer2 = args.tracer[1]
            if len(args.tracer) > 2:
                raise ValueError('Provide <= 2 tracers!')
        if tracer2 == tracer:
            tracer2 = None # otherwise counting of self-pairs
            
        #tracer2 = None
        #tracer = args.tracer[0]
        data_ = []
        randoms_ = []
        for t in args.tracer:
            
            mainp = main(t,survey='Y1')
            dchi2 = mainp.dchi2
            #outaa = args.outdir
            #flaa = args.basedir
            #outaa = outaa + "/" + tracer
            flaa = cat_dir + "/" + t
            flinr = cat_dir + "/" + t + "_"
            if t == 'BGS_BRIGHT-21.5':
                flinr = cat_dir + "/BGS_BRIGHT_"

            #rann = 1

            if t == "LRG":
                zminr = 0.4
                zmaxr = 1.1
            if t[:3] == "ELG":
                zminr = 0.8
                zmaxr = 1.6
            if t[:3] == "BGS":
                zminr = 0.1
                zmaxr = 0.4

            if t == "QSO":
                zminr = 0.8
                zmaxr = 3.5
            rcols = ['Z', 'WEIGHT', 'WEIGHT_SYS', 'WEIGHT_COMP', 'WEIGHT_ZFAIL']#,'WEIGHT_FKP']
            if 'FKP' in args.weight_type:
                rcols.append('WEIGHT_FKP')
            datai_ = ct.mkclusdat(flaa,weighttileloc=True,zmask=False,tp=t,dchi2=dchi2,tsnrcut=0,rcut=None,ntilecut=0,ccut=None,ebits=None,zmin=zminr,zmax=zmaxr,write_cat='n',return_cat='y',use_map_veto=args.use_map_veto)
            print('data columns',datai_.dtype.names)
            
            ranl =[]
            for rann in range(0,args.nran):
                rani = ct.mkclusran(flinr,flinr,rann,rcols=rcols,zmask=False,tsnrcut=0,tsnrcol='TSNR2_ELG',utlid=False,ebits=None,write_cat='n',return_cat='y', clus_arrays = datai_,use_map_veto=args.use_map_veto)
                ranl.append(np.array(rani))
            randomsi_ = np.concatenate(ranl)
            #out_dir = args.outdir
            
            data_.append(datai_)
            randoms_.append(randomsi_)

    
    elif args.use_arrays == 'n':
        print("use_arrays set to false")

        tracer, tracer2 = args.tracer[0], None
        if len(args.tracer) > 1:
            tracer2 = args.tracer[1]
            if len(args.tracer) > 2:
                raise ValueError('Provide <= 2 tracers!')
        if tracer2 == tracer:
            tracer2 = None # otherwise counting of self-pairs
        catalog_kwargs = dict(tracer=tracer, tracer2=tracer2, survey=args.survey, cat_dir=cat_dir, recon_dir=args.recon_dir,rec_type=args.rec_type) # survey required for zdone
        
    distance = TabulatedDESI().comoving_radial_distance

    regions = args.region
    if regions is None:
        regions = io.get_regions(args.survey, rec=bool(args.rec_type))

    option = None
    if args.zlim is None:
        zlims = io.get_zlims(tracer, tracer2=tracer2)
    elif not args.zlim[0].replace('.', '').isdigit():
        option = args.zlim[0]
        zlims = io.get_zlims(tracer, tracer2=tracer2, option=option)
    else:
        zlims = [float(zlim) for zlim in args.zlim]


    if args.maglim is not None:
        magmin = float(args.maglim[0])
        magmax = float(args.maglim[1])
        maglims = (magmin,magmax)
    else:
        maglims = None

    zlims = list(zip(zlims[:-1], zlims[1:])) + ([(zlims[0], zlims[-1])] if len(zlims) > 2 else []) # len(zlims) == 2 == single redshift range
    rebinning_factors = [1, 4, 5, 10] if 'lin' in args.bin_type else [1, 2, 4]
    pi_rebinning_factors = [1, 4, 5, 10] if 'log' in args.bin_type else [1]
    if mpicomm is None or mpicomm.rank == mpiroot:
        logger.info('Computing correlation functions {} in regions {} in redshift ranges {}.'.format(args.corr_type, regions, zlims))

    for zmin, zmax in zlims:
        base_file_kwargs = dict(tracer=tracer, tracer2=tracer2, zmin=zmin, zmax=zmax, recon_dir=args.recon_dir,rec_type=args.rec_type, weight_type=args.weight_type, bin_type=args.bin_type, njack=args.njack, nrandoms=args.nran, split_randoms_above=args.split_ran_above, option=option, rpcut=args.rpcut)
        for region in regions:
            if args.use_arrays == 'y':
                sel_ngc_dat = common.splitGC(data_[0])
                sel_ngc_ran = common.splitGC(randoms_[0])
                fcdn2=None; randoms_gc2=None
                
                if len(args.tracer) > 1:
                    sel_ngc_dat2 = common.splitGC(data_[1])
                    sel_ngc_ran2 = common.splitGC(randoms_[1])
                
                if region == 'NGC':
                    ffr = randoms_[0][sel_ngc_ran]
                    fcdn = data_[0][sel_ngc_dat]
                    randoms_gc = ct.clusran_resamp_arrays(ffr,fcdn,region,tracer,rcols=rcols)
                    
                    if len(args.tracer) > 1:
                        ffr2 = randoms_[1][sel_ngc_ran2]
                        fcdn2 = data_[1][sel_ngc_dat2]
                        randoms_gc2 = ct.clusran_resamp_arrays(ffr2,fcdn2,region,tracer2,rcols=rcols)
                    
                    catalog_kwargs = dict(tracer=tracer, tracer2=tracer2, recon_dir=args.recon_dir, rec_type=args.rec_type, cat_read='Y', 
                                          dat_cat=fcdn, ran_cat=randoms_gc, dat_cat2=fcdn2, ran_cat2=randoms_gc2)
                if region == 'SGC':
                    ffr = randoms_[0][~sel_ngc_ran]
                    fcdn = data_[0][~sel_ngc_dat]
                    randoms_gc = ct.clusran_resamp_arrays(ffr,fcdn,region,tracer,rcols=rcols)
                    
                    if len(args.tracer) > 1:
                        ffr2 = randoms_[1][~sel_ngc_ran2]
                        fcdn2 = data_[1][~sel_ngc_dat2]
                        randoms_gc2 = ct.clusran_resamp_arrays(ffr2,fcdn2,region,tracer2,rcols=rcols)
                        
                    catalog_kwargs = dict(tracer=tracer, tracer2=tracer2, recon_dir=args.recon_dir, rec_type=args.rec_type, cat_read='Y', 
                                          dat_cat=fcdn, ran_cat=randoms_gc, dat_cat2=fcdn2, ran_cat2=randoms_gc2)

                #if region == "N":
                #    catalog_kwargs = dict(tracer=tracer, tracer2=tracer2, recon_dir=args.recon_dir, rec_type=args.rec_type, cat_read='Y', dat_cat=data_[0], ran_cat=randoms_[0])
                #if region == "S":
                #    catalog_kwargs = dict(tracer=tracer, tracer2=tracer2, recon_dir=args.recon_dir, rec_type=args.rec_type, cat_read='Y', dat_cat=data_[1], ran_cat=randoms_[1])
                
            wang = None
            for corr_type in args.corr_type:
                if mpicomm is None or mpicomm.rank == mpiroot:
                    logger.info('Computing correlation function {} in region {} in redshift range {}.'.format(corr_type, region, (zmin, zmax)))
                edges = get_edges(corr_type=corr_type, bin_type=args.bin_type)
            
                result, wang = compute_correlation_function(corr_type, edges=edges, distance=distance, nrandoms=args.nran, split_randoms_above=args.split_ran_above, nthreads=args.nthreads, region=region, zlim=(zmin, zmax), maglim=maglims, weight_type=args.weight_type, njack=args.njack, wang=wang, mpicomm=mpicomm, mpiroot=mpiroot, option=option, rpcut=args.rpcut, **catalog_kwargs)
                print('Finished compute_correlation_function')
                # Save pair counts
                if mpicomm is None or mpicomm.rank == mpiroot:
                    result.save(corr_fn(file_type='npy', region=region, out_dir=os.path.join(out_dir, corr_type), **base_file_kwargs))
            if mpicomm is None or mpicomm.rank == mpiroot:
                 if wang is not None:
                        for name in wang:
                            if wang[name] is not None:
                                wang[name].save(corr_fn(file_type='npy', region=region, out_dir=os.path.join(out_dir, 'wang'), **base_file_kwargs, wang=name))

        # Save combination and .txt files
        for corr_type in args.corr_type:
            all_regions = regions.copy()
            if mpicomm is None or mpicomm.rank == mpiroot:
                if 'N' in regions and 'S' in regions:  # let's combine
                    result = sum([TwoPointCorrelationFunction.load(
                                  corr_fn(file_type='npy', region=region, out_dir=os.path.join(out_dir, corr_type), **base_file_kwargs)).normalize() for region in ['N', 'S']])
                    result.save(corr_fn(file_type='npy', region='NScomb', out_dir=os.path.join(out_dir, corr_type), **base_file_kwargs))
                    all_regions.append('NScomb')
                if 'NGC' in regions and 'SGC' in regions:  # let's combine
                    result = sum([TwoPointCorrelationFunction.load(
                                  corr_fn(file_type='npy', region=region, out_dir=os.path.join(out_dir, corr_type), **base_file_kwargs)).normalize() for region in ['NGC', 'SGC']])
                    result.save(corr_fn(file_type='npy', region='GCcomb', out_dir=os.path.join(out_dir, corr_type), **base_file_kwargs))
                    all_regions.append('GCcomb')

                if args.rebinning:
                    for region in all_regions:
                        txt_kwargs = base_file_kwargs.copy()
                        txt_kwargs.update(region=region, out_dir=os.path.join(out_dir, corr_type))
                        result = TwoPointCorrelationFunction.load(corr_fn(file_type='npy', **txt_kwargs))
                        for factor in rebinning_factors:
                            #result = TwoPointEstimator.load(fn)
                            rebinned = result[:(result.shape[0] // factor) * factor:factor]
                            txt_kwargs.update(bin_type=args.bin_type+str(factor))
                            if corr_type == 'smu':
                                fn_txt = corr_fn(file_type='xismu', **txt_kwargs)
                                rebinned.save_txt(fn_txt)
                                fn_txt = corr_fn(file_type='xipoles', **txt_kwargs)
                                rebinned.save_txt(fn_txt, ells=(0, 2, 4), ignore_nan=True)
                                fn_txt = corr_fn(file_type='xiwedges', **txt_kwargs)
                                rebinned.save_txt(fn_txt, wedges=(-1., -2./3, -1./3, 0., 1./3, 2./3, 1.))
                            elif corr_type == 'rppi':
                                fn_txt = corr_fn(file_type='wp', **txt_kwargs)
                                rebinned.save_txt(fn_txt, pimax=pimax)
                                for pifac in pi_rebinning_factors:
                                    rebinned = result[:(result.shape[0]//factor)*factor:factor,:(result.shape[1]//pifac)*pifac:pifac]
                                    txt_kwargs.update(bin_type=args.bin_type+str(factor)+'_'+str(pifac))
                                    fn_txt = corr_fn(file_type='xirppi', **txt_kwargs)
                                    rebinned.save_txt(fn_txt)
                            elif corr_type == 'theta':
                                fn_txt = corr_fn(file_type='theta', **txt_kwargs)
                                rebinned.save_txt(fn_txt)

                            if args.vis:
                                if corr_type == 'smu':
                                    sep, xis = rebinned(ells=(0, 2, 4), return_sep=True, return_std=False)
                                elif corr_type == 'rppi':
                                    sep, xis = rebinned(pimax=pimax, return_sep=True, return_std=False)
                                else:
                                    sep, xis = rebinned(return_sep=True, return_std=False)
                                if args.bin_type == 'log':
                                    for xi in xis: plt.loglog(sep, xi)
                                if args.bin_type == 'lin':
                                    for xi in xis: plt.plot(sep, sep**2 * xi)
                                tracers = tracer
                                if tracer2 is not None: tracers += ' x ' + tracer2
                                plt.title('{} {:.2f} < z {:.2f} in {}'.format(tracers, zmin, zmax, region))
                                plt.show()

            os.system('rm '+os.path.join(out_dir, corr_type)+'/xirp*')
            os.system('rm '+os.path.join(out_dir, corr_type)+'/allcounts*')