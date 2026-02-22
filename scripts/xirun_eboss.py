#!/usr/bin/env python
# coding: utf-8

"""
Note
----
The script can be called with multiple processes as (e.g. on 2 nodes, 64 threads for each):
```
srun -n 2 python xirunpc.py --nthreads 64 ...
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

from cosmoprimo.fiducial import BOSS
distance = BOSS(engine='class').comoving_radial_distance


logger = logging.getLogger('xirunpc')


import os
import argparse
import logging
import numpy as np
from astropy.table import Table, vstack
from matplotlib import pyplot as plt


# from pycorr import TwoPointCorrelationFunction, TwoPointEstimator, utils, project_to_multipoles, project_to_wp, setup_logging
# from Cosmo import distance
# 
# ds = distance(.31,.69) #eboss fiducial cosmo
# 
# parser = argparse.ArgumentParser()
# parser.add_argument("--type", help="tracer type to be selected",default='QSO')
# parser.add_argument("--basedir", help="where to find catalogs",default='/global/cscratch1/sd/ajross/ebosscat/')
# parser.add_argument("--outdir", help="where to output results",default='/global/cscratch1/sd/ajross/ebossxi/')
# parser.add_argument("--bintype",help="log or lin",default='lin')
# parser.add_argument("--nthreads",help="number of threads for parallel comp",default=32,type=int)
# parser.add_argument("--vis",help="set to y to plot each xi ",default='n')
# 
# #only relevant for reconstruction
# parser.add_argument("--rectype",help="IFT or MG supported so far",default='IFT')
# parser.add_argument("--convention",help="recsym or reciso supported so far",default='reciso')
# 
# setup_logging()
# args = parser.parse_args()
# 
# ttype = args.type
# lssdir = args.basedir
# 
# if args.bintype == 'log':
#     bine = np.logspace(-1.5, 2.2, 80)
# if args.bintype == 'lin':
#     bine = np.linspace(1e-4, 200, 201)
# 
# dirxi = args.outdir
# 
# 
# dirname = lssdir 
# 
# def compute_correlation_function(mode, tracer='QSO', region='NGC', zlim=(0., np.inf), nthreads=8, dtype='f8', wang=None,fnroot=''):
#     data_fn = os.path.join(dirname, 'eBOSS_{}_clustering_data-{}-vDR16.fits'.format(tracer, region))
#     data = Table.read(data_fn)
# 
#     shifted = None
# 
#     randoms_fn = os.path.join(dirname, 'eBOSS_{}_clustering_random-{}-vDR16.fits'.format(tracer, region)) 
#     randoms = Table.read(randoms_fn) 
#   
#     corrmode = mode
#     if mode == 'wp':
#         corrmode = 'rppi'
#     if mode == 'multi':
#         corrmode = 'smu'
#     if corrmode == 'smu':
#         edges = (bine, np.linspace(-1., 1., 201)) #s is input edges and mu evenly spaced between 0 and 1
#     if corrmode == 'rppi':
#         edges = (bine, bine) #transverse and radial separations are  coded to be the same here
#         if mode == 'wp':
#             edges = (bins,np.linspace(0,40.,41)) #if you want wp, only go out to pi = 40; consider setting pi_max as argument
#     
#     data_positions, data_weights = get_positions_weights(data)
#     randoms_positions, randoms_weights = get_positions_weights(randoms)
#     shifted_positions, shifted_weights = None, None
#     if shifted is not None:
#         shifted_positions, shifted_weights = get_positions_weights(shifted, name='shifted')
# 
#     kwargs = {}
# 
#     result = TwoPointCorrelationFunction(corrmode, edges, data_positions1=data_positions, data_weights1=data_weights,
#                                          randoms_positions1=randoms_positions, randoms_weights1=randoms_weights,
#                                          shifted_positions1=shifted_positions, shifted_weights1=shifted_weights,
#                                          engine='corrfunc', position_type='rdd', nthreads=nthreads, dtype=dtype, **kwargs)
#     #save paircounts
#     fn = dirxi+'paircounts_'+fnroot+'.npy'
#     result.save(fn)
#     return result
# 
# ranwt1=False
# 
# regl = ['NGC','SGC']
# 
# tcorr = ttype
# tw = ttype
# 
# bsl = [1,4,5,10]
# 
# if ttype == 'QSO':
#     zmin = 0.8
#     zmax = 2.2
# 
# for reg in regl:
#   print(reg)
#   #(sep, xiell), wang = compute_correlation_function('multi', bs, tracer=tcorr, region=reg, nrandoms=args.nran, zlim=(zmin,zmax), weight_type=weight_type,nthreads=args.nthreads)
#   fnroot = tw+'eboss'+reg+'_'+str(zmin)+str(zmax)+'DR16'+'_'+args.bintype
#   pfn = dirxi+'paircounts_'+fnroot+'.npy'
#   result = compute_correlation_function('multi', tracer=tcorr, region=reg, zlim=(zmin,zmax),nthreads=args.nthreads,fnroot=fnroot)
#   for bs in bsl:
#       result = TwoPointEstimator.load(pfn)
#       result.rebin((bs, 1))
#       sep,xiell = project_to_multipoles(result)#, wang
#       fo = open(dirxi+'xi024'+fnroot+str(bs)+'.dat','w')
#       for i in range(0,len(sep)):
#           fo.write(str(sep[i])+' '+str(xiell[0][i])+' '+str(xiell[1][i])+' '+str(xiell[2][i])+'\n')
#       fo.close()
#       if args.vis == 'y':
#           if args.bintype == 'log':
#               plt.loglog(sep,xiell[0])
#           if args.bintype == 'lin':
#               plt.plot(sep,sep**2.*xiell[0])
#           plt.title(ttype+' '+str(zmin)+'<z<'+str(zmax)+' in '+reg)
#           plt.show()    


dirname = '/global/cfs/cdirs/sdss/data/sdss/dr16/eboss/lss/catalogs/DR16/'


def compute_correlation_function(corr_type, edges, zlim,region,nthreads=8, gpu=False, dtype='f8', wang=None, split_randoms_above=30., weight_type='default', tracer='ELG', tracer2=None, recon_dir=None, rec='_rec', njack=120, option=None, mpicomm=None, mpiroot=None, cat_read=None, dat_cat=None, ran_cat=None, rpcut=None):
    data_fn = os.path.join(dirname, 'eBOSS_{}_clustering_data{}-{}-vDR16.fits'.format(tracer,rec, region))
    data = Table.read(data_fn)

    randoms_fn = os.path.join(dirname, 'eBOSS_{}_clustering_random-{}-vDR16.fits'.format(tracer, region)) 
    randoms = Table.read(randoms_fn) 
    
    if rec == '_rec':
        randoms_fn = os.path.join(dirname, 'eBOSS_{}_clustering_random{}-{}-vDR16.fits'.format(tracer, rec,region)) 
        randoms_rec = Table.read(randoms_fn) 

    def get_positions_weights(catalog):
        mask = (catalog['Z'] >= zlim[0]) & (catalog['Z'] < zlim[1])
        print('Using {:d} rows '.format(mask.sum()))
        positions = [catalog['RA'][mask], catalog['DEC'][mask], distance(catalog['Z'][mask])]
        weights = catalog['WEIGHT_FKP'][mask]*catalog['WEIGHT_SYSTOT'][mask]*catalog['WEIGHT_CP'][mask]*catalog['WEIGHT_NOZ'][mask] 
                
        return positions, weights

    autocorr = tracer2 is None


    data_positions1, data_weights1, data_samples1, data_positions2, data_weights2, data_samples2 = None, None, None, None, None, None
    randoms_positions1, randoms_weights1, randoms_samples1, randoms_positions2, randoms_weights2, randoms_samples2 = None, None, None, None, None, None
    shifted_positions1, shifted_weights1, shifted_samples1, shifted_positions2, shifted_weights2, shifted_samples2 = None, None, None, None, None, None
    jack_positions = None

    if mpicomm is None or mpicomm.rank == mpiroot:

        data_positions1, data_weights1 = get_positions_weights(data)
        randoms_positions1, randoms_weights1 = get_positions_weights(randoms)
        if rec == '_rec':
            shifted_positions1, shifted_weights1 = get_positions_weights(randoms_rec)

        jack_positions = data_positions1



    kwargs = {}
    selection_attrs = None
    if rpcut is not None: selection_attrs = {'rp': (rpcut, np.inf)}
    randoms_kwargs = dict(randoms_positions1=randoms_positions1, randoms_weights1=randoms_weights1, randoms_samples1=randoms_samples1,
                          randoms_positions2=randoms_positions2, randoms_weights2=randoms_weights2, randoms_samples2=randoms_samples2,
                          shifted_positions1=shifted_positions1, shifted_weights1=shifted_weights1, shifted_samples1=shifted_samples1,
                          shifted_positions2=shifted_positions2, shifted_weights2=shifted_weights2, shifted_samples2=shifted_samples2)


   
    if mpicomm is None:
        nran = len(randoms_positions1)
    else:
        nran = mpicomm.bcast(len(randoms_positions1) if mpicomm.rank == mpiroot else None, root=mpiroot)
    result = TwoPointCorrelationFunction(corr_type, edges, data_positions1=data_positions1,data_weights1=data_weights1,\
    randoms_positions1=randoms_positions1, randoms_weights1=randoms_weights1,shifted_positions1=shifted_positions1, shifted_weights1=shifted_weights1,\
    engine='corrfunc', position_type='rdd', nthreads=nthreads, gpu=gpu, dtype=dtype, **kwargs,\
    mpicomm=mpicomm, mpiroot=mpiroot, selection_attrs=selection_attrs)
    return result


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
        if bin_type == 'lin':
            edges = (sedges, np.linspace(-40., 40, 101)) #transverse and radial separations are coded to be the same here
        else:
            edges = (sedges, np.linspace(0., 40., 41))
    elif corr_type == 'theta':
        edges = (np.linspace(0., 4., 101),)
    else:
        raise ValueError('corr_type must be one of ["smu", "rppi", "theta"]')
    return edges


# def corr_fn(file_type='npy', region='', tracer='ELG', tracer2=None, zmin=0, zmax=np.inf, recon_dir='n',rec_type=False, weight_type='default', bin_type='lin', njack=0, nrandoms=8, split_randoms_above=10, out_dir='.', option=None, wang=None, rpcut=None):
#     if tracer2: tracer += '_' + tracer2
#     if rec_type: tracer += '_' + rec_type
#     if region: tracer += '_' + region
#     if option:
#         zmax = str(zmax) + option
#     #if recon_dir != 'n':
#     #    out_dir += recon_dir+'/'
#     split = '_split{:.0f}'.format(split_randoms_above) if split_randoms_above < np.inf else ''
#     wang = '{}_'.format(wang) if wang is not None else ''
#     root = '{}{}_{}_{}_{}_{}_njack{:d}_nran{:d}{}'.format(wang, tracer, zmin, zmax, weight_type, bin_type, njack, nrandoms, split)
#     if rpcut is not None:
#         root += '_rpcut{}'.format(rpcut)
#     if file_type == 'npy':
#         return os.path.join(out_dir, 'allcounts_{}.npy'.format(root))
#     return os.path.join(out_dir, '{}_{}.txt'.format(file_type, root))


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--tracer', help='tracer(s) to be selected - 2 for cross-correlation', type=str, nargs='+', default=['LRG'])
    # Set basedir to the actual location of catalogs when starting with full (not clustering) catalogs for mocks. No need to set survey, verspec, version then.
    #parser.add_argument('--region', help='regions; by default, run on N, S; pass NS to run on concatenated N + S', type=str, nargs='*', choices=['N', 'S', 'NS','NGC','SGC','NGCS','DES','SGCnotDES'], default=None)
    parser.add_argument('--zlim', help='z-limits, or options for z-limits, e.g. "highz", "lowz", "fullonly"', type=str, nargs='*', default=None)
    parser.add_argument('--maglim', help='absolute r-band magnitude limits', type=str, nargs='*', default=None)
    parser.add_argument('--option', help='place to put extra options for cutting catalogs', default=None)
    parser.add_argument('--corr_type', help='correlation type', type=str, nargs='*', choices=['smu', 'rppi', 'theta'], default='smu')
    parser.add_argument('--weight_type', help='types of weights to use; use "default_angular_bitwise" for PIP with angular upweighting; "default" just uses WEIGHT column', type=str, default='default')
    # Need to add support for fkp weights for use_arrays option
    parser.add_argument('--bin_type', help='binning type', type=str, choices=['log', 'lin'], default='lin')
    parser.add_argument('--njack', help='number of jack-knife subsamples; 0 for no jack-knife error estimates', type=int, default=0)
    parser.add_argument('--gpu', help='whether to run on the GPU', action='store_true')
    parser.add_argument('--nthreads', help='number of threads (defaults to 4 if --gpu else 128)', type=int, default=None)
    parser.add_argument('--outdir', help='base directory for output (default: SCRATCH)', type=str, default=None)
    parser.add_argument('--vis', help='show plot of each xi?', action='store_true', default=False)
    parser.add_argument('--rebinning', help='whether to rebin the xi or just keep the original .npy file', default='y')
    parser.add_argument('--rec', type=str, default='_rec')

    setup_logging()
    args = parser.parse_args()

    gpu, nthreads = args.gpu, args.nthreads
    if nthreads is None:
        if gpu: nthreads = 4
        else: nthreads = 128
    
    if args.rebinning == 'n':
        args.rebinning = False
    if args.rebinning == 'y':
        args.rebinning = True

    mpicomm, mpiroot = None, None
    if True:#args.mpi:
        from pycorr import mpi
        mpicomm = mpi.COMM_WORLD
        mpiroot = 0

    print("use_arrays set to false")

    if args.outdir is None:
        out_dir = os.path.join(os.getenv('SCRATCH'), 'ebossxi')
    else:
        out_dir = args.outdir
    if mpicomm is None or mpicomm.rank == mpiroot:
        logger.info('Output directory is {}.'.format(out_dir))

    tracer, tracer2 = args.tracer[0], None
    if len(args.tracer) > 1:
        tracer2 = args.tracer[1]
        if len(args.tracer) > 2:
            raise ValueError('Provide <= 2 tracers!')
    if tracer2 == tracer:
        tracer2 = None # otherwise counting of self-pairs
        
    

    regions = ['NGC','SGC']
    if 'LRG' in tracer:
        zmin = 0.6
        zmax = 1.0
    rebinning_factors = [1, 4, 5, 10] if 'lin' in args.bin_type else [1, 2, 4]
    pi_rebinning_factors = [1, 4, 5, 10] if 'log' in args.bin_type else [1]
    if mpicomm is None or mpicomm.rank == mpiroot:
        logger.info('Computing correlation functions {} in regions {} in redshift ranges {}.'.format(args.corr_type, regions, (zmin,zmax)))
    corr_type = args.corr_type
    for region in regions:
        if mpicomm is None or mpicomm.rank == mpiroot:
            root = '{}_{}_{}{}_{}'.format(tracer, zmin, zmax, args.rec,region)
            fout = os.path.join(out_dir, 'allcounts_{}.npy'.format(root))

        if mpicomm is None or mpicomm.rank == mpiroot:
                logger.info('Computing correlation function {} in region {} in redshift range {}.'.format(corr_type, region, (zmin, zmax)))
        edges = get_edges(corr_type=corr_type, bin_type=args.bin_type)
            
        result = compute_correlation_function(corr_type, edges=edges, zlim=(zmin,zmax), region=region,nthreads=nthreads, gpu=gpu, njack=args.njack, mpicomm=mpicomm, mpiroot=mpiroot,tracer=tracer)
                # Save pair counts
        if mpicomm is None or mpicomm.rank == mpiroot:
            
            result.save(fout)

    all_regions = regions.copy()
    if mpicomm is None or mpicomm.rank == mpiroot:
        if 'NGC' in regions and 'SGC' in regions:  # let's combine
            result = sum([TwoPointCorrelationFunction.load(os.path.join(out_dir, 'allcounts_{}_{}_{}{}_{}.npy'.format(tracer, zmin, zmax, args.rec,region))).normalize() for region in ['NGC', 'SGC']])
            result.save(os.path.join(out_dir, 'allcounts_{}_{}_{}{}_{}.npy'.format(tracer, zmin, zmax, args.rec,'GCcomb')))
            all_regions.append('GCcomb')

#           if args.rebinning:
#               for region in all_regions:
#                   txt_kwargs = base_file_kwargs.copy()
#                   txt_kwargs.update(region=region, out_dir=os.path.join(out_dir, corr_type))
#                   result = TwoPointCorrelationFunction.load(corr_fn(file_type='npy', **txt_kwargs))
#                   for factor in rebinning_factors:
#                       #result = TwoPointEstimator.load(fn)
#                       rebinned = result[:(result.shape[0] // factor) * factor:factor]
#                       txt_kwargs.update(bin_type=args.bin_type+str(factor))
#                       if corr_type == 'smu':
#                           fn_txt = corr_fn(file_type='xismu', **txt_kwargs)
#                           rebinned.save_txt(fn_txt)
#                           fn_txt = corr_fn(file_type='xipoles', **txt_kwargs)
#                           rebinned.save_txt(fn_txt, ells=(0, 2, 4), ignore_nan=True)
#                           fn_txt = corr_fn(file_type='xiwedges', **txt_kwargs)
#                           rebinned.save_txt(fn_txt, wedges=(-1., -2./3, -1./3, 0., 1./3, 2./3, 1.))
#                       elif corr_type == 'rppi':
#                           fn_txt = corr_fn(file_type='wp', **txt_kwargs)
#                           rebinned.save_txt(fn_txt, pimax=40.)
#                           for pifac in pi_rebinning_factors:
#                               rebinned = result[:(result.shape[0]//factor)*factor:factor,:(result.shape[1]//pifac)*pifac:pifac]
#                               txt_kwargs.update(bin_type=args.bin_type+str(factor)+'_'+str(pifac))
#                               fn_txt = corr_fn(file_type='xirppi', **txt_kwargs)
#                               rebinned.save_txt(fn_txt)
#                       elif corr_type == 'theta':
#                           fn_txt = corr_fn(file_type='theta', **txt_kwargs)
#                           rebinned.save_txt(fn_txt)
# 
#                       if args.vis:
#                           if corr_type == 'smu':
#                               sep, xis = rebinned(ells=(0, 2, 4), return_sep=True, return_std=False)
#                           elif corr_type == 'rppi':
#                               sep, xis = rebinned(pimax=40, return_sep=True, return_std=False)
#                           else:
#                               sep, xis = rebinned(return_sep=True, return_std=False)
#                           if args.bin_type == 'log':
#                               for xi in xis: plt.loglog(sep, xi)
#                           if args.bin_type == 'lin':
#                               for xi in xis: plt.plot(sep, sep**2 * xi)
#                           tracers = tracer
#                           if tracer2 is not None: tracers += ' x ' + tracer2
#                           plt.title('{} {:.2f} < z {:.2f} in {}'.format(tracers, zmin, zmax, region))
#                           plt.show()
