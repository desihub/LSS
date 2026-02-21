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
import fitsio
from astropy.table import Table, vstack, join
from matplotlib import pyplot as plt

from pycorr import TwoPointCorrelationFunction, TwoPointEstimator, KMeansSubsampler, utils, setup_logging

from LSS.tabulated_cosmo import TabulatedDESI
import LSS.common_tools as common
#import LSS.main.cattools as ct



logger = logging.getLogger('cosmodesi_io')


def get_scratch_dir():
    if os.environ['NERSC_HOST'] == 'cori':
        scratch_dir = os.environ['CSCRATCH']
        os.system('export OMP_NUM_THREADS=64')
    elif os.environ['NERSC_HOST'] == 'perlmutter':
        scratch_dir = os.environ['PSCRATCH']
        os.system('export OMP_NUM_THREADS=128')
    else:
        msg = 'NERSC_HOST is not cori or permutter but is {};\n'.format(os.environ['NERSC_HOST'])
        msg += 'NERSC_HOST not known (code only works on NERSC), not proceeding'
        raise ValueError(msg)
    return scratch_dir


def get_zlims(tracer, tracer2=None, option=None):

    if tracer2 is not None:
        zlims1 = get_zlims(tracer, option=option)
        zlims2 = get_zlims(tracer2, option=option)
        return [zlim for zlim in zlims1 if zlim in zlims2]

    if tracer.startswith('LRG'):
        zlims = [0.4, 0.6, 0.8, 1.1]

    if tracer.startswith('ELG'):# or type == 'ELG_HIP':
        zlims = [0.8, 1.1, 1.6] #[1.5,1.6]
        if option:
            if option == 'safez':
                zlims = [0.9, 1.48]
            if 'extended' in option:
                logger.warning('extended is no longer a meaningful option')
                #zlims = [0.8, 1.1, 1.6]
            if 'smallshells' in option:
                zlims = [0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6] 

    if tracer.startswith('QSO'):
        #zlims = [0.8, 1.1, 1.6, 2.1]
        zlims = [0.8,2.1]
        if option == 'highz':
            zlims = [2.1, 3.5]
        if option == 'lowz':
            zlims = [0.8, 2.1]

    if tracer.startswith('BGS'):
        zlims = [0.1, 0.4]
        #if option == 'lowz':
        #    zlims = [0.1, 0.3]
        #if option == 'highz':
        #    zlims = [0.3, 0.5]

    if option == 'fullonly':
        zlims = [zlims[0], zlims[-1]]

    return zlims


def get_regions(survey, rec=False):
    regions = ['NGC', 'SGC']#, '']
    #if survey in ['main', 'DA02']:
    #    regions = ['DN', 'DS', 'N', 'S']
    #    if rec: regions = ['DN', 'N']
    return regions


def select_region(ra, dec, region):
    mask_ra = (ra > 100 - dec)
    mask_ra &= (ra < 280 + dec)
    #mask = None
    if region == 'DN':
        mask = dec < 32.375
        mask &= mask_ra
    elif region == 'DS':
        mask = dec > -25
        mask &= ~mask_ra
    elif 'GC' in region:
        from astropy.coordinates import SkyCoord
        import astropy.units as u
        c = SkyCoord(ra* u.deg,dec* u.deg,frame='icrs')
        gc = c.transform_to('galactic')
        sel_ngc = gc.b > 0
        if region == 'NGC':
            mask = sel_ngc
        if region == 'SGC':
            mask = ~sel_ngc        
    else:
        raise ValueError('Input region must be one of ["DN", "DS","NGC","SGC"].')
    return mask


def catalog_dir(survey='Y1', verspec='iron', version='v1.2', base_dir='/global/cfs/cdirs/desi/survey/catalogs'):
    return os.path.join(base_dir, survey, 'LSS', verspec, 'LSScats', version)



def catalog_fn(tracer='ELG', region='', ctype='clustering', name='data', ran_sw='',recon_dir='n',rec_type=False, nrandoms=4, cat_dir=None, survey='Y1', **kwargs):
    #print(kwargs)
    if cat_dir is None:
        cat_dir = catalog_dir(survey=survey, **kwargs).replace('global','dvs_ro')
    #if survey in ['main', 'DA02']:
    #    tracer += 'zdone'
    if 'edav1' in cat_dir:
        cat_dir += ctype
    #if ctype == 'clustering':
    #    cat_dir += '/unblinded/'           
    if 'full' in ctype:# == 'full':
        # For subsamples, the tracer name is augmented with details about the subsample
        # That will be needed in most places, but for angular upweighting we use the parent sample, not the subsample to find the angular weights.
        # So strip it off just here
        # You may need to adjust this
        region = ''
        cat_dir = cat_dir.replace('/unblinded','')
        cat_dir = cat_dir.replace('/blinded','')
        cat_dir = cat_dir.replace('PIP','')
        logger.info('cat_dir for full file is '+cat_dir)
        if 'BGS_BRIGHT' in tracer:
            tracer = 'BGS_BRIGHT'
            logger.info('reset tracer name to BGS_BRIGHT for reading full file')
        if 'LRG' in tracer:
            tracer = 'LRG'
            logger.info('reset tracer name to LRG for reading full file')
        
    dat_or_ran = name[:3]
    if name == 'randoms' and tracer == 'LRG_main' and ctype == 'full':
        tracer = 'LRG'
    if region: region = '_' + region
    #recon_dir = kwargs['recon_dir']
    if recon_dir != 'n':
        tracer = recon_dir+'/'+tracer
    if rec_type:
        dat_or_ran = '{}.{}'.format(rec_type, dat_or_ran)
    if name == 'data':
        fname = os.path.join(cat_dir, '{}{}_{}.{}.fits'.format(tracer, region, ctype, dat_or_ran))
        fname5 = os.path.join(cat_dir, '{}{}_{}.{}.h5'.format(tracer, region, ctype, dat_or_ran))
        if os.path.isfile(fname):
            return fname
        elif os.path.isfile(fname5):
            return fname5
        else:
            return fname +' (or .h5) file_not_found!!!'
    #print(nrandoms)
    test_fname = os.path.join(cat_dir, '{}{}{}_0_{}.{}.fits'.format(tracer, ran_sw, region, ctype, dat_or_ran))
    test_fname5 = os.path.join(cat_dir, '{}{}{}_0_{}.{}.h5'.format(tracer, ran_sw, region, ctype, dat_or_ran))
    if os.path.isfile(test_fname):
        return [os.path.join(cat_dir, '{}{}{}_{:d}_{}.{}.fits'.format(tracer, ran_sw, region, iran, ctype, dat_or_ran)) for iran in range(nrandoms)]
    elif os.path.isfile(test_fname5):
        return [os.path.join(cat_dir, '{}{}{}_{:d}_{}.{}.h5'.format(tracer, ran_sw, region, iran, ctype, dat_or_ran)) for iran in range(nrandoms)]


def _format_bitweights(bitweights):
    if bitweights.ndim == 2: return list(bitweights.T)
    return [bitweights]


def get_clustering_positions_weights(catalog, distance, zlim=(0., np.inf),fac_ntmp=None,maglim=None, weight_type='default', name='data', return_mask=False, option=None,P0=None):
    #logger.info('get pos P0 is '+str(P0))
    if maglim is None:
        mask = (catalog['Z'] >= zlim[0]) & (catalog['Z'] < zlim[1])
    if maglim is not None:
        mask = (catalog['Z'] >= zlim[0]) & (catalog['Z'] < zlim[1]) & (catalog['ABSMAG_R'] >= maglim[0]) & (catalog['ABSMAG_R'] < maglim[1])
    if 'bitwise' in weight_type and 'default' in weight_type:
        mask &= (catalog['FRAC_TLOBS_TILES'] != 0)
    if option:
        if 'noNorth' in option:
            decmask = catalog['DEC'] < 32.375
            mask &= decmask
            
        if 'noDES' in option:
            #from LSS import common_tools as common
            isDES = common.select_regressis_DES(catalog)
            mask &= ~isDES
        
        if 'elgzmask' in option:
            zmask = ((catalog['Z'] >= 1.49) & (catalog['Z'] < 1.52))
            mask &= ~zmask

        if 'elgzcatas' in option:
            zmask = ((catalog['Z'] >= 1.31) & (catalog['Z'] < 1.33))
            mask &= ~zmask
            
        if 'ntile' in option:
            if '=' in option:
                opsp = option.split('=')
                nt = int(opsp[1])
                mask &= catalog['NTILE'] == nt
            if '>' in option:
                opsp = option.split('>')
                nt = int(opsp[1])
                mask &= catalog['NTILE'] >= nt
        #if 'zcmb' in option:
        #    #fix to do random properly
        #    from LSS import common_tools as common
        #    zcmb = common.get_zcmbdipole(catalog['RA'],catalog['DEC'])
        #    newz = (1+catalog['Z'])*(1+zcmb)-1
        #    catalog['Z'] = newz
    logger.info('Using {:d} rows for {}.'.format(mask.sum(), name))
    positions = [catalog['RA'][mask], catalog['DEC'][mask], distance(catalog['Z'][mask])]
    weights = np.ones_like(positions[0])

    if 'completeness_only' in weight_type and 'bitwise' in weight_type:
        raise ValueError('inconsistent choices were put into weight_type')

    #if name == 'data':
    if 'zfail' in weight_type:
        weights *= catalog['WEIGHT_ZFAIL'][mask]
        print('multiplying weights by WEIGHT_ZFAIL')
    if 'default' in weight_type and 'bitwise' not in weight_type:
        weights *= catalog['WEIGHT'][mask]
        print('multiplying weights by WEIGHT')     
    if 'removeIMSYS' in weight_type:
        #assumes default already added the rest of the weights and that SYS was used as default weight
        weights /= catalog['WEIGHT_SYS'][mask]
        print('dividing weights by WEIGHT_SYS') 
    #if 'RF' in weight_type:
    #    weights *= catalog['WEIGHT_RF'][mask]
    #    print('multiplying weights by WEIGHT_RF')
    #if 'SN' in weight_type:
    #    weights *= catalog['WEIGHT_SN'][mask]
    #    print('multiplying weights by WEIGHT_SN')
    if 'swapinLIN' in weight_type or 'IMLIN_FID' in weight_type:
        #assumes default already added the rest of the weights and that SN was used as default weight
        weights *=  catalog['WEIGHT_IMLIN'][mask]/catalog['WEIGHT_SYS'][mask]


    if 'IMLIN_FINEZBIN_ALLEBVCMB' in weight_type:
        weights *=  catalog['WEIGHT_IMLIN_FINEZBIN_ALLEBVCMB'][mask]/catalog['WEIGHT_SYS'][mask]
    
    elif 'IMLIN_FINEZBIN_ALLEBV' in weight_type:
        weights *=  catalog['WEIGHT_IMLIN_FINEZBIN_ALLEBV'][mask]/catalog['WEIGHT_SYS'][mask]

    elif 'IMLIN_FINEZBIN' in weight_type:
        weights *=  catalog['WEIGHT_IMLIN_FINEZBIN'][mask]/catalog['WEIGHT_SYS'][mask]

    if 'IMLIN_ALLEBVCMB' in weight_type:
        weights *=  catalog['WEIGHT_IMLIN_ALLEBVCMB'][mask]/catalog['WEIGHT_SYS'][mask]

    elif 'IMLIN_ALLEBV' in weight_type:
        weights *=  catalog['WEIGHT_IMLIN_ALLEBV'][mask]/catalog['WEIGHT_SYS'][mask]

    elif 'IMLIN_ALL' in weight_type:
        weights *=  catalog['WEIGHT_IMLIN_ALL'][mask]/catalog['WEIGHT_SYS'][mask]

    if 'IMLIN_NODEBV' in weight_type:
        weights *=  catalog['WEIGHT_IMLIN_NODEBV'][mask]/catalog['WEIGHT_SYS'][mask]


    if 'IMLIN_1ZBIN' in weight_type:
        weights *=  catalog['WEIGHT_IMLIN_1ZBIN'][mask]/catalog['WEIGHT_SYS'][mask]


    if 'swapinREGLIN' in weight_type:
        #assumes default already added the rest of the weights and that SN was used as default weight
        weights *=  catalog['WEIGHT_Linear'][mask]/catalog['WEIGHT_IMLIN'][mask]

    if 'swapinLNCLUS' in weight_type:
        #assumes default already added the rest of the weights and that SN was used as default weight
        weights *=  catalog['WEIGHT_IMLIN_CLUS'][mask]/catalog['WEIGHT_SYS'][mask]



    if 'swapinRF' in weight_type:
        #assumes default already added the rest of the weights and that SN was used as default weight
        weights *=  catalog['WEIGHT_RF'][mask]/catalog['WEIGHT_SN'][mask]
    if 'removeSN' in weight_type:
        #assumes default already added the rest of the weights and that SN was used as default weight
        weights /=  catalog['WEIGHT_SN'][mask]

    if 'addRF' in weight_type:
        #assumes no imaging systematic weights were in default
        weights *=  catalog['WEIGHT_RF'][mask]
    if 'addSN' in weight_type:
        #assumes no imaging systematic weights were in default
        weights *=  catalog['WEIGHT_SN'][mask]
    if 'addIMLIN' in weight_type:
        #assumes no imaging systematic weights were in default
        weights *=  catalog['WEIGHT_IMLIN'][mask]

    if 'data_compup' in weight_type and name == 'data':
        weights = catalog['WEIGHT_COMP'][mask]/catalog['FRAC_TLOBS_TILES'][mask]
    if 'completeness_only' in weight_type:
        weights = catalog['WEIGHT_COMP'][mask]
        print('weights set to WEIGHT_COMP')
    #if 'EB' in weight_type:
    #    weights *=  catalog['WEIGHT_SYSEB'][mask]
    #    print('multiplying weights by WEIGHT_SYSEB')
    if 'nofail' in weight_type:
        weights /= catalog['WEIGHT_ZFAIL'][mask]
        print('dividing weights by WEIGHT_ZFAIL')
    if 'addGFLUX' in weight_type:
        weights *= catalog['WEIGHT_FIBERFLUX'][mask]
        print('multiplying weights by WEIGHT_FIBERFLUX')
    if 'addSSR' in weight_type:
        weights *= catalog['WEIGHT_focal'][mask]
        print('multiplying weights by WEIGHT_focal')
    if 'FKP' in weight_type:
        logger.info('P0 when applying is '+str(P0))
        if P0 is None:
            weights *= catalog['WEIGHT_FKP'][mask]
            logger.info('multiplying weights by WEIGHT_FKP')
        else:
            P0 = float(P0)
            weights *= 1/(1+catalog['NX'][mask]*P0)
            logger.info('multiplied weights wfkp using custom P0 '+str(P0))
    if name == 'data' and 'bitwise' in weight_type:
        if 'default' in weight_type:
            #weights /= catalog['WEIGHT_COMP'][mask]
            weights = catalog['WEIGHT_SYS'][mask]*catalog['WEIGHT_ZFAIL'][mask]#/(129/(1+128*catalog['PROB_OBS'][mask]))
            

            #print('dividing weights by WEIGHT_COMP')
        if 'pip' in weight_type:
            weights = catalog['WEIGHT'][mask]/catalog['WEIGHT_COMP'][mask]
        if 'FKP' in weight_type:
            weights *= catalog['WEIGHT_FKP'][mask]
            print('multiplying weights by WEIGHT_FKP')

        weights = _format_bitweights(catalog['BITWEIGHTS'][mask]) + [weights]
    if name == 'data' and 'IIP' in weight_type:
        weights = 129/(1+128*catalog['PROB_OBS'][mask])
        if 'default' in weight_type:
            #weights /= catalog['WEIGHT_COMP'][mask]
            weights *= catalog['WEIGHT_SYS'][mask]*catalog['WEIGHT_ZFAIL'][mask]
        if 'FKP' in weight_type:
            weights *= catalog['WEIGHT_FKP'][mask]
            print('multiplying weights by WEIGHT_FKP')

    if name == 'randoms':
        #if 'default' in weight_type:
        #    weights *= catalog['WEIGHT'][mask]
        if 'bitwise' in weight_type:# and 'default' in weight_type:
            if 'default' in weight_type:
                weights = np.ones_like(positions[0])#catalog['WEIGHT_SYS'][mask]*catalog['WEIGHT_ZFAIL'][mask]
                logger.info('all random weights set to 1')
            if 'pip' in weight_type:
                weights = catalog['WEIGHT'][mask]
            if 'FKP' in weight_type:
                weights *= catalog['WEIGHT_FKP'][mask]
                print('multiplying weights by WEIGHT_FKP')

        if 'IIP' in weight_type and 'default' in weight_type:
            weights = np.ones_like(positions[0])#catalog['WEIGHT_SYS'][mask]*catalog['WEIGHT_ZFAIL'][mask]
            logger.info('all random weights set to 1')
            if 'FKP' in weight_type:
                weights *= catalog['WEIGHT_FKP'][mask]
                print('multiplying weights by WEIGHT_FKP')

            #weights /= catalog['FRAC_TLOBS_TILES'][mask]
            #print('dividing weights by FRAC_TLOBS_TILES')
        if 'NTMP' in weight_type:
            if fac_ntmp is not None:
                wts = common.apply_wntmp(catalog['NTILE'][mask], fac_ntmp[0], fac_ntmp[1])
                weights *= wts[0]
                logger.info('multiplied randoms by NTMP weights')
            else:
                logger.info('fac_ntmp was None, so nothing happening with it')

#         if 'RF' in weight_type:
#             weights *= catalog['WEIGHT_RF'][mask]*catalog['WEIGHT_COMP'][mask]
#         if 'zfail' in weight_type:
#             weights *= catalog['WEIGHT_ZFAIL'][mask]
#         if 'completeness_only' in weight_type:
#             weights = catalog['WEIGHT_COMP'][mask]
#         if 'EB' in weight_type:
#             weights *=  catalog['WEIGHT_SYSEB'][mask]*catalog['WEIGHT_COMP'][mask]   
#         if 'FKP' in weight_type:
#             weights *= catalog['WEIGHT_FKP'][mask]
#         if 'nofail' in weight_type:
#             weights /= catalog['WEIGHT_ZFAIL'][mask]
#         if 'fluxfail' in weight_type:
#             weights *= (catalog['WEIGHT_ZFAIL_FIBERFLUX'][mask]/catalog['WEIGHT_ZFAIL'][mask])



    if return_mask:
        return positions, weights, mask
    return positions, weights


def _concatenate(arrays):
    if isinstance(arrays[0], (tuple, list)):  # e.g., list of bitwise weights for first catalog
        array = [np.concatenate([arr[iarr] for arr in arrays], axis=0) for iarr in range(len(arrays[0]))]
    else:
        array = np.concatenate(arrays, axis=0)  # e.g. individual weights for first catalog
    return array


def read_clustering_positions_weights(distance, zlim =(0., np.inf), maglim=None, weight_type='default', name='data', concatenate=False, option=None, region=None, cat_read=None, dat_cat=None, ran_cat=None,P0=None, **kwargs):
    #print(kwargs)
    #if 'GC' in region:
    if type(region) is not list:
        region = [region]
    
    if cat_read == None:

        def read_positions_weights(name):
            cat_full = catalog_fn(ctype='full_HPmapcut', name='data', **kwargs)
#                    cat_full = catalog_fn(ctype='full', name=name, **kwargs)
            fac_ntmp =  None
            if name != 'data' and 'NTMP' in weight_type:
                getntmp = 'getntmp'
                logger.info('getting NTMP info')
                ff = fitsio.read(cat_full.replace('PIP/',''),columns=['BITWEIGHTS','PROB_OBS','LOCATION_ASSIGNED','NTILE'])
                fac_ntmp = common.compute_wntmp(ff['BITWEIGHTS'], ff['PROB_OBS'], ff['LOCATION_ASSIGNED'], ff['NTILE'])
                del ff
            else:
                logger.info(str(name)+' '+str(weight_type))

            positions, weights = [], []
            for reg in region:
                cat_fns = catalog_fn(ctype='clustering', name=name, region=reg, **kwargs)
                #if name=='data':
                logger.info('Loading {}.'.format(cat_fns))
                logger.info(name+' fac_ntmp is '+str(fac_ntmp))
                isscalar = not isinstance(cat_fns, (tuple, list))
   
                
                if isscalar:
                    cat_fns = [cat_fns]
                if name=='data':
                    def _get_tab(cat_fn):
                        if '.fits' in cat_fn:
                            tab = Table.read(cat_fn)
                        if '.h5' in cat_fn:
                            tab = common.read_hdf5_blosc(cat_fn)
                        if 'bitwise' in weight_type:
                            if 'BITWEIGHTS' in list(tab.dtype.names):
                                pass
                            else:   
                                tab = join(tab, Table.read(cat_full)['TARGETID', 'BITWEIGHTS'], keys='TARGETID', join_type='left')
                        #else:
                        if option is not None:
                            if 'RSDZ' in option:
                                if 'Z' in list(tab.dtype.names):
                                    tab.remove_column('Z')
                                tab.rename_column('RSDZ', 'Z')    
                        return tab
                    positions_weights = [get_clustering_positions_weights(_get_tab(cat_fn), distance, zlim=zlim, maglim=maglim, weight_type=weight_type, name=name, option=option,P0=P0) for cat_fn in cat_fns]
                else:
                    if '.fits' in cat_fns[0]:
                        read_func = Table.read
                    elif '.h5' in cat_fns[0]:
                        read_func = common.read_hdf5_blosc
                    positions_weights = [get_clustering_positions_weights(read_func(cat_fn), distance, zlim=zlim, maglim=maglim, weight_type=weight_type, name=name, option=option,P0=P0,fac_ntmp=fac_ntmp) for cat_fn in cat_fns]
                
                if isscalar:
                    positions.append(positions_weights[0][0])
                    weights.append(positions_weights[0][1])
                else:
                    p, w = [tmp[0] for tmp in positions_weights], [tmp[1] for tmp in positions_weights]
                    if concatenate:
                        p, w = _concatenate(p), _concatenate(w)
                    positions.append(p)
                    weights.append(w)
            
            return positions, weights

    if cat_read != None:
        def read_positions_weights(name):
            positions, weights = [], []
            for reg in region:
                logger.info('Using arrays.')
                
                if name == 'data':
                    cat_read = dat_cat
                if name == 'randoms':
                    cat_read = ran_cat
                   
                    
                positions_weights = [get_clustering_positions_weights(cat_read, distance, zlim=zlim, maglim=maglim, weight_type=weight_type, name=name, option=option,P0=P0)]
                if name == 'data':
                    positions.append(positions_weights[0][0])
                    weights.append(positions_weights[0][1])
                
                if name == 'randoms':
                    p, w = [tmp[0] for tmp in positions_weights], [tmp[1] for tmp in positions_weights]
                    positions.append(p)
                    weights.append(w)
            
            return positions, weights
        
    
    if isinstance(name, (tuple, list)):
        return [read_positions_weights(n) for n in name]
    return read_positions_weights(name)


def get_full_positions_weights(catalog, name='data', weight_type='default', fibered=False, region='', return_mask=False, weight_attrs=None):
    
    from pycorr.twopoint_counter import get_inverse_probability_weight
    if weight_attrs is None: weight_attrs = {}
    mask = np.ones(len(catalog), dtype='?')
    #if region in ['DS', 'DN']:
    if region != 'full':
        mask &= select_region(catalog['RA'], catalog['DEC'], region)
    #elif region:
    #    #mask &= catalog['PHOTSYS'] == region.strip('_')
    #    mask &= catalog['PHOTSYS'] == region.strip('GC')
        

    if fibered: mask &= catalog['LOCATION_ASSIGNED']
    positions = [catalog['RA'][mask], catalog['DEC'][mask], catalog['DEC'][mask]]
    weights = np.ones_like(positions[0])
    if 'pip' in weight_type:
        weights = catalog['WEIGHT_NTILE'][mask]
    if 'FKP' in weight_type:
        weights *= catalog['WEIGHT_FKP_NTILE'][mask]

    if name == 'data' and fibered:
        #if 'default' in weight_type or 'completeness' in weight_type:
        #    weights = get_inverse_probability_weight(_format_bitweights(catalog['BITWEIGHTS'][mask]), **weight_attrs)
        if 'bitwise' in weight_type:
            weights = _format_bitweights(catalog['BITWEIGHTS'][mask])+[weights]
    #else: weights = np.ones_like(positions[0])
    if return_mask:
        return positions, weights, mask
    return positions, weights


def read_full_positions_weights(name='data', weight_type='default', fibered=False, region='', weight_attrs=None, **kwargs):
    #if 'GC' in region:
    #    region = [region]
    if type(region) is not list:
        region = [region]

    def read_positions_weights(name):
        positions, weights = [], []
        for reg in region:
            cat_fn = catalog_fn(ctype='full_HPmapcut', name=name, **kwargs)
            #cat_fn = catalog_fn(ctype='full', name=name, **kwargs)
            logger.info('Loading {}.'.format(cat_fn))
            if isinstance(cat_fn, (tuple, list)):
                #catalog = vstack([Table.read(fn.replace('PIP/','')) for fn in cat_fn])
                catalog = vstack([Table.read(fn) for fn in cat_fn])
            else:
                #catalog = Table.read(cat_fn.replace('PIP/',''))
                catalog = Table.read(cat_fn)
            p, w = get_full_positions_weights(catalog, name=name, weight_type=weight_type, fibered=fibered, region=reg, weight_attrs=weight_attrs)
            positions.append(p)
            weights.append(w)
            if fibered:
                logger.info('loaded fibered full for '+name + ' for region '+reg + ' with '+str(len(p))+' entries')
            else:
                logger.info('loaded parent full for '+name+ ' for region '+reg + ' with '+str(len(p))+' entries')    
        return positions, weights

    if isinstance(name, (tuple, list)):
        return [read_positions_weights(n) for n in name]
    return read_positions_weights(name)


def normalize_data_randoms_weights(data_weights, randoms_weights, weight_attrs=None):
    # Renormalize randoms / data for each input catalogs
    # data_weights should be a list (for each N/S catalogs) of weights
    import inspect
    from pycorr.twopoint_counter import _format_weights, get_inverse_probability_weight
    if weight_attrs is None: weight_attrs = {}
    weight_attrs = {k: v for k, v in weight_attrs.items() if k in inspect.getargspec(get_inverse_probability_weight).args}
    wsums, weights = {}, {}
    for name, catalog_weights in zip(['data', 'randoms'], [data_weights, randoms_weights]):
        wsums[name], weights[name] = [], []
        for w in catalog_weights:
            w, nbits = _format_weights(w, copy=True)  # this will sort bitwise weights first, then single individual weight
            iip = get_inverse_probability_weight(w[:nbits], **weight_attrs) if nbits else 1.
            iip = iip * w[nbits]
            wsums[name].append(iip.sum())
            weights[name].append(w)
    wsum_data, wsum_randoms = sum(wsums['data']), sum(wsums['randoms'])
    for icat, w in enumerate(weights['randoms']):
        factor = wsums['data'][icat] / wsums['randoms'][icat] * wsum_randoms / wsum_data
        w[-1] *= factor
        logger.info('Rescaling randoms weights of catalog {:d} by {:.4f}.'.format(icat, factor))
    return weights['data'], weights['randoms']


def concatenate_data_randoms(data, randoms=None, **kwargs):

    if randoms is None:
        positions, weights = data
        return _concatenate(positions), _concatenate(weights)

    positions, weights = {}, {}
    for name in ['data', 'randoms']:
        positions[name], weights[name] = locals()[name]
    for name in positions:
        concatenated = not isinstance(positions[name][0][0], (tuple, list))  # first catalog, unconcatenated [RA, DEC, distance] (False) or concatenated RA (True)?
        if concatenated:
            positions[name] = _concatenate(positions[name])
        else: 
            positions[name] = [_concatenate([p[i] for p in positions[name]]) for i in range(len(positions['randoms'][0]))]
    data_weights, randoms_weights = [], []
    if concatenated:
        wd, wr = normalize_data_randoms_weights(weights['data'], weights['randoms'], weight_attrs=kwargs.get('weight_attrs', None))
        weights['data'], weights['randoms'] = _concatenate(wd), _concatenate(wr)
    else:
        for i in range(len(weights['randoms'][0])):
            wd, wr = normalize_data_randoms_weights(weights['data'], [w[i] for w in weights['randoms']], weight_attrs=kwargs.get('weight_attrs', None))
            data_weights.append(_concatenate(wd))
            randoms_weights.append(_concatenate(wr))
        weights['data'] = data_weights[0]
        for wd in data_weights[1:]:
            for w0, w in zip(weights['data'], wd): assert np.all(w == w0)
        weights['randoms'] = randoms_weights
    return [(positions[name], weights[name]) for name in ['data', 'randoms']] 


