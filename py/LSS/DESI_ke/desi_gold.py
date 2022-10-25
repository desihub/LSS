import os
import runtime
import numpy as np
import astropy.units as u
import argparse

from   config              import Configuration
from   findfile            import findfile
from   astropy.coordinates import SkyCoord
from   astropy.table       import Table, vstack, hstack, unique, join
from   ros_tools           import tile2rosette, calc_rosr, ros_limits
from   gama_limits         import gama_field
from   cartesian           import cartesian, rotate
from   cosmo               import cosmo, distmod
from   lss                 import fetch_lss
from   bitmask             import lumfn_mask


def desi_gold(args, survey='sv3', release='fuji'):
    from   desiutil.dust                 import mwdust_transmission
    from   desitarget.sv3.sv3_targetmask import desi_mask, bgs_mask


    dryrun   = args.dryrun

    if survey != 'sv3':
        raise  NotImplementedError('Only SV3 is currently supported.')

    releases = {'sv3': 'fuji', 'main': 'guadalupe'}
    release  = releases[survey] 
    
    root     = os.environ['DESI_ROOT'] + f'/spectro/redux/{release}/healpix/'
    fpath    = root + 'tilepix.fits'

    print(f'Fetching {fpath}')

    tpix     = Table.read(fpath)

    # Note:  max. Bright SV3 tileid is 595.
    tiles    = np.arange(1000)
    ros      = np.array([tile2rosette(x) for x in tiles])
    
    # https://desi.lbl.gov/trac/wiki/SurveyOps/OnePercent
    # G12: [1,2]; G15: [8,9,10, 17]

    # Unique rosettes in SV3.
    uros   = np.unique(ros)
    uros   = uros[uros > -1]
    
    gama   = np.isin(ros, uros)

    tiles  = tiles[gama]
    
    tpix   = tpix[np.isin(tpix['TILEID'].data, tiles)]
    hps    = np.unique(tpix['HEALPIX'].data)

    root  += '/{}/bright/'.format(survey)

    fpaths = [root + '{}/{}/redrock-{}-bright-{}.fits'.format(str(x)[:3], x, survey, x) for x in hps]
    fpaths = [x for x in fpaths if os.path.exists(x)]

    print('Fetching {}'.format(fpaths[0]))
    
    # e.g. 280/28027/redrock-sv3-bright-28027.fits
    tabs   = []

    print('Gathering DESI zs.')

    for x in fpaths:
        zbest = Table.read(x, hdu='REDSHIFTS')
        fmap  = Table.read(x, hdu='FIBERMAP')
        efmap = Table.read(x, hdu='EXP_FIBERMAP')

        # row ordered.                                                                                                                                                                                     
        assert  np.all(zbest['TARGETID'] == fmap['TARGETID'])
        
        tids  = np.unique(zbest['TARGETID'])

        efmap['ROS'] = tile2rosette(efmap['TILEID'].data)
    
        efmap_tid = efmap['TARGETID', 'ROS'] 
        efmap_tid = unique(efmap_tid, keys=['TARGETID'])
        efmap_tid.sort('TARGETID')

        assert  np.all(efmap_tid['TARGETID'].data == tids)
    
        del  fmap['TARGETID']
        
        zbest = hstack([zbest, fmap])
    
        zbest = join(zbest, efmap_tid, join_type='left', keys='TARGETID')
        zbest['ROS_DIST'] = 1.e4 * np.ones_like(zbest['Z'])

        # Distance to the closest rosette center. 
        for rosn in np.unique(efmap['ROS'].data):
            # TODO: Small rosette overlap.                                                                                                                                                         
            isin  = zbest['ROS'].data == rosn
            zbest['ROS_DIST'][isin] = calc_rosr(rosn, zbest['TARGET_RA'].data[isin], zbest['TARGET_DEC'].data[isin])
    
        tabs.append(zbest)
    
    desi_zs = vstack(tabs)

    # remove skies.                                                                                                                                                                                 
    print('Sky frac: 1-{:.6f}'.format(np.mean(desi_zs['TARGETID'].data >= 0)))

    desi_zs = desi_zs[desi_zs['TARGETID'].data >= 0]    
    desi_zs = desi_zs[(desi_zs['SV3_BGS_TARGET'].data & bgs_mask['BGS_BRIGHT']) != 0]

    desi_zs['GOOD_Z'] = (desi_zs['ZWARN'].data == 0) & (desi_zs['DELTACHI2'].data > 40)

    desi_zs['ZDESI'] = desi_zs['Z']

    del  desi_zs['Z']

    desi_zs['FIELD'] = [f'R{xx}' for xx in desi_zs['ROS'].data]

    ##  Cut DESI to good redshifts.                                                                                                                                                                
    desi_zs['GAMA_FIELD'] = gama_field(desi_zs['TARGET_RA'].data, desi_zs['TARGET_DEC'].data)

    ##  TODO/HACK: PHOTSYS ASSUMED S.                                                                                                                                                                     
    desi_zs['GMAG_DRED']   = 22.5 - 2.5 * np.log10(desi_zs['FLUX_G'].data / mwdust_transmission(desi_zs['EBV'].data, 'G', 'S', match_legacy_surveys=True))
    desi_zs['RMAG_DRED']   = 22.5 - 2.5 * np.log10(desi_zs['FLUX_R'].data / mwdust_transmission(desi_zs['EBV'].data, 'R', 'S', match_legacy_surveys=True))

    desi_zs['W1MAG_DRED']  = 22.5 - 2.5 * np.log10(desi_zs['FLUX_W1'].data / mwdust_transmission(desi_zs['EBV'].data, 'W1', 'S', match_legacy_surveys=True))
    desi_zs['W2MAG_DRED']  = 22.5 - 2.5 * np.log10(desi_zs['FLUX_W2'].data / mwdust_transmission(desi_zs['EBV'].data, 'W2', 'S', match_legacy_surveys=True))
    
    desi_zs['GMR']        = desi_zs['GMAG_DRED'] - desi_zs['RMAG_DRED']
    desi_zs['DETMAG']     = desi_zs['RMAG_DRED']
    
    desi_zs.meta['PET_OFFSET'] = 0.12
    desi_zs['LEGACYPET']       = desi_zs['RMAG_DRED'] + desi_zs.meta['PET_OFFSET']
    
    desi_zs['IN_GOLD']    = desi_zs['GOOD_Z'].data & (desi_zs['ZDESI'] > 0.039)  & (desi_zs['ZDESI'] < 0.263)

    clustering, full      = fetch_lss(pprint=False, sort=False) 

    clustering_ids        = np.unique(clustering['TARGETID'].data)
    full_ids              = np.unique(full['TARGETID'].data)

    desi_zs['IN_CLUSTERING'] = np.isin(desi_zs['TARGETID'].data, clustering_ids)
    desi_zs['IN_FULL']       = np.isin(desi_zs['TARGETID'].data, full_ids)

    clustering_cols = ['WEIGHT_ZFAIL', 'WEIGHT', 'NZ']
    full_cols       = ['NTILE', 'TILES', 'TILELOCIDS', 'LOCATION_ASSIGNED', 'TILELOCID_ASSIGNED', 'COMP_TILE', 'FRACZ_TILELOCID', 'BITWEIGHTS', 'PROB_OBS']
                       
    for cols, cat in zip([clustering_cols, full_cols], [clustering, full]):
        cols       += ['TARGETID']

        # TODO/HACK multiple targetid occurences in clustering catalog?
        desi_zs     = join(desi_zs, cat[cols], keys='TARGETID')
                       
    ##  Archive step. 
    archive = Table(desi_zs, copy=True)

    xyz     = cartesian(desi_zs['TARGET_RA'].data, desi_zs['TARGET_DEC'].data, desi_zs['ZDESI'].data)

    desi_zs['CARTESIAN_X'] = xyz[:,0]
    desi_zs['CARTESIAN_Y'] = xyz[:,1]
    desi_zs['CARTESIAN_Z'] = xyz[:,2]

    xyz     = rotate(desi_zs['TARGET_RA'].data, desi_zs['TARGET_DEC'].data, xyz)

    desi_zs['ROTCARTESIAN_X'] = xyz[:,0]
    desi_zs['ROTCARTESIAN_Y'] = xyz[:,1]
    desi_zs['ROTCARTESIAN_Z'] = xyz[:,2]

    desi_zs['LUMDIST']        = cosmo.luminosity_distance(desi_zs['ZDESI'].data)
    desi_zs['DISTMOD']        = distmod(desi_zs['ZDESI'].data)

    desi_zs['IN_D8LUMFN']     = np.zeros_like(desi_zs['FIELD'], dtype=int)
    desi_zs['CONSERVATIVE']   = np.zeros_like(desi_zs['FIELD'], dtype=int)

    desi_zs.meta['IMMUTABLE'] = 'TRUE'
    
    desi_zs.pprint()

    survey = 'desi'
    '''
    fpath  = findfile(ftype='gold', dryrun=False, survey=survey)
    opath  = fpath.replace('desi_gold', 'desi_sv3_gold')

    print('Writing {}'.format(opath))

    desi_zs.write(opath, format='fits', overwrite=True)
    
    ##  ----  GAMA GOLD
    gold  = Table.read(fpath)

    del  gold['CARTESIAN_X']
    del  gold['CARTESIAN_Y']
    del  gold['CARTESIAN_Z']
    del  gold['ROTCARTESIAN_X']
    del  gold['ROTCARTESIAN_Y']
    del  gold['ROTCARTESIAN_Z']
    del  gold['DETMAG']
    del  gold['DISTMOD']
    del  gold['LUMDIST']
    del  gold['GMR']

    gold.pprint()

    # DESI
    uros          = [1,2,8,9,10,17] 
    desi_zs       = desi_zs[np.isin(desi_zs['ROS'], uros)]
    archive       = archive[np.isin(archive['ROS'], uros)]
    
    c             = SkyCoord(ra=desi_zs['TARGET_RA']*u.degree, dec=desi_zs['TARGET_DEC']*u.degree)

    print('Matching DESI to GAMA Gold.')

    # GAMA
    catalog       = SkyCoord(ra=gold['RA'], dec=gold['DEC'])
    idx, d2d, d3d = c.match_to_catalog_3d(catalog)
    
    # Now idx are indices into catalog that are the closest objects to each of the coordinates in c, d2d are the on-sky distances between them
    gold_match   = gold[idx]

    del  gold_match['FIELD']

    to_join               = Table(gold_match, copy=True) 
    to_join.meta          = {}
    
    desi_zs               = hstack([desi_zs, to_join])
    desi_zs['GAMA_SEP']   = d2d.to(u.arcsec)
    desi_zs['GOOD_MATCH'] = desi_zs['GAMA_SEP'] < 0.3
    
    desi_zs.meta['IMMUTABLE']  = 'TRUE'

    max_sep      = 0.5 * u.arcsec

    print('Fraction desi matched to gold at 0.5 arcseconds: {:.6f}'.format(np.mean(desi_zs['GAMA_SEP'] < max_sep)))

    opath = fpath.replace('desi_gold', 'desi_gama')
    
    print('Writing {}'.format(opath))
    
    desi_zs.write(opath, format='fits', overwrite=True)
    
    ## --------------------
    idx, d2d, d3d = catalog.match_to_catalog_3d(c)
    
    desi_match    = archive[idx]
    
    to_join       = Table(desi_match, copy=True)
    to_join.meta  = {}

    gold          = hstack([gold, to_join])
    gold['DESI_SEP'] = d2d.to(u.arcsec)
    gold['GOOD_MATCH'] = gold['DESI_SEP'] < 0.3
    
    gold['ROS_DIST'] = 1.e99
    
    for rosn in uros:
        new_dist = calc_rosr(rosn, gold['RA'].data, gold['DEC'].data)
        
        gold['ROS_DIST'] = np.minimum(gold['ROS_DIST'].data, new_dist)
    
    max_sep = 0.5 * u.arcsec

    print('Fraction gold matched to desi at 0.5 arcseconds: {:.6f}'.format(np.mean(gold['DESI_SEP'] < max_sep)))

    gold.meta['IMMUTABLE'] = 'TRUE'

    opath = fpath.replace('desi_gold', 'gama_desi')
    
    print('Writing {}'.format(opath))
    
    gold.write(opath, format='fits', overwrite=True)
    '''

    in_gold                   =  desi_zs['GOOD_Z'].data & (desi_zs['ZDESI'] > 0.039)  & (desi_zs['ZDESI'] < 0.263)

    # No cut to GAMA rosettes.
    # in_gold                &=  np.isin(desi_zs['ROS'].data, [1,2,8,9,10,17])
    
    desi_zs                   = desi_zs[in_gold]
    desi_zs['RA']             = desi_zs['TARGET_RA']
    desi_zs['DEC']            = desi_zs['TARGET_DEC']
    desi_zs['ZSURV']          = desi_zs['ZDESI']
    desi_zs['DETMAG']         = desi_zs['RMAG_DRED']
    desi_zs['DISTMOD']        = distmod(desi_zs['ZDESI'].data)

    limits                    = ros_limits(dryrun)

    hi_comp                   = (desi_zs['ROS_DIST'].data > limits[0]) & (desi_zs['ROS_DIST'].data < limits[1])
    area                      = np.pi * (limits[1]**2. - limits[0]**2.)

    desi_zs['IN_D8LUMFN']    += ~hi_comp * lumfn_mask.DESI_HICOMP

    desi_zs                   = desi_zs[desi_zs['IN_D8LUMFN'].data == 0]

    desi_zs.meta['AREA']       = area * len(np.unique(desi_zs['FIELD'].data))
    desi_zs.meta['IMMUTABLE']  = 'TRUE'
    desi_zs.meta['RELEASE']    = release.upper()
    
    desi_zs.meta['RLIM']       = 19.5
    desi_zs.meta['RMAX']       = 12.0
    desi_zs.meta['MAX_SEP']    = 10.0 # Expected max. angular separation in a rosette.
    
    opath                      = findfile(ftype='gold', dryrun=dryrun, survey=survey)

    print('Writing {}'.format(opath))

    desi_zs.write(opath, format='fits', overwrite=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Gen desi gold cat.')
    parser.add_argument('--log',          help='Create a log file of stdout.', action='store_true')
    parser.add_argument('--config',       help='Path to configuration file', type=str, default=findfile('config'))
    parser.add_argument('--dryrun',       help='Dryrun of 5k galaxies', action='store_true')
    parser.add_argument('--nooverwrite',  help='Do not overwrite outputs if on disk', action='store_true')

    args   = parser.parse_args()
    '''
    config = Configuration(args.config)
    config.update_attributes('gold', args)
    config.write()
    '''
    desi_gold(args)

    print('Done.')


