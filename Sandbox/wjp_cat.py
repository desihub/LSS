# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
=========================
lsscatalog
=========================
Output LSS catalog, child of a zcatalog.fits, a mtl.fits and
the fiber assignment output.

WJP data model.
https://docs.google.com/document/d/1GNqJUvAugIhsVI7gpeZHOCcpXp4DkaOtrHqRHr_Cgt8/edit#

By MJW (Dec. 28 2019)

https://portal.nersc.gov/project/desi/users/mjwilson/lss_cat.py
"""

import os
import fitsio
import numpy                as      np
import astropy.io.fits      as      fits
import astropy.table        as      table

from   astropy.table        import  Table, Column
from   desiutil.log         import  get_logger, DEBUG
from   desimodel.footprint  import  is_point_in_desi


MAX_NFIBER  = 10
MAX_NTILE   = 10

log         = get_logger()

def blind_distances(final):
    '''
    Given an LSS catalog instance with Z filled, 
    return the same instance with (blinded) cosmo
    distances. 

    Parameters                                                                                                                                                                                                                                                                                                                                                                                                 
    ----------                                                                                                                                                                                                                                                                                                                                                                                                 
    final:  :class:`astropy.table.Table` 
        LSS catalog Table.                                                                                                                                                                                                                                                                                                                                                                                                                
    Returns                                                                                                                                                                                                                                                                                                                                                                                                    
    -------                                                                                                                                                                                                                                                                                                                                                                                                    
    lsscatalog: :class:`astropy.table.Table`                                                                                                                                                                                                                                                                                                                                                                   
        LSS catalog Table. 
    '''

    from astropy.cosmology import WMAP9 as cosmo

    
    isin                          = final['Z'] > 0. 
    final['COSMO_BLINDCHI'][isin] = (cosmo.H(0) / 100.) * cosmo.comoving_distance(final['Z'][isin])
        
    return  final
    
def good_tile(tiles):
    # Assume a tile is good if an exposure achieves a cumulative SNR2 > 1.0
    is_good = tiles['SNR2'] > 1.0

    return  is_good

def datestamped_tiles(exps, MJD0, MJD1, printit=False):
    '''                                                                                                                                                                                                                                                                                                                                                                                                       
    Given an two MJDs and an exposures list, 
    returns a Table of the "good" tiles observed
    between these dates, together with a header. 

    Parameters                                                                                                                                                                                                                                                                                                                                                                                        
    ----------                                                                                                                                                                                                                                                                                                                                                                                            
                                                                                                                                                                                                                                                                                                                                                                                                           
    final:  :class:`astropy.table.Table`                                                                                                                                                                                                                                                                                                                                                                       
        LSS catalog Table.                                                                                                                                                                                                                                                                                                                                                                                   
                                                                                                                                                                                                                                                                                                                                                                                                          
    Returns                                                                                                                                                                                                                                                                                                                                                                                                   
    ----------
    _tiles: :class:`astropy.table.Table` 
    hdr:  : FITS table header.                                                                                                                                                                                                                                                                                                                                                                             
    -------                                                                                                                                                                                                                                                                                                                                                                                                 


    '''

    _tiles  = exps['TILEID', 'MJD', 'SNR2']
    
    # Remove spurious TILEIDs                                                                                                                                                                                      
    _tiles  = _tiles[_tiles['TILEID'] > -1]
              
    # Keep between input MJD limits. 
    _tiles  = _tiles[_tiles['MJD'] >= MJD0]
    _tiles  = _tiles[_tiles['MJD'] <= MJD1]

    # Define what tiles are acceptable.
    is_good = good_tile(_tiles)
    _tiles  = _tiles[is_good] 
    
    # Catch tiles for which N > 1 exposures were "good".                                                    
    _tiles  = table.unique(_tiles, keys='TILEID')
    
    # Remove MJD column. 
    del _tiles['MJD']
    del _tiles['SNR2']
   
    # Create header. 
    hdr         = fits.Header()

    hdr['MJD0'] = MJD0
    hdr['MJD1'] = MJD1
    
    if printit:
        print(_tiles)
    
    return  _tiles, hdr

def initialise_lsscatalog(final):
    final['Z']               =  -99.
    final['ZERR']            =  -99.
    final['ZWARN']           =  -99

    final['COSMO_BLINDCHI']  =  -99.
    final['COSMO_BLINDNZ']   =  -99.
    final['COSMO_BLINDWFKP'] =  -99.

    final['NGOOD_FIBERS']    =    0
    final['NGOOD_TILES']     =    0

    final['IMAGE_WGHT']      =  1.0
    final['SPEC_WGHT']       =  1.0

    final['ASSIGN_IIP']      =  1.0

    return  final

def lss_catalog(nobj=1):
    """
    Create an empty 'lsscatalog' table.
    
    Parameters
    ----------
    ntarget : :class:`int`
        Number of targets.

    Returns
    -------
    lsscatalog: :class:`astropy.table.Table`
        LSS catalog Table.    
    """
    from astropy.table import Table, Column

    # One row per target.
    lsscatalog = Table()
    
    lsscatalog.add_column(Column(name='TARGETID',        length=nobj, dtype='int64'))

    lsscatalog.add_column(Column(name='RA',              length=nobj, dtype='float64'))
    lsscatalog.add_column(Column(name='DEC',             length=nobj, dtype='float64'))
    
    lsscatalog.add_column(Column(name='DESI_TARGET',     length=nobj, dtype='int64'))
    lsscatalog.add_column(Column(name='BGS_TARGET',      length=nobj, dtype='int64'))
    
    # --  Targeting --
    lsscatalog.add_column(Column(name='IN_IMAGING',      length=nobj, dtype='i2'))    # INSIDE IMAGING FOOTPRINT
    lsscatalog.add_column(Column(name='IN_DESI',         length=nobj, dtype='i2'))    # INSIDE SPEC.   FOOTPRINT (DATE STAMPED ACCEPTABLE TILES)
    lsscatalog.add_column(Column(name='IN_COSMO',        length=nobj, dtype='int64')) # Inside cosmo footprint, e.g. EBV << 1.
    
    lsscatalog.add_column(Column(name='ANG_VETO_FLAG',   length=nobj, dtype='int64'))
    lsscatalog.add_column(Column(name='Z_VETO_FLAG',     length=nobj, dtype='int64'))
    
    # Assignment history.  Priority changs?
    lsscatalog.add_column(Column(name='PRIORITY_INIT',   length=nobj, dtype='float64'))
    lsscatalog.add_column(Column(name='SUBPRIORITY',     length=nobj, dtype='float64'))
        
    lsscatalog.add_column(Column(name='GOOD_FIBERS',     length=nobj, dtype='int64', shape=MAX_NFIBER))
    lsscatalog.add_column(Column(name='GOOD_TILES',      length=nobj, dtype='int64', shape=MAX_NTILE))

    # For each of GOOD_FIBERS, what was the PRIORITY_INIT of the target to which it was assigned.  
    lsscatalog.add_column(Column(name='FIBPRIORITY',     length=nobj, dtype='float64', shape=MAX_NFIBER))
    
    lsscatalog.add_column(Column(name='NGOOD_FIBERS',    length=nobj, dtype='int64'))
    lsscatalog.add_column(Column(name='NGOOD_TILES',     length=nobj, dtype='int64'))
    
    # Effective assignment completeness, non-pairwise.                                                                                                                                                                                 
    lsscatalog.add_column(Column(name='ASSIGN_IIP',      length=nobj, dtype='float64'))
    
    # --  Spectroscopic --
    lsscatalog.add_column(Column(name='Z',               length=nobj, dtype='float64'))
    lsscatalog.add_column(Column(name='ZERR',            length=nobj, dtype='float64'))
    lsscatalog.add_column(Column(name='ZWARN',           length=nobj, dtype='int64'))
    lsscatalog.add_column(Column(name='SPECTYPE',        length=nobj, dtype='S16'))

    lsscatalog.add_column(Column(name='COSMO_BLINDCHI',  length=nobj, dtype='float64'))
    lsscatalog.add_column(Column(name='COSMO_BLINDNZ',   length=nobj, dtype='float64'))
    lsscatalog.add_column(Column(name='COSMO_BLINDWFKP', length=nobj, dtype='float64'))
    
    # 1. / Imaging completeness [0., 1.]  
    lsscatalog.add_column(Column(name='IMAGE_WGHT',      length=nobj, dtype='float64'))

    # 1. / Spectroscopic completeness [0., 1.]
    lsscatalog.add_column(Column(name='SPEC_WGHT',      length=nobj, dtype='float64'))

    initialise_lsscatalog(lsscatalog)

    return  lsscatalog

def get_randassign(fpath, gen=True, mtl=None, randoms=None, lite=False):
    # Potential assignments for randoms.
    if gen:
        scratch         = os.environ['CSCRATCH']

        if not os.path.exists(scratch + '/desi/fiberassign/'):
          os.makedirs(scratch + '/desi/fiberassign/')
        
        if lite:
            cmd  = 'fiberassign --mtl {}'.format(scratch + '/desi/fiberassign/lite.fits')
            cmd += ' --sky {}'.format(scratch + '/desi/fiberassign/lite_skies.fits')
            cmd += ' --outdir {}/desi/fiberassign/'.format(scratch)
            cmd += ' --overwrite'

            print('System call:  {}'.format(cmd))

            os.system(cmd)

            return  None, None

        else:
            nretain         = len(mtl)

            rand_mtl        = Table(mtl, copy=True)
            rand_mtl['RA']  = randoms['RA'][:nretain]
            rand_mtl['DEC'] = randoms['DEC'][:nretain]
        
            rand_mtl.write(scratch + '/desi/fiberassign/rand_mtl.fits', format='fits', overwrite=True)
        
            cmd  = 'fiberassign --mtl {}'.format(scratch + '/desi/fiberassign/rand_mtl.fits')
            cmd += ' --sky /project/projectdirs/desi/target/catalogs/dr8/0.31.0/skies/skies-dr8-0.31.0.fits'        
            cmd += ' --outdir {}/desi/fiberassign/'.format(scratch)
                    
            print('System call:  {}'.format(cmd))
        
            os.system(cmd)

            return  None, None
        
    else:
        return  Table(fits.open(fpath)[4].data), Table(fits.open(fpath)[1].data)  

def pop_assigned(final, printit=False):
    ##  Match ID list to the date stamped tiles.
    ids               = ['072015']

    for id in ids:
        assigned        = Table(fits.open(root + 'fiberassign/tile-{}.fits'.format(id))[1].data)
        
        ##
        assignable      = Table(fits.open(root + 'fiberassign/tile-{}.fits'.format(id))[4].data)

        ##  Unique assignable targets in this tile.  Note:  uassigned, ** NOT ** unassined.
        uassigned, cnts = np.unique(assignable['TARGETID'], return_counts=True)
        repeated        = cnts > 1
        
        nassignable     = len(uassigned)
        ##  ntargetable = is_point_in_desi(tiles[tiles['TILEID'] == id], ...)

        ##  assign_cpt  = nassignable / ntargetable 
      
        print('Tile ID {} has {} and {} single target-fiber and many target-fiber respectively.'.format(id, len(repeated) - np.count_nonzero(repeated.astype(np.int)), np.count_nonzero(repeated.astype(np.int))))

      
        count           = 0
      
        # Set of targets reachable by a single fiber.
        for entry in uassigned[~repeated]:
            # This shouldn't be necessary ... 
            if entry in final['TARGETID']:
                print(entry)
            
                # Tiles.
                col                                                   = np.array(final[final['TARGETID'] == entry]['NGOOD_TILES'])[0]
            
                final['GOOD_TILES'][final['TARGETID'] == entry, col]  = id
                final['NGOOD_TILES'][final['TARGETID'] == entry]     += 1
          
                # Fibers.
                col                                                   = np.array(final[final['TARGETID'] == entry]['NGOOD_FIBERS'])[0] 
          
                final['GOOD_FIBERS'][final['TARGETID'] == entry, col] = assignable[assignable['TARGETID'] == entry]['FIBER']
                final['NGOOD_FIBERS'][final['TARGETID'] == entry]    += 1

                # Priority of the target to which the potential fiber was actually assigned. 
                final['FIBPRIORITY'][final['TARGETID'] == entry, col] = assigned[assigned == assignable[assignable['TARGETID'] == entry]['FIBER']]['PRIORITY_INIT']
                
                count                                                += 1

                if count > 25:
                    break

        count = 0
        
        # Set of targets reachable by multiple fibers, for this tile.
        for entry in uassigned[repeated]:
            # This shouldn't be necessary ...
            if entry in final['TARGETID']:
                print(entry)
            
                # Tiles.                                                                                                                                                                                                                     
                col                                                   = np.array(final[final['TARGETID'] == entry]['NGOOD_TILES'])[0]

                final['GOOD_TILES'][final['TARGETID'] == entry, col]  = id
                final['NGOOD_TILES'][final['TARGETID'] == entry]     += 1

                # Fibers.
                fentries    = assignable['FIBER'][assignable['TARGETID'] == entry]

                for fentry in fentries:
                    col                                               = np.array(final['NGOOD_FIBERS'][final['TARGETID'] == entry])[0]

                    final['GOOD_FIBERS'][final['TARGETID'] == entry, col] = fentry
                    final['NGOOD_FIBERS'][final['TARGETID'] == entry]    += 1

                    # Priority of the target to which the potential fiber was actually assigned.  Necessary for the ANG_VETO_FLAG = 3 class.
                    final['FIBPRIORITY'][final['TARGETID'] == entry, col] = assigned[assigned == assignable[assignable['TARGETID'] == entry]['FIBER']]['PRIORITY_INIT']
                
                    count                                                += 1

                    if count > 25:
                        break
          
        # Catch if the initialization memory assignment was not sufficient for the max. # of fibers available to any target.
        # May well have thrown an error by this point if that's the case. 
        assert  np.all(final['NGOOD_FIBERS'] < MAX_NFIBER)
        assert  np.all(final['NGOOD_TILES']  < MAX_NTILE)
            
        if printit:
            toprint = final[final['NGOOD_FIBERS'] >= 1]

            print(toprint['NGOOD_FIBERS', 'NGOOD_TILES', 'GOOD_FIBERS', 'GOOD_TILES', 'FIBPRIORITY'])
      
    return  final

def set_ang_veto(final):
    # Set for all targets. 
    final['ANG_VETO_FLAG']            = 0

    # Target can reached by a fiber in the survey. 
    reachable = final['NGOOD_FIBERS'] > 0
    final['ANG_VETO_FLAG'][reachable] = 1

    # All reachable fibers for this target were blocked by a higher (initial) priority class.
    blocked    = []
    
    for i, x in enumerate(final['PRIORITY_INIT']):
        blocks = x < final['FIBPRIORITY'][i,:]
        blocked.append(np.all(blocks))
        
    final['ANG_VETO_FLAG'][blocked] = 3
    
    return  final

def _in_zlimits(zee):
    ##  e.g. ELG-specific?  Basic in the meantime. 
    return  (zee > 0.05) & (zee < 2.1)

def set_z_veto(final):
    ##  No secure classification:  DCHISQ. DIFF RAISED ON DEGENERATE TYPES, e.g. QSO & STAR.  
    ##  Requires... zbest.  E.g. 
    ##  /project/projectdirs/desi/datachallenge/svdc-summer2019/svdc2019c/spectro/redux/v1/spectra-64/101/10151/zbest-64-10151.fits
    ##  
    ##  Easiest via redrock.  Remap ZWARN in the meantime.   
    badz                                        = final['ZWARN'] > 0
    final['Z_VETO_FLAG'][badz]                  = 4

    ##  i.e. non-stellar                                                                                                                                                                                                                                                                                                                                                                                          
    cosmological                                = np.array([x in ['GALAXY', 'QSO'] for x in final['SPECTYPE']])
    final['Z_VETO_FLAG'][~badz & ~cosmological] = 3
    
    in_zlimits                                  = _in_zlimits(final['Z'])
    final['Z_VETO_FLAG'][~badz &  in_zlimits]   = 1
    final['Z_VETO_FLAG'][~badz & ~in_zlimits]   = 2
        
    return  final

def in_cosmo(final):
    no_zveto =   final['Z_VETO_FLAG'] < 2
    no_aveto = final['ANG_VETO_FLAG'] < 2

    no_veto  = no_zveto & no_aveto

    final['IN_COSMO'] = no_veto.astype(np.int)
    
    return  final

def set_assigncomplete(dst_tiles, final):
    return  dst_tiles


def lsscat_gen(root, prod, odir):
    '''
    # Get required parent pipeline files. 
    _mtl            = root + 'targets/mtl.fits'

    # Needed for fiberassign runs. 
    # _skies        = '/project/projectdirs/desi/target/catalogs/dr8/0.31.0/skies/skies-dr8-0.31.0.fits'

    # targets       = Table(fits.open(root + 'targets/targets.fits')[1].data)
    
    mtl             = Table(fits.open(root + 'targets/mtl.fits')[1].data)
    # exps          = Table(fits.open(root + 'survey/sv_exposures.fits')[1].data)

    # HACK:  Cumulative SNR2 achieved on a given tile by this exposure. Weird that this was missing?                                                                                                               
    # exps['SNR2']  = np.random.uniform(0.9, 1.1, len(exps))
    
    _tiles          = fits.open(root + 'survey/SV-tiles.fits')[1]
    tiles           = Table(fits.open(root + 'survey/SV-tiles.fits')[1].data)

    # zcat          = Table(fits.open(root + 'spectro/redux/{}/zcatalog-{}.fits'.format(prod, prod))[1].data)
    
    # DR8 randoms. 
    rows          = np.arange(len(mtl))
    randoms       = fitsio.read('/project/projectdirs/desi/target/catalogs/dr8/0.31.0/randoms/randoms-inside-dr8-0.31.0-4.fits', rows=rows)

    # A list of targets that could have been assigned to each fiber;
    # desidatamodel.readthedocs.io/en/latest/DESI_TARGET/fiberassign/tile-TILEID-FIELDNUM.html#hdu2

    rand_assignable, rand_assigned = get_randassign(root + 'fiberassign/tile-072015.fits', mtl=randoms, gen=True, randoms=randoms)
    '''
    rand_assignable, rand_assigned = get_randassign(root + 'fiberassign/tile-072015.fits', mtl=None, gen=True, randoms=None, lite=True)
    
    exit(0)
    
    # TARGETID | FIBERID | LOCATION
    # print(rand_assignable)
    
    # Date stamped tiles meeting quality requirements.  
    MJD0            = 58853.51129915193
    MJD1            = 58881.40164263005
    
    dst_tiles, hdr           = datestamped_tiles(exps, MJD0, MJD1, printit=False)
    
    # Write here (should add header).  Assign tile assignment completeness below. 
    dst_tiles.write(odir + 'dst_tiles.fits', format='fits', overwrite=True)
    
    '''
    # Sort Targets by TARGETID
    zcat.sort('TARGETID')
    targets.sort('TARGETID')

    # Unique TARGETID lists.
    assert  np.all(np.unique(zcat['TARGETID'])    == zcat['TARGETID'])
    assert  np.all(np.unique(targets['TARGETID']) == targets['TARGETID'])
    '''
    # Catalogue generation.
    ntarget                        = len(targets)
    
    final                          = lss_catalog(ntarget)

    final['TARGETID']              = targets['TARGETID']
    
    final['RA']                    = targets['RA']
    final['DEC']                   = targets['DEC']

    final['IN_IMAGING']            = 1.0

    ##  Neglects area lost to e.g. GFA.
    final['IN_DESI']               = is_point_in_desi(tiles, final['RA'], final['DEC']).astype(np.int)

    final['DESI_TARGET']           = targets['SV1_DESI_TARGET']
    final['BGS_TARGET']            = targets['SV1_BGS_TARGET']

    final['PRIORITY_INIT']         = targets['PRIORITY_INIT']
    final['SUBPRIORITY']           = targets['SUBPRIORITY']

    in_zcat                        = np.in1d( zcat['TARGETID'], final['TARGETID'])
    in_final                       = np.in1d(final['TARGETID'],  zcat['TARGETID'])
    
    final['Z'][in_final]           = zcat['Z'][in_zcat]
    final['ZERR'][in_final]        = zcat['ZERR'][in_zcat]
    final['ZWARN'][in_final]       = zcat['ZWARN'][in_zcat]
    final['SPECTYPE'][in_final]    = zcat['SPECTYPE'][in_zcat]

    ##  Tag all galaxies that appeared in the redshift cats., even if bad etc. 
    final['Z_VETO_FLAG']           = 0

    final                          = pop_assigned(final, printit=True)

    final                          = set_ang_veto(final)          
    
    final                          = set_z_veto(final)
    
    final                          = blind_distances(final)
        
    final                          = in_cosmo(final)
    
    final                          = final[final['Z'] > 0.]
    
    final.write(odir + 'lsscat_v1.fits', format='fits', overwrite=True)

    final.pprint(max_width=-1)
    
    # Add tile assignment completeness based on reachable fibers from final. 
    dst_tiles                      = set_assigncomplete(dst_tiles, final)
    
    
if __name__ == '__main__':
    odir            = '/project/projectdirs/desi/www/users/mjwilson/'

    prod            = 'v1'
    root            = '/project/projectdirs/desi/datachallenge/svdc-summer2019/svdc2019c/'    

    lsscat_gen(root, prod, odir)

    print('\n\nDone.\n\n')
