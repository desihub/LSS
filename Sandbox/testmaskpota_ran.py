'''
This is a little demo script for the Assignment.check_avail_collisions() function.
'''

import numpy as np

from astropy.table import Table,join,setdiff

from fiberassign.hardware import load_hardware
from fiberassign.tiles import load_tiles
from fiberassign.targets import Targets, TargetsAvailable, LocationsAvailable, create_tagalong, load_target_file, targets_in_tiles
from fiberassign.assign import Assignment

from fiberassign.utils import Logger

import fitsio

t = Table.read('/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/tiles-DARK.fits')
#print('tiles:', t)

margins = dict(pos=0.05,
                   petal=0.4,
                   gfa=0.4)


#def main():

    # from LSS.mkCat_singletile.fa4lsscat import getfatiles
    # getfatiles()
    # return
log = Logger.get()
rann=0
n = 0
for tile in t['TILEID']:    


    #tile = 1230
    ts = '%06i' % tile

    fbah = fitsio.read_header('/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'+ts[:3]+'/fiberassign-'+ts+'.fits.gz')
    dt = fbah['RUNDATE']#[:19]
    hw = load_hardware(rundate=dt, add_margins=margins)
    pr = 'DARK'
    t['OBSCONDITIONS'] = 516
    t['IN_DESI'] = 1
    t['MTLTIME'] = fbah['MTLTIME']
    t['FA_RUN'] = fbah['FA_RUN']
    t['PROGRAM'] = pr
    #t['FA_HA'] = fbah['FA_HA']
    #t['FIELDROT'] = fbah['FIELDROT']
    obsha = fbah['FA_HA']
    obstheta = fbah['FIELDROT']

    t.write('tiles.fits', overwrite=True)

    tiles = load_tiles(
        tiles_file='tiles.fits',
        select=[tile])

    tids = tiles.id
    print('Tile ids:', tids)
    I = np.flatnonzero(np.array(tids) == tile)
    assert(len(I) == 1)
    i = I[0]
    tile_ra  = tiles.ra[i]
    tile_dec = tiles.dec[i]

    # Create empty target list
    tgs = Targets()
    # Create structure for carrying along auxiliary target data not needed by C++.
    plate_radec=True
    tagalong = create_tagalong(plate_radec=plate_radec)

    # Load target files...
    load_target_file(tgs, tagalong, '/global/cfs/cdirs/desi/survey/catalogs/main/LSS/random0/tilenofa-%i.fits' % tile)
    #loading it again straight to table format because I can't quickly figure out exactly where targetid,ra,dec gets stored
    tar_tab = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/main/LSS/random'+str(rann)+'/tilenofa-%i.fits' % tile,columns =['TARGETID','RA','DEC'])

    # Find targets within tiles, and project their RA,Dec positions
    # into focal-plane coordinates.
    tile_targetids, tile_x, tile_y, tile_xy_cs5 = targets_in_tiles(hw, tgs, tiles, tagalong)
    # Compute the targets available to each fiber for each tile.
    tgsavail = TargetsAvailable(hw, tiles, tile_targetids, tile_x, tile_y)
    # Compute the fibers on all tiles available for each target and sky
    favail = LocationsAvailable(tgsavail)
    
    # FAKE stucksky
    stucksky = {}

    # Create assignment object
    asgn = Assignment(tgs, tgsavail, favail, stucksky)
    tgsavail = asgn.targets_avail()
    avail = tgsavail.tile_data(tile)
    navail = np.sum([len(avail[x]) for x in avail.keys()])
    fibers = dict(hw.loc_fiber)
    fdata = Table()
    fdata['LOCATION'] = np.zeros(navail,dtype=int)
    fdata['FIBER'] = np.zeros(navail,dtype=int)
    fdata['TARGETID'] = np.zeros(navail,dtype=int)
    
    off = 0
    # The "FAVAIL" (available targets) HDU is sorted first by LOCATION,
    # then by TARGETID.
    for lid in sorted(avail.keys()):
        # lid (location id) is a scalar, tg (target ids) is an array
        tg = avail[lid]
        fdata['LOCATION'][off:off+len(tg)] = lid
        fdata['FIBER']   [off:off+len(tg)] = fibers[lid]
        fdata['TARGETID'][off:off+len(tg)] = sorted(tg)
        off += len(tg)
    #print(avail.keys())
    print(len(fdata))
    fdata = join(fdata,tar_tab,keys=['TARGETID'],join_type='left')
    print(len(fdata))
    coll = asgn.check_avail_collisions(tile)
    kl = np.array(list(coll.keys())).transpose()
    locs = kl[0]
    ids = kl[1]
    locids = ids*10000+locs
    #print('collisions:', coll)
    print('N collisions:', len(coll))
    # coll: dict (loc, targetid) -> bitmask
    forig = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/main/LSS/random0/fba-'+ts+'.fits',ext='FAVAIL')
    #print(coll)
    locidsin = np.isin(fdata['LOCATION']+10000*fdata['TARGETID'],locids)
    locidsino = np.isin(forig['LOCATION']+10000*forig['TARGETID'],locids)
    print(np.sum(locidsin),np.sum(locidsino),len(fdata))
    #jt = setdiff(fdata,Table(forig),keys=['TARGETID','FIBER','LOCATION'])#,join_type='inner')
    #jto = setdiff(Table(forig),fdata,keys=['TARGETID','FIBER','LOCATION'])
    #print(len(jt),len(jto),len(forig),len(fdata))
    n += 1
    if n >= 1:
        break

