'''
This is a little demo script for the Assignment.check_avail_collisions() function.
'''

import numpy as np
import argparse
from astropy.table import Table

from fiberassign.hardware import load_hardware
from fiberassign.tiles import load_tiles
from fiberassign.targets import Targets, TargetsAvailable, LocationsAvailable, create_tagalong, load_target_file, targets_in_tiles
from fiberassign.assign import Assignment

from fiberassign.utils import Logger

import fitsio

import LSS.common_tools as common


parser = argparse.ArgumentParser()
parser.add_argument("--prog", choices=['DARK','BRIGHT'])
parser.add_argument("--survey", choices=['Y1','DA2'],default='DA2')

args = parser.parse_args()


t = Table.read('/global/cfs/cdirs/desi/survey/catalogs/'+args.survey+'/LSS/tiles-'+args.prog+'.fits')
#print('tiles:', t)

margins = dict(pos=0.05,
                   petal=0.4,
                   gfa=0.4)


log = Logger.get()

n = 0
colls = []

#for tile in t['TILEID']:    
#for tile in tls:

    #tile = 1230
def getcoll(tile):
    ts = '%06i' % tile

    fbah = fitsio.read_header('/dvs_ro/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'+ts[:3]+'/fiberassign-'+ts+'.fits.gz')
    #if fbah['FA_VER'][0] == '5':
    dt = fbah['RUNDATE']#[:19]
    hw = load_hardware(rundate=dt, add_margins=margins)
    obsha = fbah['FA_HA']
    obstheta = fbah['FIELDROT']

    tiles = load_tiles(
        tiles_file='/dvs_ro/cfs/cdirs/desi/survey/fiberassign/main/'+ts[:3]+'/'+ts+'-tiles.fits',obsha=obsha,obstheta=obstheta,
        select=[tile])

    tids = tiles.id
    #print('Tile ids:', tids)
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
    load_target_file(tgs, tagalong, '/dvs_ro/cfs/cdirs/desi/survey/fiberassign/main/'+ts[:3]+'/'+ts+'-targ.fits',rundate=dt)
    load_target_file(tgs, tagalong, '/dvs_ro/cfs/cdirs/desi/survey/fiberassign/main/'+ts[:3]+'/'+ts+'-scnd.fits',rundate=dt)
    load_target_file(tgs, tagalong, '/dvs_ro/cfs/cdirs/desi/survey/fiberassign/main/'+ts[:3]+'/'+ts+'-sky.fits',rundate=dt)

    ttids = fitsio.read('/dvs_ro/cfs/cdirs/desi/survey/fiberassign/main/'+ts[:3]+'/'+ts+'-targ.fits')['TARGETID']


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

    coll = asgn.check_avail_collisions(tile)
    kl = np.array(list(coll.keys())).transpose()
    locs = kl[0]
    ids = kl[1]
    locids = ids*10000+locs
    #print('N collisions:', len(coll))
    # coll: dict (loc, targetid) -> bitmask
    forig = fitsio.read('/dvs_ro/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'+ts[:3]+'/fiberassign-'+ts+'.fits.gz',ext='POTENTIAL_ASSIGNMENTS')
    #print(coll)
    selo = np.isin(forig['TARGETID'],ttids)
    forig = forig[selo]
    check = np.isin(forig['TARGETID'],fdata['TARGETID'])
    return (tile,len(check),np.sum(check),dt)
    #locidsin = np.isin(forig['LOCATION']+10000*forig['TARGETID'],locids)
    #colltab = Table(forig[locidsin])
    #colltab['TILEID'] = tile
    #return colltab
    #colls.append(colltab)
        

    #n += 1
    #else:
    #    print(ts,fbah['FA_VER'])
    #print(n,len(t))
    #if n >= 100:
    #    break

#colltot = np.concatenate(colls)

if __name__ == '__main__':
    from multiprocessing import Pool
    tls = list(t['TILEID'])#[:10])
    with Pool(processes=128) as pool:
        res = pool.map(getcoll, tls)
    tiles_notm = []
    nfail = 0
    for i in range(0,len(res)):
        if res[i][2] != res[i][1]:
            print(res[i])
            tiles_notm.append(res[i][0])
            nfail += 1
    print('the number of mismatches is '+str(nfail))
    
    #colltot = np.concatenate(res)
    #common.write_LSS(colltot,'/global/cfs/cdirs/desi/survey/catalogs/'+args.survey+'/LSS/collisions-'+args.prog+'.fits')

