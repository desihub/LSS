'''
find and record potential assignments that would never have actually been assigned, due to collisions
'''
import os
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
parser.add_argument("--outmode", choices=['test','prod'],default='test')

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

    coll = asgn.check_avail_collisions(tile)
    kl = np.array(list(coll.keys())).transpose()
    locs = kl[0]
    ids = kl[1]
    locids = ids*10000+locs
    print('N collisions:', len(coll))
    # coll: dict (loc, targetid) -> bitmask
    forig = fitsio.read('/dvs_ro/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'+ts[:3]+'/fiberassign-'+ts+'.fits.gz',ext='POTENTIAL_ASSIGNMENTS')
    #print(coll)
    selo = np.isin(forig['TARGETID'],ttids)
    forig = forig[selo]
    locidsin = np.isin(forig['LOCATION']+10000*forig['TARGETID'],locids)
    colltab = Table(forig[locidsin])
    colltab['TILEID'] = tile
    return colltab
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
    colltot = np.concatenate(res)
    if args.outmode == 'prod':
        outroot = '/global/cfs/cdirs/desi/survey/catalogs/'
    if args.outmode == 'test':
        outroot = os.getenv('SCRATCH')+'/'
    if not os.path.exists(outroot+args.survey+'/LSS'):
        os.makedirs(outroot+args.survey+'/LSS')

    common.write_LSS(colltot,outroot+args.survey+'/LSS/collisions-'+args.prog+'.fits')

