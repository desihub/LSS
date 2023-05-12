'''
Find all of the collisions for Y1 randoms
'''

import numpy as np
import os
from astropy.table import Table
import argparse
from fiberassign.hardware import load_hardware
from fiberassign.tiles import load_tiles
from fiberassign.targets import Targets, TargetsAvailable, LocationsAvailable, create_tagalong, load_target_file, targets_in_tiles
from fiberassign.assign import Assignment

from fiberassign.utils import Logger

import fitsio

import LSS.common_tools as common

parser = argparse.ArgumentParser()
parser.add_argument("--prog", choices=['DARK','BRIGHT'])

args = parser.parse_args()


tiletab = Table.read('/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/tiles-'+args.prog+'.fits')

margins = dict(pos=0.05,
                   petal=0.4,
                   gfa=0.4)


#def main():

    # from LSS.mkCat_singletile.fa4lsscat import getfatiles
    # getfatiles()
    # return
log = Logger.get()
rann = 0
n = 0
#for tile in t['TILEID']:    
def getcoll(ind):

    #tile = 1230
    tile = tiletab[ind]['TILEID']
    ts = '%06i' % tile

    fbah = fitsio.read_header('/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'+ts[:3]+'/fiberassign-'+ts+'.fits.gz')
    dt = fbah['RUNDATE']#[:19]
    pr = args.prog
    t = Table(tiletab[ind])
    t['OBSCONDITIONS'] = 516
    t['IN_DESI'] = 1
    t['MTLTIME'] = fbah['MTLTIME']
    t['FA_RUN'] = fbah['FA_RUN']
    t['PROGRAM'] = pr
    obsha = fbah['FA_HA']
    obstheta = fbah['FIELDROT']

    hw = load_hardware(rundate=dt, add_margins=margins)

    t.write(os.environ['SCRATCH']+'/rantiles/'+str(tile)+'-'+str(rann)+'-tiles.fits', overwrite=True)

    tiles = load_tiles(
        tiles_file=os.environ['SCRATCH']+'/rantiles/'+str(tile)+'-'+str(rann)+'-tiles.fits',obsha=obsha,obstheta=obstheta,
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
    load_target_file(tgs, tagalong, '/global/cfs/cdirs/desi/survey/catalogs/main/LSS/random'+str(rann)+'/tilenofa-%i.fits' % tile)

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
    #print('collisions:', coll)
    print('N collisions:', len(coll))
    # coll: dict (loc, targetid) -> bitmask
    forig = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/main/LSS/random'+str(rann)+'/fba-'+ts+'.fits',ext='FAVAIL')
    #print(coll)
    locidsin = np.isin(forig['LOCATION']+10000*forig['TARGETID'],locids)
    print('N collisions original:',np.sum(locidsin))
    colltab = Table(forig[locidsin])
    colltab['TILEID'] = tile
    return colltab
    
if __name__ == '__main__':
    from multiprocessing import Pool
    tls = list(tiletab['TILEID'])#[:10])
    inds = np.arange(len(tls))
    for rann in range(0,18):
        with Pool(processes=128) as pool:
            res = pool.map(getcoll, inds)
        colltot = np.concatenate(res)
        common.write_LSS(colltot,'/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/random'+str(rann)+'collisions-'+args.prog+'.fits')

