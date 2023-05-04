'''
This is a little demo script for the Assignment.check_avail_collisions() function.
'''

import numpy as np

from astropy.table import Table

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

n = 0
badl = []
badtot = 0
tls = [1100,1112,1705,2125,2143,2187, 2414, 2559, 2564, 2836, 3160, 4226, 4582, 5774, 5981, 6153, 6289, \
6293, 7169, 7206, 7216, 7478, 7814, 7816, 8371, 8772, 8956, 9406, 9828, 9856, 9931, 10384, \
11239, 11505, 11709, 11709, 11943]
for tile in t['TILEID']:    
#for tile in tls:

    #tile = 1230
    ts = '%06i' % tile

    fbah = fitsio.read_header('/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'+ts[:3]+'/fiberassign-'+ts+'.fits.gz')
    #if fbah['FA_VER'][0] == '5':
    dt = fbah['RUNDATE']#[:19]
    hw = load_hardware(rundate=dt, add_margins=margins)
    t = Table(fitsio.read('/global/cfs/cdirs/desi/survey/fiberassign/main/'+ts[:3]+'/'+ts+'-tiles.fits'))
    t['MTLTIME'] = fbah['MTLTIME']
    t['FA_RUN'] = fbah['FA_RUN']
    obsha = fbah['FA_HA']
    obstheta = fbah['FIELDROT']


    t.write('tiles.fits', overwrite=True)

    tiles = load_tiles(
        tiles_file='tiles.fits',obsha=obsha,obstheta=obstheta,
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
    load_target_file(tgs, tagalong, '/global/cfs/cdirs/desi/survey/fiberassign/main/'+ts[:3]+'/'+ts+'-targ.fits',rundate=dt)
    load_target_file(tgs, tagalong, '/global/cfs/cdirs/desi/survey/fiberassign/main/'+ts[:3]+'/'+ts+'-scnd.fits',rundate=dt)
    load_target_file(tgs, tagalong, '/global/cfs/cdirs/desi/survey/fiberassign/main/'+ts[:3]+'/'+ts+'-sky.fits',rundate=dt)

    ttids = fitsio.read('/global/cfs/cdirs/desi/survey/fiberassign/main/'+ts[:3]+'/'+ts+'-targ.fits')['TARGETID']


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
    #locs = []
    #ids = []
    #for key in coll.keys():
    #    locs.append(key[0])
    #    ids.append(key[1])
    #locs = np.array(locs)
    #ids = np.array(ids)
    #for loc,id in zip(locs,ids):
    #    try:
    #        bit = coll[(loc,id)]
    #    except:
    #        print(loc,id,'key error')
    #        break    
    #sel = np.isin(ids,ttids)
    #locs = locs[sel]
    #ids = ids[sel]
    #print('collisions:', coll)
    print('N collisions:', len(coll))
    # coll: dict (loc, targetid) -> bitmask
    forig = fitsio.read('/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'+ts[:3]+'/fiberassign-'+ts+'.fits.gz',ext='POTENTIAL_ASSIGNMENTS')
    #print(coll)
    selo = np.isin(forig['TARGETID'],ttids)
    forig = forig[selo]
    locidsin = np.isin(forig['LOCATION']+10000*forig['TARGETID'],locids)
    #locsin = np.isin(forig['LOCATION'],locs)
    #idsin = np.isin(forig['TARGETID'],ids)
    masked = locidsin#locsin&idsin
    #print(np.sum(locsin),np.sum(idsin),np.sum(masked),len(forig))
    print(np.sum(masked),len(forig))
    print('checking actual assignments, should not find any are masked')
    forig = fitsio.read('/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'+ts[:3]+'/fiberassign-'+ts+'.fits.gz',ext='FIBERASSIGN')
    #print(coll)
    #locsin = np.isin(forig['LOCATION'],locs)
    #idsin = np.isin(forig['TARGETID'],ids)
    locidsin = np.isin(forig['LOCATION']+10000*forig['TARGETID'],locids)
    masked = locidsin#locsin&idsin
    print(np.sum(masked))
    if np.sum(masked) != 0:
        print('BAD, assigned id/location is in mask')
        #print(forig[masked])
        print(forig['LOCATION'][masked],forig['TARGETID'][masked])
        for i in range(0,len(forig[masked])):
            loc = forig[masked][i]['LOCATION']
            id = forig[masked][i]['TARGETID']
            print(loc,id,coll[(loc,id)])
            badl.append((tile,loc,id,coll[(loc,id)]))
        badtot += np.sum(masked)
        

    n += 1
    #else:
    #    print(ts,fbah['FA_VER'])
    print(n,badtot)
    #if n >= 100:
    #    break
print('#the tileid, location, targetid, bits that were bad are ')
print(badl)
