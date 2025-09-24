'''
Find all potential assignments (and collisions) for DA2 randoms
'''

import numpy as np
import os
from astropy.table import Table, join
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
parser.add_argument("--survey", choices=['Y1','DA2'],default='DA2')
parser.add_argument("--minr", default=0,type=int)
parser.add_argument("--maxr", default=18,type=int)
parser.add_argument("--outmode", choices=['test','prod'],default='test')

args = parser.parse_args()


tiletab = Table.read('/dvs_ro/cfs/cdirs/desi/survey/catalogs/'+args.survey+'/LSS/tiles-'+args.prog+'.fits')

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

    fbah = fitsio.read_header('/dvs_ro/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'+ts[:3]+'/fiberassign-'+ts+'.fits.gz')
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
    load_target_file(tgs, tagalong, '/dvs_ro/cfs/cdirs/desi/survey/catalogs/main/LSS/random'+str(rann)+'/tilenofa-%i.fits' % tile)
    #loading it again straight to table format because I can't quickly figure out exactly where targetid,ra,dec gets stored
    tar_tab = fitsio.read('/dvs_ro/cfs/cdirs/desi/survey/catalogs/main/LSS/random'+str(rann)+'/tilenofa-%i.fits' % tile,columns =['TARGETID','RA','DEC'])

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
    fdata = join(fdata,tar_tab,keys=['TARGETID'],join_type='left')
    coll = asgn.check_avail_collisions(tile)
    kl = np.array(list(coll.keys())).transpose()
    locs = kl[0]
    ids = kl[1]
    locids = ids*10000+locs
    #print('collisions:', coll)
    print('N collisions:', len(coll))
    # coll: dict (loc, targetid) -> bitmask
    #forig = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/main/LSS/random'+str(rann)+'/fba-'+ts+'.fits',ext='FAVAIL')
    #print(coll)
    locidsin = np.isin(fdata['LOCATION']+10000*fdata['TARGETID'],locids)
    print('N collisions original:',np.sum(locidsin),len(fdata))
    #colltab = Table(forig[locidsin])
    fdata['TILEID'] = tile
    fdata['COLLISION'] = locidsin
    return fdata
    
if __name__ == '__main__':
    from multiprocessing import Pool
    tls = list(tiletab['TILEID'])#[:10])
    inds = np.arange(len(tls))
    if args.outmode == 'prod':
        outroot = '/global/cfs/cdirs/desi/survey/catalogs/'
    if args.outmode == 'test':
        outroot = os.getenv('SCRATCH')+'/'
    if not os.path.exists(outroot+args.survey+'/LSS'):
        os.makedirs(outroot+args.survey+'/LSS')

    for rann in range(int(args.minr),int(args.maxr)):
        with Pool(processes=128) as pool:
            res = pool.map(getcoll, inds)
        colltot = np.concatenate(res)
        print(len(colltot),np.sum(colltot['COLLISION']))
        outdir = outroot+args.survey+'/LSS/random'+str(rann)
        if not os.path.exists(outdir):        
            os.makedirs(outdir)
        common.write_LSS_scratchcp(colltot,outdir+'/pota-'+args.prog+'.fits')

