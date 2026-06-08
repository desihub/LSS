'''
Find all potential assignment and counts tiles for Y1 mocks
Use the following environment
source /global/common/software/desi/desi_environment.sh main
'''

import numpy as np
import os
from astropy.table import Table, join, vstack
import argparse
from fiberassign.hardware import load_hardware, get_default_exclusion_margins
from fiberassign._internal import Hardware
from fiberassign.tiles import load_tiles
from fiberassign.targets import Targets, TargetsAvailable, LocationsAvailable, create_tagalong, load_target_file, targets_in_tiles
from fiberassign.assign import Assignment

from fiberassign.utils import Logger

from desitarget.io import read_targets_in_tiles
import desimodel.focalplane
import desimodel.footprint
trad = desimodel.focalplane.get_tile_radius_deg()*1.1 #make 10% greater just in case


import fitsio

import LSS.common_tools as common
#from LSS.imaging import get_nobsandmask
from LSS.main.cattools import count_tiles_better
from LSS.globals import main

import bisect
import time
from datetime import datetime
import multiprocessing

t_start = time.time()

log = Logger.get()

parser = argparse.ArgumentParser()
parser.add_argument("--prog", choices=['DARK','BRIGHT'],default='DARK')
parser.add_argument("--mock", default='ab2ndgen')
parser.add_argument("--mock_version",default='')
parser.add_argument("--realization")
parser.add_argument("--getcoll",default='y')
parser.add_argument("--base_output", help="base directory for output",default='/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/')
parser.add_argument("--tracer", help="tracer for CutSky EZ mocks", default=None)
parser.add_argument("--base_input", help="base directory for input for EZ mocks 6Gpc", default = None)
parser.add_argument("--tile-temp-dir", help="Directory for temp tile files, default %(default)s",
                    default=os.path.join(os.environ['SCRATCH'], 'rantiles'))
parser.add_argument("--counttiles", default = 'n')
parser.add_argument("--secgen_ver", default = None)
parser.add_argument("--nprocs", help="Number of multiprocessing processes to use, default %(default)i",
                    default=multiprocessing.cpu_count()//2, type=int)

# On Perlmutter, this read-only access point can be *much* faster thanks to aggressive caching.
#   If you didn't want this for some reason, you could revert '/dvs_ro/cfs/cdirs/desi' to '/global/cfs/cdirs/desi' in the following.
desi_input_dir = os.getenv('DESI_ROOT_READONLY', default='/dvs_ro/cfs/cdirs/desi')

args = parser.parse_args()
print(args)

if args.mock == 'ab2ndgen':
    #infn = args.base_output+'FirstGenMocks/AbacusSummit/forFA'+args.realization+'_matched_input_full_masknobs.fits'
    #infn = args.base_output+'SecondGenMocks/AbacusSummit/forFA'+args.realization+'.fits'
    infn = os.path.join(args.base_output+'SecondGenMocks', 'AbacusSummit'+args.mock_version, 'forFA'+args.realization+'.fits')
    log.info('Reading %s' % infn)
    tars = fitsio.read(infn)
    tarcols = list(tars.dtype.names)
    #tileoutdir = args.base_output+'SecondGenMocks/AbacusSummit/tartiles'+args.realization+'/'
    tileoutdir = os.path.join(os.getenv('SCRATCH'), 'SecondGenMocks', 'AbacusSummit', 'tartiles'+args.realization)
    paoutdir = os.path.join(args.base_output+'SecondGenMocks', 'AbacusSummit'+args.mock_version, 'mock'+args.realization)
elif args.mock == 'ezmocks6':
    # #tr = args.tracer
    # rz = args.realization
    # print("Doing %s"%tr)
    #
    # if  tr == "LRG":
    #     infn1 = "/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/EZmock/CutSky_6Gpc/LRG/z0.800/cutsky_LRG_z0.800_EZmock_B6000G1536Z0.8N216424548_b0.385d4r169c0.3_seed%s_NGC.fits"%rz
    #     infn2 = "/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/EZmock/CutSky_6Gpc/LRG/z0.800/cutsky_LRG_z0.800_EZmock_B6000G1536Z0.8N216424548_b0.385d4r169c0.3_seed%s_SGC.fits"%rz
    # elif tr == "ELG":
    #     infn1 = "/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/EZmock/CutSky_6Gpc/ELG/z1.100/cutsky_ELG_z1.100_EZmock_B6000G1536Z1.1N648012690_b0.345d1.45r40c0.05_seed%s_NGC.fits"%rz
    #     infn2 = "/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/EZmock/CutSky_6Gpc/ELG/z1.100/cutsky_ELG_z1.100_EZmock_B6000G1536Z1.1N648012690_b0.345d1.45r40c0.05_seed%s_SGC.fits"%rz
    # elif tr == "QSO":
    #     infn1 = "/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/EZmock/CutSky_6Gpc/QSO/z1.400/cutsky_QSO_z1.400_EZmock_B6000G1536Z1.4N27395172_b0.053d1.13r0c0.6_seed%s_NGC.fits"%rz
    #     infn2 = "/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/EZmock/CutSky_6Gpc/QSO/z1.400/cutsky_QSO_z1.400_EZmock_B6000G1536Z1.4N27395172_b0.053d1.13r0c0.6_seed%s_SGC.fits"%rz
    ## infn1 = "/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/EZmock/CutSky_6Gpc/LRG/z0.800/cutsky_LRG_z0.800_EZmock_B6000G1536Z0.8N216424548_b0.385d4r169c0.3_seed1_NGC.fits"
    ## infn2 = "/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/EZmock/CutSky_6Gpc/LRG/z0.800/cutsky_LRG_z0.800_EZmock_B6000G1536Z0.8N216424548_b0.385d4r169c0.3_seed1_SGC.fits"
    # tars1 = Table.read(infn1)#fitsio.read(infn1)
    # tars2 = Table.read(infn2)#fitsio.read(infn2)
    # tars1["GALCAP"] = "N"
    # tars2["GALCAP"] = "S"
    # tars = vstack([tars1, tars2])
    # tars['TARGETID'] = np.arange(len(tars))

    infn = os.path.join(args.base_input + 'EZMocks_6Gpc', 'EZMocks_6Gpc_' + args.realization + '.fits')
    tars = fitsio.read(infn)
    tarcols = list(tars.dtype.names)#['TARGETID','RA','DEC', 'Z','Z_COSMO','GALCAP', 'NZ', 'RAW_NZ']

    tileoutdir = os.path.join(args.base_output+'EZMocks_6Gpc', 'tartiles'+args.realization)
    paoutdir = os.path.join(args.base_output+'EZMocks_6Gpc', 'EzMocks', 'mock'+args.realization)
    if args.tracer is not None:
        tileoutdir = os.path.join(tileoutdir, args.tracer)
        paoutdir = os.path.join(paoutdir, args.tracer)

# Ensure that the targets file is sorted by Dec.
t0 = time.time()
is_sorted = np.all(tars['DEC'][:-1] <= tars['DEC'][1:])
if not is_sorted:
    I = np.argsort(tars['DEC'])
    tars = tars[I]
t1 = time.time()
log.info('Sorting/verifying mocks: %.1f' % (t1-t0))

if not os.path.exists(tileoutdir):
    os.makedirs(tileoutdir)
    #print('made '+tileoutdir)
if not os.path.exists(paoutdir):
    os.makedirs(paoutdir)
    #print('made '+paoutdir)

tiletab = Table.read(os.path.join(desi_input_dir, 'survey', 'catalogs', 'Y1', 'LSS', 'tiles-'+args.prog+'.fits'))
log.info('Reading startup globals: %.3f' % (time.time() - t_start))

def get_tile_targ(tile):
    '''
    Creates an astropy Table of (mock) targets within the given `tile`.
    '''
    tdec = tile['DEC']
    decmin = tdec - trad
    decmax = tdec + trad
    dec = tars['DEC']
    # The `tars` global table of targets is sorted by Dec.  We therefore only need to look at
    # indices that can possibly be within range given just the Dec distance (decmin to decmax).
    # "bisect_left" is way faster than "np.searchsorted"!
    #i0,i1 = np.searchsorted(dec, [np.float32(decmin), np.float32(decmax)])
    i0 = bisect.bisect_left(dec, decmin)
    i1 = bisect.bisect_left(dec, decmax, lo=i0)
    Idec = slice(i0, i1+1)
    inds = desimodel.footprint.find_points_radec(tile['RA'], tdec,
                                                 tars['RA'][Idec], tars['DEC'][Idec])
    rtw = tars[i0 + np.array(inds)]
    rmtl = Table(rtw)
    #print('made table')
    del rtw
    if 'DESI_TARGET' not in tarcols:
        rmtl['DESI_TARGET'] = np.ones(len(rmtl),dtype=int)*2
    if 'NUMOBS_INIT' not in tarcols:
        rmtl['NUMOBS_INIT'] = np.zeros(len(rmtl),dtype=int)
    if 'NUMOBS_MORE' not in tarcols:
        rmtl['NUMOBS_MORE'] = np.ones(len(rmtl),dtype=int)
    if 'PRIORITY' not in tarcols:
        rmtl['PRIORITY'] = np.ones(len(rmtl),dtype=int)*3400
    #if 'OBSCONDITIONS' not in tarcols:
    rmtl['OBSCONDITIONS'] = np.ones(len(rmtl),dtype=int)*516#forcing it to match value assumed below
    if 'SUBPRIORITY' not in tarcols:
        rmtl['SUBPRIORITY'] = np.random.random(len(rmtl))
    return rmtl

def write_tile_targ(ind):
    '''
    Write the targets file for the single tile table index "ind".
    '''
    tile = tiletab[ind]
    fname = os.path.join(tileoutdir, 'tilenofa-'+str(tile['TILEID'])+'.fits')
    log.info('creating %s' % fname)
    rmtl = get_tile_targ(tile)
    #print('added columns for', fname)
    rmtl.write(fname, format='fits', overwrite=True)
    #print('added columns, wrote to', fname)

margins = get_default_exclusion_margins()
rann = 0
n = 0

def getpa(ind):
    #tile = 1230
    tile = tiletab[ind]['TILEID']
    ts = '%06i' % tile

    fbah = fitsio.read_header(os.path.join(desi_input_dir, 'target', 'fiberassign', 'tiles', 'trunk', ts[:3], 'fiberassign-'+ts+'.fits.gz'))
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

    tt = parse_datetime(dt)
    hw = get_hardware_for_time(tt)
    assert(hw is not None)

    tilefn = os.path.join(args.tile_temp_dir, str(tile)+'-'+str(rann)+'-tiles.fits')
    t.write(tilefn, overwrite=True)

    tiles = load_tiles(
        tiles_file=tilefn, obsha=obsha, obstheta=obstheta,
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
    
    #print(tile)
    # Load target files...
    tilenofafn = os.path.join(tileoutdir, 'tilenofa-%i.fits' % tile)
    load_target_file(tgs, tagalong, tilenofafn)
    #loading it again straight to table format because I can't quickly figure out exactly where targetid,ra,dec gets stored
    tar_tab = fitsio.read(tilenofafn, columns=tarcols)

    # Find targets within tiles, and project their RA,Dec positions
    # into focal-plane coordinates.
    tile_targetids, tile_x, tile_y, tile_xy_cs5 = targets_in_tiles(hw, tgs, tiles, tagalong)
    # Compute the targets available to each fiber for each tile.
    tgsavail = TargetsAvailable(hw, tiles, tile_targetids, tile_x, tile_y)
    # Compute the fibers on all tiles available for each target and sky
    favail = LocationsAvailable(tgsavail)

    # FAKE stucksky (positioners that happen to be stuck on good sky positions)
    stucksky = {}

    # Create assignment object
    asgn = Assignment(tgs, tgsavail, favail, stucksky)
    tgsavail = asgn.targets_avail()
    avail = tgsavail.tile_data(tile)
    navail = np.sum([len(avail[x]) for x in avail.keys()])
    fibers = dict(hw.loc_fiber)
    fdata = Table()
    fdata['LOCATION'] = np.zeros(navail, dtype=int)
    fdata['FIBER']    = np.zeros(navail, dtype=int)
    fdata['TARGETID'] = np.zeros(navail, dtype=int)

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
    if args.getcoll == 'y':
        coll = asgn.check_avail_collisions(tile)
        kl = np.array(list(coll.keys())).transpose()
        locs = kl[0]
        ids = kl[1]
        locids = ids*10000+locs
        log.info('N collisions: %i' % len(coll))
        locidsin = np.isin(fdata['LOCATION']+10000*fdata['TARGETID'],locids)
        log.info('N collisions original: %i %i' % (np.sum(locidsin),len(fdata)))
        fdata['COLLISION'] = locidsin
    #colltab = Table(forig[locidsin])
    fdata['TILEID'] = tile
    return fdata

def run_one_tile(ind):
    t0 = time.time()
    write_tile_targ(ind)
    res = getpa(ind)
    res = np.array(res)
    t1 = time.time()
    log.info('Tile %i took %.3f sec' % (tiletab[ind]['TILEID'], t1-t0))
    return res

def read_fba_header(ind):
    '''
    Read the fiberassign header for one tile index.
    '''
    tile = tiletab['TILEID'][ind]
    ts = '%06i' % tile
    fbah = fitsio.read_header(os.path.join(desi_input_dir, 'target', 'fiberassign', 'tiles', 'trunk', ts[:3], 'fiberassign-'+ts+'.fits.gz'))
    return dict([(k, fbah[k]) for k in ['RUNDATE', 'MTLTIME', 'FA_RUN', 'FA_HA', 'FIELDROT']])

def parse_datetime(s):
    try:
        return datetime.strptime(s, "%Y-%m-%dT%H:%M:%S%z")
    except ValueError:
        d = datetime.strptime(s, "%Y-%m-%dT%H:%M:%S")
        # msg = "Requested run date '{}' is not timezone-aware.  Assuming UTC.".format(runtime)
        d = d.replace(tzinfo=timezone.utc)

hardware_times = []
def get_hardware_for_time(t):
    global hardware_times
    for tlo,thi,hw in hardware_times:
        if (tlo <= t) and (thi is None or thi > t):
            #print('Match to time range', tlo, 'to', thi)
            return hw
    return None

def main():
    from multiprocessing import Pool
    tls = list(tiletab['TILEID'])
    #inds = np.flatnonzero(np.array(tls) == 1230)
    #inds = np.arange(256)
    inds = np.arange(len(tls))

    t0 = time.time()
    # Read all fiberassign headers to get the RUNDATES.
    with Pool(processes=args.nprocs) as pool:
        headers = pool.map(read_fba_header, inds)
    rundates = set([h['RUNDATE'] for h in headers])
    rundates = sorted(list(rundates))
    log.info('Unique rundates: %i of %i' % (len(rundates), len(headers)))
    t1 = time.time()
    log.info('Reading fiberassign headers in parallel: %.3f sec' % (t1-t0))

    # Read all hardware configurations for our RUNDATES.
    global hardware_times
    for t in rundates:
        dt = parse_datetime(t)
        cached = get_hardware_for_time(dt)
        if cached is not None:
            continue
        hw,time_lo,time_hi = load_hardware(rundate=t, add_margins=margins, get_time_range=True)
        hardware_times.append((time_lo, time_hi, hw))

    t2 = time.time()
    log.info('Loading hardware in series: %.3f sec' % (t2-t1))

    # Keeping this old code because it's a little easier to understand than what we're doing
    # below (streaming results to disk).
    #
    # # Run fiber assignment on tiles in parallel
    # with Pool(processes=128) as pool:
    #     res = pool.map(run_one_tile, inds)
    # t3 = time.time()
    # log.info('Running tiles in parallel: %.3f sec' % (t3-t2))
    #
    # # Merge and write results
    # colltot = np.concatenate(res)
    # if args.getcoll == 'y':
    #     print(len(colltot),np.sum(colltot['COLLISION']))
    # t3b = time.time()
    # 
    # common.write_LSS(colltot,paoutdir+'/pota-'+args.prog+'.fits')
    # t4 = time.time()
    # log.info('Merging results and writing: %.3f sec (%.3f + %.3f)' % (t4-t3, t3b-t3, t4-t3b))

    # Write output *while* retrieving results in parallel
    outfn = os.path.join(paoutdir, 'pota-'+args.prog+'.fits')
    tempout = outfn + '.tmp'
    fits = fitsio.FITS(tempout, 'rw', clobber=True)
    first = True
    ntot = 0
    ncoll = 0
    with Pool(processes=args.nprocs) as pool:
        it = pool.imap_unordered(run_one_tile, inds)
        # fetch results as they complete
        for res in it:
            ntot += len(res)
            ncoll += np.sum(res['COLLISION'])
            # First result: write to output file.
            if first:
                fits.write(res, extname='LSS')
                first = False
            # Subsequent results: append to output file.
            else:
                fits[-1].append(res)
            del res
    fits.close()
    os.rename(tempout, outfn)
    log.info('Wrote %s' % outfn)
    t3 = time.time()
    if args.getcoll == 'y':
        log.info('%i %i' % (ntot, ncoll))
    log.info('Running tiles and writing results: %.3f sec' % (t3-t2))

if __name__ == '__main__':
    main()
