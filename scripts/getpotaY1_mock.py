'''
Find all potential assignment and counts tiles for Y1 mocks
'''

import numpy as np
import os
from astropy.table import Table, join, vstack
import argparse
from fiberassign.hardware import load_hardware
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


parser = argparse.ArgumentParser()
parser.add_argument("--prog", choices=['DARK','BRIGHT'],default='DARK')
parser.add_argument("--mock", default='ab2ndgen')
parser.add_argument("--realization")
parser.add_argument("--getcoll",default='y')
parser.add_argument("--base_output", help="base directory for output",default='/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/')
parser.add_argument("--tracer", help="tracer for CutSky EZ mocks", default=None)
parser.add_argument("--base_input", help="base directory for input for EZ mocks 6Gpc", default = None)
parser.add_argument("--counttiles", default = 'n')

args = parser.parse_args()
if args.mock == 'ab2ndgen':
    #infn = args.base_output+'FirstGenMocks/AbacusSummit/forFA'+args.realization+'_matched_input_full_masknobs.fits'
    infn = args.base_output+'SecondGenMocks/AbacusSummit/forFA'+args.realization+'.fits'
    tars = fitsio.read(infn)
    tarcols = list(tars.dtype.names)
    #tileoutdir = args.base_output+'SecondGenMocks/AbacusSummit/tartiles'+args.realization+'/'
    tileoutdir = os.getenv('SCRATCH')+'/SecondGenMocks/AbacusSummit/tartiles'+args.realization+'/'
    if not os.path.exists(tileoutdir):
        os.makedirs(tileoutdir)
    paoutdir = args.base_output+'SecondGenMocks/AbacusSummit/mock'+args.realization+'/'
elif args.mock == 'ezmocks6':
#     #tr = args.tracer
#     rz = args.realization
#     print("Doing %s"%tr)

#     if  tr == "LRG":
#         infn1 = "/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/EZmock/CutSky_6Gpc/LRG/z0.800/cutsky_LRG_z0.800_EZmock_B6000G1536Z0.8N216424548_b0.385d4r169c0.3_seed%s_NGC.fits"%rz
#         infn2 = "/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/EZmock/CutSky_6Gpc/LRG/z0.800/cutsky_LRG_z0.800_EZmock_B6000G1536Z0.8N216424548_b0.385d4r169c0.3_seed%s_SGC.fits"%rz
#     elif tr == "ELG":
#         infn1 = "/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/EZmock/CutSky_6Gpc/ELG/z1.100/cutsky_ELG_z1.100_EZmock_B6000G1536Z1.1N648012690_b0.345d1.45r40c0.05_seed%s_NGC.fits"%rz
#         infn2 = "/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/EZmock/CutSky_6Gpc/ELG/z1.100/cutsky_ELG_z1.100_EZmock_B6000G1536Z1.1N648012690_b0.345d1.45r40c0.05_seed%s_SGC.fits"%rz
#     elif tr == "QSO":
#         infn1 = "/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/EZmock/CutSky_6Gpc/QSO/z1.400/cutsky_QSO_z1.400_EZmock_B6000G1536Z1.4N27395172_b0.053d1.13r0c0.6_seed%s_NGC.fits"%rz
#         infn2 = "/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/EZmock/CutSky_6Gpc/QSO/z1.400/cutsky_QSO_z1.400_EZmock_B6000G1536Z1.4N27395172_b0.053d1.13r0c0.6_seed%s_SGC.fits"%rz
#    # infn1 = "/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/EZmock/CutSky_6Gpc/LRG/z0.800/cutsky_LRG_z0.800_EZmock_B6000G1536Z0.8N216424548_b0.385d4r169c0.3_seed1_NGC.fits"
#    # infn2 = "/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/EZmock/CutSky_6Gpc/LRG/z0.800/cutsky_LRG_z0.800_EZmock_B6000G1536Z0.8N216424548_b0.385d4r169c0.3_seed1_SGC.fits"
#     tars1 = Table.read(infn1)#fitsio.read(infn1)
#     tars2 = Table.read(infn2)#fitsio.read(infn2)
#     tars1["GALCAP"] = "N"
#     tars2["GALCAP"] = "S"
#     tars = vstack([tars1, tars2])
#     tars['TARGETID'] = np.arange(len(tars))
    
    infn = args.base_input + 'EZMocks_6Gpc/EZMocks_6Gpc_' + args.realization + '.fits'
    tars = fitsio.read(infn)
    tarcols = list(tars.dtype.names)#['TARGETID','RA','DEC', 'Z','Z_COSMO','GALCAP', 'NZ', 'RAW_NZ']

    tileoutdir = args.base_output+'EZMocks_6Gpc/tartiles'+args.realization+'/'
    paoutdir = args.base_output+'EZMocks_6Gpc/EzMocks/mock'+args.realization+'/'
    if args.tracer is not None:
        tileoutdir += args.tracer+'/'
        paoutdir += args.tracer+'/'

print(tars.dtype.names)

if not os.path.exists(tileoutdir):
    os.makedirs(tileoutdir)
    print('made '+tileoutdir)
if not os.path.exists(paoutdir):
    os.makedirs(paoutdir)
    print('made '+paoutdir)


tiletab = Table.read('/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/tiles-'+args.prog+'.fits')


def write_tile_targ(inds ):
    tiles = tiletab[inds]
    #for i in range(0,len(tiles)):
    fname = tileoutdir+'/tilenofa-'+str(tiles['TILEID'])+'.fits'
    print('creating '+fname)
    tdec = tiles['DEC']
    decmin = tdec - trad
    decmax = tdec + trad
    wdec = (tars['DEC'] > decmin) & (tars['DEC'] < decmax)
    #print(len(rt[wdec]))
    inds = desimodel.footprint.find_points_radec(tiles['RA'], tdec,tars[wdec]['RA'], tars[wdec]['DEC'])
    print('got indexes')
    rtw = tars[wdec][inds]
    rmtl = Table(rtw)
    print('made table for '+fname)
    del rtw
    #n=len(rmtl)
    #rmtl['TARGETID'] = np.arange(1,n+1)+10*n*rannum 
    #rmtl['TARGETID'] = np.arange(len(rmtl))
    #print(len(rmtl['TARGETID'])) #checking this column is there
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
    print('added columns for '+fname)
    rmtl.write(fname,format='fits', overwrite=True)
    del rmtl
    print('added columns, wrote to '+fname)
    #nd += 1
    #print(str(nd),len(tiles))



margins = dict(pos=0.05,
                   petal=0.4,
                   gfa=0.4)


log = Logger.get()
rann = 0
n = 0


def getpa(ind):

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
    
    print(tile)
    # Load target files...
    load_target_file(tgs, tagalong, tileoutdir+'/tilenofa-%i.fits' % tile)
    #loading it again straight to table format because I can't quickly figure out exactly where targetid,ra,dec gets stored
    tar_tab = fitsio.read(tileoutdir+'/tilenofa-%i.fits' % tile,columns =tarcols)

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
    if args.getcoll == 'y':
        coll = asgn.check_avail_collisions(tile)
        kl = np.array(list(coll.keys())).transpose()
        locs = kl[0]
        ids = kl[1]
        locids = ids*10000+locs
        print('N collisions:', len(coll))
        locidsin = np.isin(fdata['LOCATION']+10000*fdata['TARGETID'],locids)
        print('N collisions original:',np.sum(locidsin),len(fdata))
        fdata['COLLISION'] = locidsin
    #colltab = Table(forig[locidsin])
    fdata['TILEID'] = tile
    return fdata
    
if __name__ == '__main__':
    from multiprocessing import Pool
    tls = list(tiletab['TILEID'])#[:10])
    inds = np.arange(len(tls))
    #write_tile_targ(inds[0])
    with Pool(processes=128) as pool:
        res = pool.map(write_tile_targ, inds)

    with Pool(processes=128) as pool:
        res = pool.map(getpa, inds)
    colltot = np.concatenate(res)
    if args.getcoll == 'y':
        print(len(colltot),np.sum(colltot['COLLISION']))
    
    common.write_LSS(colltot,paoutdir+'/pota-'+args.prog+'.fits')
        
    