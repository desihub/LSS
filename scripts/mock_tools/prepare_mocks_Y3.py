#make sure add the LSS repo to your python path
from astropy.io import fits # Access to FITS (Flexible Image Transport System) files.
from astropy.table import Table, hstack, vstack, Column # A class to represent tables of heterogeneous data.
import fitsio
import numpy as np
import os, sys
import argparse
import random
import json
from desitarget.targetmask import desi_mask, bgs_mask, obsconditions
from desimodel.footprint import is_point_in_desi
from multiprocessing import Pool
import time
import LSS.common_tools as common
from LSS.imaging import get_pixel_bitmasknobs as bitmask #get_nobsandmask
from LSS.main.cattools import count_tiles_better
from LSS.globals import main



parser = argparse.ArgumentParser()
parser.add_argument("--survey", help="e.g., Y1, DA2",default='DA2')
parser.add_argument("--specdata", help="mountain range for spec prod",default='loa-v1')
parser.add_argument("--dataversion", help="version of LSS catalogs",default='v1.1')
parser.add_argument("--mockname", help="name of mocks", default='EZmock')
parser.add_argument("--input_mockpath", help="full directory path to input mocks",default='')
parser.add_argument("--input_mockfile", help="mock file name",default='')
parser.add_argument("--output_fullpathfn", help="output mock file and full path",default='')
parser.add_argument("--nproc", help="number of processors for multiprocessing",default=128)
parser.add_argument("--tracer", help="LRG, ELG or QSO",default='LRG')
parser.add_argument("--ztruecol", help="name of column with true redshift in the input catalog", default='Z_COSMO')
parser.add_argument("--zrsdcol", help="name of column with redshift, including RSD",default='Z')
parser.add_argument("--ELGsplit", help="Are the ELGs split into LOP and VLO? If 'n', assuming all LOP",default='y')
parser.add_argument("--ELGtpcol", help="column distinguishing the ELG type; assumed boolean with True being LOP",default='LOP')
parser.add_argument("--ran_seed", help="seed for randoms; make sure this is different if running many in parallel",default=10)
parser.add_argument("--nzmask", help="apply mask on galaxy number density n(z)",default='n')

args = parser.parse_args()

rng = np.random.default_rng(seed=int(args.ran_seed))

if args.tracer in ['LRG', 'QSO', 'ELG']:
    tile = 'DARK'
elif args.tracer == 'BGS':
    tile = 'BRIGHT'

tiletab = Table.read(f'/global/cfs/cdirs/desi/survey/catalogs/{args.survey}/LSS/tiles-{tile}.fits')


data = Table.read(args.input_mockpath+args.input_mockfile)

# Adding the WEIGHT column
if 'WEIGHT' not in data.colnames:
    data['WEIGHT'] = np.ones(data['RA'].shape[0])

desitar = {'LRG':desi_mask.LRG, 'QSO': desi_mask.QSO, 'ELG':desi_mask.ELG + desi_mask.ELG_LOP, 'BGS': desi_mask.BGS_ANY}
#priority = {'LRG':3200, 'QSO':3400, 'ELG':3100,'ELG_VOL':3000,'ELG_HIP':3200,'BGS':2100}
priority = {'LRG': desi_mask.LRG.priorities['UNOBS'],
            'QSO': desi_mask.QSO.priorities['UNOBS'],
            'ELG': desi_mask.ELG_LOP.priorities['UNOBS'],
            'ELG_VOL': desi_mask.ELG.priorities['UNOBS'],
            'ELG_HIP': desi_mask.LRG.priorities['UNOBS'], # same priority as LRG
            'BGS': bgs_mask.BGS_BRIGHT.priorities['UNOBS']}

numobs = {'LRG': desi_mask.LRG.numobs,
          'ELG': desi_mask.ELG_LOP.numobs,
          'QSO': desi_mask.QSO.numobs,
          'BGS': bgs_mask.BGS_BRIGHT.numobs}
norm = {'LRG':1, 'ELG':1, 'QSO':1, 'BGS':2**60}

type_ = args.tracer
                
data['DESI_TARGET'] = desitar[type_]
data['PRIORITY_INIT'] = priority[type_]
data['PRIORITY'] = priority[type_]
if type_ == 'ELG':
    if args.ELGsplit == 'y':
        sel_LOP = data[args.ELGtpcol] == 1
        data['DESI_TARGET'][~sel_LOP] = 2+2**7
        data['PRIORITY_INIT'][~sel_LOP] = 3000
        data['PRIORITY'][~sel_LOP] = 3000
    rans = rng.random(len(data))
    sel_HIP = rans < 0.1 #10% of ELG get promoted to HIP
    data['DESI_TARGET'][sel_HIP] += 2**6
    data['PRIORITY_INIT'][sel_HIP] = 3200
    data['PRIORITY'][sel_HIP] = 3200
    print('ELG priorities',str(np.unique(data['PRIORITY'],return_counts=True)))
	
data['NUMOBS_MORE'] = numobs[type_]
data['NUMOBS_INIT'] = numobs[type_]
if type_ == 'QSO':
    sel_highz = data[args.zrsdcol] > 2.1
    data['NUMOBS_MORE'][sel_highz] = 4
    data['NUMOBS_INIT'][sel_highz] = 4
    print('numobs counts',str(np.unique(data['NUMOBS_MORE'],return_counts=True)))
targets = data
n=len(targets)  ##A Ashley le falta estoo!

del data

targets['TARGETID'] = (np.random.permutation(np.arange(1,n+1))+1e8*desitar[type_]/norm[type_]).astype(int) #different tracer types need to have different targetids
print(len(targets),' in Y5 area')
selY3 = is_point_in_desi(tiletab,targets['RA'],targets['DEC'])
targets = targets[selY3]
if args.mockname == 'EZmock' and args.nzmask == 'y':
    y3_nz_mask = (targets['STATUS']&2)>0   # match Y3 LRG/ELG/QSO n(z)
    targets = targets[y3_nz_mask]

print(len(targets),' in Y3 area')
#print('getting nobs and mask bits')

# compare the projected number density between the mock and data (full_HPmap) 
survey_lss_path = f"/dvs_ro/cfs/cdirs/desi/survey/catalogs/{args.survey}/LSS/{args.specdata}/LSScats/{args.dataversion}/"
survey_data_file = survey_lss_path + f"{type_}_full_HPmapcut.dat.fits"
survey_rand_file = survey_lss_path + f"{type_}_0_full_HPmapcut.ran.fits"
data_header = fitsio.read_header(survey_data_file, ext=1)
data_len = data_header['NAXIS2']

rand_header = fitsio.read_header(survey_rand_file, ext=1)
rand_den = 2500.0         # 2500 deg^-2 is the projected number density of randoms for the catalog
survey_skyarea = rand_header['NAXIS2']/rand_den
print("Survey sky area [deg]:", survey_skyarea)
data_nden = data_len/survey_skyarea

def get_mock_skyarea(tiletab, rand_den=2500.0):
    survey_randfile = os.path.join(os.path.dirname(args.output_fullpathfn), f"randoms-{args.survey}-{args.specdata}-nden{rand_den}.fits")
    try: 
        print("Load the pre-selected randoms.")
        Nrand = len(Table.read(survey_randfile))
        
    except FileNotFoundError:    
        print("Load the full-sky randoms and select them within survey footprint.")
        rand_file = "/dvs_ro/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/randoms-allsky-1-0.fits"
        allsky_rand = fitsio.read(rand_file, columns=['RA', 'DEC'])
        survey_mask = is_point_in_desi(tiletab, allsky_rand['RA'], allsky_rand['DEC'])  # select points in the survey footprint
        Nrand = np.sum(survey_mask)
        Table(allsky_rand[survey_mask]).write(survey_randfile, overwrite=True)
    
    sky_area = Nrand/rand_den
    print("Mock skey area [deg]:", sky_area)
    return sky_area
    
mock_skyarea = get_mock_skyarea(tiletab)
mock_nden = len(targets)/mock_skyarea

nden_ratio = mock_nden/data_nden
print("Mock/Data projected density ratio:", nden_ratio)

#if nden_ratio > 1.01:
#    print("Downsample the mock.")
#    random.seed(17)
#    rand_dist = np.random.rand(len(targets))
#    downsample_rate = 1.0/nden_ratio
#    downsample_mask = (rand_dist < downsample_rate) 
#    targets = targets[downsample_mask]
    
#assert nden_ratio > 0.99, "The mock's projected number density is >1% lower than that of data."

def wrapper(bid_index):

    idx = bidorder[bidcnts[bid_index]:bidcnts[bid_index+1]]
    brickid = bid_unique[bid_index]

    ra, dec = targets['RA'][idx], targets['DEC'][idx]
    tid = targets['TARGETID'][idx]
    bitmask2,nobsg,nobsr,nobsz = bitmask.bitmask_radec(brickid, ra, dec)

    data = Table()
    data['idx'] = idx
    data['MASKBITS'] = bitmask2
    data['NOBS_G'] = nobsg
    data['NOBS_R'] = nobsr
    data['NOBS_Z'] = nobsz
    data['TARGETID'] = tid

    return data

if 'MASKBITS' not in targets.colnames:
    if 'BRICKID' not in targets.colnames:
        from desiutil import brick
        tmp = brick.Bricks(bricksize=0.25)
        targets['BRICKID'] = tmp.brickid(targets['RA'], targets['DEC'])
    # Just some tricks to speed up things
    bid_unique, bidcnts = np.unique(targets['BRICKID'], return_counts=True)
    bidcnts = np.insert(bidcnts, 0, 0)
    bidcnts = np.cumsum(bidcnts)
    bidorder = np.argsort(targets['BRICKID'])
    
    # start multiple worker processes
    with Pool(processes=int(args.nproc)) as pool: ##hay que poner un int para que funcione!
        res = pool.map(wrapper, np.arange(len(bid_unique)))
    
    res = vstack(res)
    res.sort('idx')
    res.remove_column('idx')
    print('mask columns added')
    
    maskcols = ['NOBS_G','NOBS_R','NOBS_Z','MASKBITS']
    if np.array_equal(res['TARGETID'],targets['TARGETID']):
        for col in maskcols:
            targets[col] = res[col]
        del res

mainp = main(tp = type_, specver = args.specdata)
targets = common.cutphotmask(targets, bits=mainp.imbits)


print('cut targets based on photometric mask')
n=len(targets)
if 'TRUEZ' not in targets.colnames:
    targets.rename_column(args.ztruecol, 'TRUEZ')
if 'RSDZ' not in targets.colnames:
    targets.rename_column(args.zrsdcol, 'RSDZ')
if args.tracer == 'BGS':
    targets['BGS_TARGET'] = 2
else:	
    targets['BGS_TARGET'] = np.zeros(n, dtype='i8')
targets['MWS_TARGET'] = np.zeros(n, dtype='i8')
targets['SUBPRIORITY'] = np.random.uniform(0, 1, n)
targets['BRICKNAME'] = np.full(n, '000p0000')    #- required !?!
targets['OBSCONDITIONS'] = obsconditions.mask(tile) #np.zeros(n, dtype='i8')+int(3) 
targets['SCND_TARGET'] = np.zeros(n, dtype='i8')+int(0)
targets['ZWARN'] = np.zeros(n, dtype='i8')+int(0)

#change the name of the output ...
out_file_name = args.output_fullpathfn
common.write_LSS_scratchcp(targets, out_file_name, extname='TARGETS')
fits.setval(out_file_name, 'EXTNAME', value='TARGETS', ext=1)
fits.setval(out_file_name, 'OBSCON', value=tile, ext=1)


