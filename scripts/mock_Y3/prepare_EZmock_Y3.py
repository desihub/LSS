#make sure add the LSS repo to your python path
from astropy.io import fits # Access to FITS (Flexible Image Transport System) files.
from astropy.table import Table, hstack, vstack, Column # A class to represent tables of heterogeneous data.
import fitsio
import numpy as np
import os
import argparse
import sys
import json
from desitarget.targetmask import obsconditions
from desimodel.footprint import is_point_in_desi
from multiprocessing import Pool
import time
import LSS.common_tools as common
from LSS.imaging import get_pixel_bitmasknobs as bitmask #get_nobsandmask
from LSS.main.cattools import count_tiles_better
from LSS.globals import main



parser = argparse.ArgumentParser()
parser.add_argument("--survey", help="e.g., Y1, DA2",default='DA2')
parser.add_argument("--input_mockpath", help="full directory path to input mocks",default='')
parser.add_argument("--input_mockfile", help="mock file name",default='')
parser.add_argument("--output_fullpathfn", help="output mock file and full path",default='')
parser.add_argument("--nproc", help="number of processors for multiprocessing",default=128)
parser.add_argument("--tracer", help="LRG, ELG or QSO",default='LRG')
parser.add_argument("--zcol", help="name of column with redshift, including RSD",default='Z')
parser.add_argument("--ELGsplit", help="Are the ELGs split into LOP and VLO? If 'n', assuming all LOP",default='y')
parser.add_argument("--ELGtpcol", help="column distinguishing the ELG type; assumed boolean with True being LOP",default='LOP')
parser.add_argument("--ran_seed", help="seed for randoms; make sure this is different if running many in parallel",default=10)


args = parser.parse_args()

rng = np.random.default_rng(seed=int(args.ran_seed))

if args.tracer in ['LRG', 'QSO', 'ELG']:
    tile = 'DARK'
elif args.tracer == 'BGS':
    tile = 'BRIGHT'

tiletab = Table.read(f'/global/cfs/cdirs/desi/survey/catalogs/{args.survey}/LSS/tiles-{tile}.fits')


##tars = Table.read(args.input_mockpath+args.input_mockfile+".fits")
data = Table.read(args.input_mockpath+args.input_mockfile, format="fits")

# Adding the WEIGHT column
#tars['WEIGHT'] = np.ones(tars['RA'].shape[0])
data['WEIGHT'] = np.ones(data['RA'].shape[0])

# Conditions for NGC and SGC; AJR does not think this is necessary
#condN = common.splitGC(tars)#(tars['RA'] > 85) & (tars['RA'] < 302)
#condS = ~condN#(tars['RA'] < 85) | (tars['RA'] > 302)

# Splitting the DataFrame; AJR: is there a reason for this here? 
#tarsN = tars[condN]
#tarsS = tars[condS]

#tarsN["GALCAP"] = "N"
#tarsS["GALCAP"] = "S"

#data = vstack([tarsN, tarsS])
#data = Table(tars)#Table(data)
#del tars

desitar = {'LRG':1, 'QSO': 4, 'ELG':34,'BGS': int(2**60)}
priority = {'LRG':3200, 'QSO':3400, 'ELG':3100,'ELG_VOL':3000,'ELG_HIP':3200,'BGS':2100}
numobs = {'LRG':1, 'ELG':1, 'QSO':1, 'BGS':1}
norm = {'LRG':1, 'ELG':1, 'QSO':1, 'BGS':2**60}

type_ = args.tracer

mockname = 'EZmock'


                
data['DESI_TARGET'] = desitar[type_]
data['PRIORITY_INIT'] = priority[type_]
data['PRIORITY'] = priority[type_]
# if type_ == 'ELG' and args.ELGsplit == 'y':
#     sel_LOP = data[args.ELGtpcol] == 1
#     data['DESI_TARGET'][~sel_LOP] = 2+2**7
#     data['PRIORITY_INIT'][~sel_LOP] = 3000
#     data['PRIORITY'][~sel_LOP] = 3000
#     rans = rng.random(len(data))
#     sel_HIP = rans < 0.1 #10% of ELG get promoted to HIP
#     data['DESI_TARGET'][sel_HIP] += 2**6
#     data['PRIORITY_INIT'][sel_HIP] = 3200
#     data['PRIORITY'][sel_HIP] = 3200
#     print('ELG priorities',str(np.unique(data['PRIORITY'],return_counts=True)))

if type_ == 'ELG' and mockname == 'EZmock':
    rans = rng.random(len(data))
    sel_HIP = rans < 0.1 #10% of ELG get promoted to HIP
    data['DESI_TARGET'][sel_HIP] += 2**6
    data['PRIORITY_INIT'][sel_HIP] = 3200
    data['PRIORITY'][sel_HIP] = 3200
    print('ELG priorities',str(np.unique(data['PRIORITY'],return_counts=True)))
    


data['NUMOBS_MORE'] = numobs[type_]
data['NUMOBS_INIT'] = numobs[type_]
if type_ == 'QSO':
    sel_highz = data[args.zcol] > 2.1
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

# ZD add EZmock Y3 mask
ez_y3_mask = (targets['STATUS']&2)>0
targets = targets[ez_y3_mask]
print(len(targets),' in Y3 area')
print('getting nobs and mask bits')

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



    # Just some tricks to speed up things up
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

mainp = main(tp = type_, specver = 'kibo')
targets = common.cutphotmask(targets, bits=mainp.imbits)

print('cut targets based on photometric mask')

## for EZmock dark tracers, select targets with STATUS==2  !!!!
# status_mask = (targets['STATUS'] == 2)     ## Mark subsamples matching Y1 n(z)
# targets = targets[status_mask]


n=len(targets)
#targets.rename_column('Z_COSMO', 'TRUEZ') 
targets.rename_column(args.zcol, 'RSDZ') 
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
common.write_LSS_scratchcp(targets,out_file_name,extname='TARGETS')
fits.setval(out_file_name, 'EXTNAME', value='TARGETS', ext=1)
fits.setval(out_file_name, 'OBSCON', value=tile, ext=1)


