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

import LSS.common_tools as common
from LSS.imaging import get_pixel_bitmasknobs as bitmask #get_nobsandmask
from LSS.main.cattools import count_tiles_better
from LSS.globals import main


parser = argparse.ArgumentParser()
parser.add_argument("--survey", help="e.g., Y1, DA2",default='DA2')
parser.add_argument("--input_mockpath", help="full directory path to input mocks",default='')
parser.add_argument("--input_mockfile", help="mock file name",default='')
parser.add_argument("--output_fullpathfn", help="output mock file and full path",default='')
parser.add_argument("--nproc", help="number of processors for multiprocessing",default=128,type=int)

args = parser.parse_args()
#DARK tiles from Y3
tiletab = Table.read('/global/cfs/cdirs/desi/survey/catalogs/'+args.survey+'/LSS/tiles-DARK.fits')

##put in cases for particular input file names here
if args.input_mockfile == 'uchuu-desi-y1_v2_0p65_':
    regl = ['NGC','SGC']
    ##output from James code (for LRGs)
    
    tars1 = Table.read(args.input_mockpath+args.input_mockfile+"NGC.fits")
    tars2 = Table.read(args.input_mockpath+args.input_mockfile+"SGC.fits")
    tars1["GALCAP"] = "N"
    tars2["GALCAP"] = "S"
    data = vstack([tars1, tars2])
    print(data.dtype.names)
    data = Table(data)


##this has to be changed for BGSs or QSOs (see original code to check the numbers for each tracer)
desitar = {'LRG':1}
type_ = 'LRG'
priority = {'LRG':3200}
numobs = {'LRG':2}

                
data['DESI_TARGET'] = desitar[type_]
data['PRIORITY_INIT'] = priority[type_]
data['PRIORITY'] = priority[type_]
data['NUMOBS_MORE'] = numobs[type_]
data['NUMOBS_INIT'] = numobs[type_]
targets = data
del data
targets['TARGETID'] = np.random.permutation(np.arange(1,len(targets)+1))
print(len(targets),' in Y5 area')
selY3 = is_point_in_desi(tiletab,targets['RA'],targets['DEC'])
targets = targets[selY3]
print(len(targets),' in Y3 area')
print('getting nobs and mask bits')

if 'BRICKID' not in targets.colnames:
	from desiutil import brick
	tmp = brick.Bricks(bricksize=0.25)
	targets['BRICKID'] = tmp.brickid(targets['RA'], targets['DEC'])

def wrapper(bid_index):

    idx = bidorder[bidcnts[bid_index]:bidcnts[bid_index+1]]
    brickid = bid_unique[bid_index]

    ra, dec = targets['RA'][idx], targets['DEC'][idx]
    tid = targets['TARGETID'][idx]
    bitmask,nobsg,nobsr,nobsz = bitmask.bitmask_radec(brickid, ra, dec)

    data = Table()
    data['idx'] = idx
    data['MASKBITS'] = bitmask
    data['NOBS_G'] = nobsg
    data['NOBS_R'] = nobsr
    data['NOBS_Z'] = nobsz
    data['TARGETID'] = tid

    return data

# Just some tricks to speed up things up
bid_unique, bidcnts = np.unique(targets['BRICKID'], return_counts=True)
bidcnts = np.insert(bidcnts, 0, 0)
bidcnts = np.cumsum(bidcnts)
bidorder = np.argsort(targets['BRICKID'])

# start multiple worker processes
with Pool(processes=int(args.nproc)) as pool:
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
mainp = main(tp = 'LRG', specver = 'kibo')
targets = common.cutphotmask(targets, bits=mainp.imbits)
    
print('cut targets based on photometric mask')
n=len(targets)
#targets.rename_column('Z_COSMO', 'TRUEZ') 
targets.rename_column('Z', 'RSDZ') 
targets['BGS_TARGET'] = np.zeros(n, dtype='i8')
targets['MWS_TARGET'] = np.zeros(n, dtype='i8')
targets['SUBPRIORITY'] = np.random.uniform(0, 1, n)
targets['BRICKNAME'] = np.full(n, '000p0000')    #- required !?!
targets['OBSCONDITIONS'] = obsconditions.mask('DARK') #np.zeros(n, dtype='i8')+int(3) 
targets['SCND_TARGET'] = np.zeros(n, dtype='i8')+int(0)
targets['ZWARN'] = np.zeros(n, dtype='i8')+int(0)


print('4')
#change the name of the output ...
out_file_name = args.output_fullpathfn

common.write_LSS_scratchcp(targets,out_file_name,extname='TARGETS')
print('5')
fits.setval(out_file_name, 'EXTNAME', value='TARGETS', ext=1)
fits.setval(out_file_name, 'OBSCON', value='DARK', ext=1)



