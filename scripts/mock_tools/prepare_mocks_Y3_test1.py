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
#from multiprocessing import Pool
#import time
import LSS.common_tools as common
from LSS.imaging import get_pixel_bitmasknobs as bitmask #get_nobsandmask
#from LSS.main.cattools import count_tiles_better
from LSS.globals import main

from calibrate_nz_prep import calibrate_nz


from numpy.random import Generator, PCG64
rng = Generator(PCG64())




parser = argparse.ArgumentParser()
parser.add_argument("--survey", help="e.g., Y1, DA2",default='DA2')
parser.add_argument("--specdata", help="mountain range for spec prod",default='loa-v1')
parser.add_argument("--mockname", help="name of mocks: holimock, EZmock, abacushf, uchuu") #, default='EZmock')
parser.add_argument("--input_mockpath", help="full directory path to input mocks",default='')
parser.add_argument("--input_mockfile", help="mock file name",default='')
parser.add_argument("--output_fullpathfn", help="output mock file and full path",default='')
#parser.add_argument("--nproc", help="number of processors for multiprocessing",default=128)
parser.add_argument("--tracer", help="LRG, ELG or QSO") #,default='LRG')
parser.add_argument("--ztruecol", help="name of column with true redshift in the input catalog")
parser.add_argument("--zrsdcol", help="name of column with redshift, including RSD")
parser.add_argument("--need_footprint", help="Do we need to prune by rounded tile footprint?", default='y')

parser.add_argument("--need_nz_calib", help="Do we need to calibrate n(z) to match that of the survey?", default='n')


#parser.add_argument("--ELGsplit", help="Are the ELGs already splitted into LOP and VLO? If 'n', make the division",default='y')
parser.add_argument("--ELGtpcol", help="column distinguishing the ELG type; assumed boolean with True being LOP, if mockname not abacushf",default='LOP')

#parser.add_argument("--ran_seed", help="seed for randoms; make sure this is different if running many in parallel",default=10)
#parser.add_argument("--EZnzmask", help="apply mask on galaxy number density n(z)",default='n')

args = parser.parse_args()


if args.tracer in ['LRG', 'QSO', 'ELG']:
    tile = 'DARK'
elif args.tracer == 'BGS':
    tile = 'BRIGHT'


def read_parent_mock(mockname, filename):

    if mockname == 'holimock':
        import h5py
        
        data = Table(h5py.File(filename, 'r+'))


#        data = np.load(filename)

#        data_dict = dict(data)
#        data_dict['RA'] = data_dict.pop('ra')  # Rename key
#        data_dict['DEC'] = data_dict.pop('dec')  # Rename key
#        data = Table(data_dict)
#        del data_dict
    else:
        data = Table.read(filename)

    return data


if args.mockname == 'abacushf':
    args.need_footprint = 'n'
    args.need_nz_calib = 'n'

desitar = {'LRG':desi_mask.LRG, 'QSO': desi_mask.QSO, 'ELG':desi_mask.ELG, 'BGS': desi_mask.BGS_ANY}
#AURE PRIORITY for BGS, is BRIGHT or what??? BGS_ANY? Shall we make same as ELG?

priority = {'LRG': desi_mask.LRG.priorities['UNOBS'],
            'QSO': desi_mask.QSO.priorities['UNOBS'],
            'ELG': desi_mask.ELG.priorities['UNOBS'],
            'BGS': bgs_mask.BGS_FAINT.priorities['UNOBS']}

#AURE            'ELG_VLO': desi_mask.ELG.priorities['UNOBS'],
#    AURE        'ELG_LOP': desi_mask.ELG_LOP.priorities['UNOBS'],
#    AURE        'ELG_HIP': desi_mask.LRG.priorities['UNOBS'], # same priority as LRG
#    AURE        'BGS': bgs_mask.BGS_FAINT.priorities['UNOBS']}  #AURE change this to FAINT instead of bright

numobs = {'LRG': desi_mask.LRG.numobs,
          'ELG': desi_mask.ELG.numobs,
          'QSO': desi_mask.QSO.numobs,
          'BGS': bgs_mask.BGS_BRIGHT.numobs}

norm = {'LRG':1, 'ELG':1, 'QSO':1, 'BGS':2**60}   #AURE CONFIRM IS CORRECT


def add_basic_columns(input_data, tracercol):
    print('Adding basic columns')
    # Adding the WEIGHT column
    if 'WEIGHT' not in input_data.colnames:
        input_data['WEIGHT'] = np.ones(input_data['RA'].shape[0])

    input_data['DESI_TARGET'] = desitar[tracercol]
    input_data['PRIORITY_INIT'] = priority[tracercol]
    input_data['PRIORITY'] = priority[tracercol]
    input_data['NUMOBS_MORE'] = numobs[tracercol]
    input_data['NUMOBS_INIT'] = numobs[tracercol]

    if ('TRUEZ' not in input_data.colnames) and (args.ztruecol != None):
        input_data.rename_column(args.ztruecol, 'TRUEZ')
    if ('RSDZ' not in input_data.colnames) and (args.zrsdcol != None):
        input_data.rename_column(args.zrsdcol, 'RSDZ')


    n = len(input_data)
    
    input_data['MWS_TARGET'] = np.zeros(n, dtype='i8')
    input_data['SUBPRIORITY'] = np.random.uniform(0, 1, n)
    input_data['BRICKNAME'] = np.full(n, '000p0000')    #- required !?!
    input_data['OBSCONDITIONS'] = obsconditions.mask(tile) #np.zeros(n, dtype='i8')+int(3) 
    input_data['SCND_TARGET'] = np.zeros(n, dtype='i8')+int(0)
    input_data['ZWARN'] = np.zeros(n, dtype='i8')+int(0)
    input_data['BGS_TARGET'] = np.zeros(n, dtype='i8')

    return input_data, n


def select_footprint(input_data, todo):

    if todo == 'n':
        size = len(input_data)
        print('no footprint selection')
    elif todo == 'y':
        print('imposing approximate tile footprint', args.survey)
        from desimodel.footprint import is_point_in_desi
        
        tiletab = Table.read(f'/global/cfs/cdirs/desi/survey/catalogs/{args.survey}/LSS/tiles-{tile}.fits')

        selFOOT = is_point_in_desi(tiletab, input_data['RA'], input_data['DEC'])
        input_data = input_data[selFOOT]
        size = len(input_data)
    else:
        raise Exception('You should tell me if I need to apply footprint or not')

    return input_data, size


tracer = args.tracer


print('STARTING PREPATION')
print('tracer type', tracer)
print('mock flavour is ', args.mockname)
print('it will read the file', os.path.join(args.input_mockpath, args.input_mockfile))

data = read_parent_mock(args.mockname, os.path.join(args.input_mockpath, args.input_mockfile))
data, initial_size = add_basic_columns(data, tracer)

print('size of sample at start of the pipeline', initial_size)

data, after_footprint_size = select_footprint(data, args.need_footprint)

print('size of sample after footprint pruning (if needed)', after_footprint_size)

if args.need_nz_calib == 'y':
    print('Start nz calibration to match density and nz to the data adding STATUS columns')
    data = calibrate_nz(data, redshift_column = 'RSDZ', tracer_type=tracer, survey=args.survey)
    
if tracer == 'LRG':
    print('start specific LRG preparation')
    if args.mockname == 'abacushf':
        sel = (data['STATUS'] & 33 == 33)
    elif args.mockname == 'EZmock':
        sel = (data['STATUS'] & 2 == 2)
    else:
        if 'STATUS' not in data.colnames:
            data['STATUS'] = np.ones(len(data), dtype=np.int32)
        sel = (data['STATUS'] & 1 == 1)
    
    idx = np.arange(len(data))
    
    data = data[idx[sel]]
    size_final = len(data)
    print('After selecting LRG according to nz matching, size =', size_final)
    del data['STATUS']

if tracer == 'QSO':
    print('start specific QSO preparation')
    if args.mockname == 'abacushf':
        sel = (data['STATUS'] & 33 == 33)

    elif args.mockname == 'EZmock':
        sel = (data['STATUS'] & 2 == 2)

    else:
        if 'STATUS' not in data.colnames:
            data['STATUS'] = np.ones(len(data), dtype=np.int32)

        sel = (data['STATUS'] & 1 == 1)

    idx = np.arange(len(data))
    
    data = data[idx[sel]]
    size_final = len(data)

    print('After selecting LRG according to nz matching, size =', size_final)
    sel_highz = data['RSDZ'] > 2.1
    data['NUMOBS_MORE'][sel_highz] = 4
    data['NUMOBS_INIT'][sel_highz] = 4
    del data['STATUS']

if tracer == 'ELG':
    print('start specific ELG preparation')
    datat = []
    idx = np.arange(len(data))
    if args.mockname == 'abacushf':

        from calibrate_nz_prep import mask_abacusHF
        status = data['STATUS'][()]
        mask_main = mask_abacusHF(nz=1, foot=args.survey)
        idx_main = idx[(status & (mask_main))==mask_main]

        mask_LOP = mask_abacusHF(nz=1, foot=args.survey, nz_lop=1)
        idx_LOP = idx[(status & (mask_LOP))==mask_LOP]
            
        idx_VLO = np.setdiff1d(idx_main, idx_LOP)
        data_lop = Table(data[idx_LOP])
        data_vlo = Table(data[idx_VLO])


    elif args.mockname == 'EZmock':
        import sample_elg_ezmock as se
        lop, vlo = se.create_subsample(data)
        
        data_lop = Table.from_pandas(lop)
        data_vlo = Table.from_pandas(vlo)


    else:
        if 'STATUS' not in data.colnames:
            idx_LOP = idx[data[args.ELGtpcol] == 1]
            idx_VLO = np.setdiff1d(idx, idx_LOP)

        else:
            status = data['STATUS'][()]
            idx_main = idx[(status & 1)==1]
            idx_LOP = idx[(status & 4)==4]
            idx_VLO = np.setdiff1d(idx_main, idx_LOP)
    
        data_lop = Table(data[idx_LOP])
        data_vlo = Table(data[idx_VLO])


    data_lop['DESI_TARGET'] += 2**5

    data_vlo['PRIORITY_INIT'] = 3000
    data_vlo['PRIORITY'] = 3000
    data_vlo['DESI_TARGET'] += 2**7 
    
    datat.append(data_lop)
    datat.append(data_vlo)
    data = vstack(datat)
    size_lop = len(data_lop)
    size_vlo = len(data_vlo)

    rans = rng.random(len(data))
    sel_HIP = rans < 0.1 #10% of ELG get promoted to HIP
    data['DESI_TARGET'][sel_HIP] += 2**6
    data['PRIORITY_INIT'][sel_HIP] = 3200
    data['PRIORITY'][sel_HIP] = 3200

    size_final = len(data)

    print('After selecting ELG according to nz matching, size =', size_final)
    print('By subtypes, LOP =', size_lop, 'VLO =', size_vlo)
    del data_lop
    del data_vlo
    del data['STATUS']

#BGS missing
if tracer == 'BGS':

    if args.mockname.lower() == 'uchuu':
        mask_bright = (data['BGS_TYPE'] == 'BRIGHT')
        mask_faint  = (data['BGS_TYPE'] == 'FAINT')
        dat_bright  = data[mask_bright]
        dat_faint   = data[mask_faint]
        print('size of BRIGHT', len(dat_bright))
        print('size of FAINT', len(dat_faint))

        dat_bright['BGS_TARGET'] = 2**1
                
        dat_faint['BGS_TARGET'] = 2**0
        
        PromoteFracBGSFaint=0.2
        ran_hip = np.random.uniform(size = len(dat_faint))
        faint_hip_mask = (ran_hip <= PromoteFracBGSFaint)
        
        dat_faint['BGS_TARGET'][faint_hip_mask] += 2**3   # for high-priority BGS faint
    
        dat_faint['PRIORITY_INIT'][~faint_hip_mask] = 2000
        dat_faint['PRIORITY'][~faint_hip_mask] = 2000

        data = vstack([dat_faint, dat_bright])
        print("Unique PRIORITY_INIT: ", np.unique(data['PRIORITY_INIT']))
        print("Unique PRIORITY: ", np.unique(data['PRIORITY']))
        print("Unique BGS_TARGET: ", np.unique(data['BGS_TARGET']))
        print("High_priority_BGS_faint/BGS_faint: ", np.sum(dat_faint['BGS_TARGET']==9)/len(dat_faint))



    if args.mockname == 'EZmock' and args.EZnzmask == 'y':
        nz_mask = (targets['STATUS']&2)>0   # match Y3 LRG/ELG/QSO n(z)
        targets = targets[nz_mask]




n=len(data)  ##A Ashley le falta estoo!
data['TARGETID'] = (np.arange(1,n+1)+1e8*desitar[tracer]/norm[tracer]).astype(int) #different tracer types need to have different targetids
print('creating TARGETID', np.min(data['TARGETID']), np.max(data['TARGETID']))
if (tracer == 'BGS') and (args.mockname.lower() == 'uchuu'):
    targets.rename_column('ABSMAG_R', 'R_MAG_ABS')
    
out_file_name = args.output_fullpathfn

#change the name of the output ...



#Change type of RA,DEC
if data['RA'].dtype != np.float64:
    print('Imposing float64 to RA')
    data['RA'] = Column(data['RA'], dtype=np.float64)

if data['DEC'].dtype != np.float64:
    print('Imposing float64 to DEC')
    data['DEC'] = Column(data['DEC'], dtype=np.float64)


common.write_LSS_scratchcp(data, out_file_name, extname='TARGETS')
fits.setval(out_file_name, 'EXTNAME', value='TARGETS', ext=1)
fits.setval(out_file_name, 'OBSCON', value=tile, ext=1)


'''
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

dophot = False
out_file_name = args.output_fullpathfn
if dophot:

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

    mainp = main(tp = tracer, specver = args.specdata)
    targets = common.cutphotmask(targets, bits=mainp.imbits)
else:

    ending = out_file_name.split('.')[-1]
    out_file_name = out_file_name.split('.' + ending)[0] + '_noimagingmask_applied.fits'
'''


