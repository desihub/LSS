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
parser.add_argument("--mockname", help="name of mocks", default='ab_secondgen')
parser.add_argument("--input_mockpath", help="full directory path to input mocks",default='')
parser.add_argument("--input_mockfile", help="mock file name",default='')
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
elif args.tracer[:3] == 'BGS':
    tile = 'BRIGHT'

tiletab = Table.read(f'/global/cfs/cdirs/desi/survey/catalogs/{args.survey}/LSS/tiles-{tile}.fits')

#just using 1 random file for now
ranf = f'/global/cfs/cdirs/desi/survey/catalogs/{args.survey}/LSS/rands_intiles_{tile}_nomask_0.fits'
if not os.path.isfile(ranf):
    print('did not find '+ranf+', will make it')
    input_ran = fitsio.read('/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/randoms-allsky-1-0.fits',columns=['RA','DEC'])
    sel_tiles = is_point_in_desi(tiletab,input_ran['RA'],input_ran['DEC'])
    input_ran = input_ran[sel_tiles]
    common.write_LSS_scratchcp(input_ran,ranf)

nran = fitsio.read_header(ranf,ext=1)['NAXIS2']
area = nran/2500
print('the area is '+str(area))

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


if args.mockname == 'ab_secondgen':
    zs = {'ELG':{'z0.950':[0.,1.1], 'z1.325':[1.1,99.]}, 'LRG':{'z0.500':[0.,0.6], 'z0.800':[0.6,99.]}, 'QSO':{'z1.400':[0.,99.]}}
    downsampling = {'ELG':0.7345658717688022, 'LRG':0.708798313382828, 'QSO':0.39728966594530174}
    percentage_elg_hip = 0.1
    mockpath = '/global/cfs/cdirs/desi/cosmosim/SecondGenMocks/AbacusSummit/CutSky'
    file_name = 'cutsky_{TYPE}_{Z}_AbacusSummit_base_c000_ph{PH}.fits'
    real = 0
    def mask_secondgen(nz=0, foot=None, nz_lop=0):
        if foot == 'Y1':
            Y5 = 0
            Y1 = 1
        elif foot == 'Y5':
            Y5 = 1
            Y1 = 0
        else:
            Y5 = 0
            Y1 = 0
        return nz * (2**0) + Y5 * (2**1) + nz_lop * (2**2) + Y1 * (2**3)
    datas = []

    for bins in zs[type_]:
        print(bins)
        thepath = os.path.join(mockpath, type_, bins, file_name.format(TYPE = type_, Z = bins, PH = "%03d" % real))
        print('thepath')
        print(thepath)
        dat = fitsio.read(thepath, columns=['RA','DEC','Z','Z_COSMO','STATUS'])#f[1].data
        mask = (dat['Z']>= zs[type_][bins][0])&(dat['Z']< zs[type_][bins][1])
        datas.append(Table(dat[mask]))
    data = vstack(datas)
#elif conditions could be added here to properly process other kinds of inputs
elif args.mockname == 'uchuu' and args.tracer == 'BGS_ANY':
    data = Table(fitsio.read(args.input_mockpath+args.input_mockfile,columns=['ra','dec']))
    data.rename_column('ra','RA')
    data.rename_column('dec','DEC')
elif args.input_mockfile[:-2] == 'h5':
    import h5py
    import hdf5plugin #need to be in the cosmodesi test environment, as of Sep 4th 25
    data = Table()
    with h5py.File(args.input_mockpath+args.input_mockfile) as fn:
        columns = fn.keys()
        for col in columns:
            data[col] = fn[col][:]
else:
    data = Table.read(args.input_mockpath+args.input_mockfile)

#cut data to area in tiles file
sel_tiles = is_point_in_desi(tiletab,data['RA'],data['DEC'])
data = data[sel_tiles]
print(len(data),' in tiles area')

#downsampling needed for abacus
if args.mockname == 'ab_secondgen':
    status = data['STATUS'][()]
    idx = np.arange(len(status))

    mask_main = mask_secondgen(nz=1, foot='Y5')
    idx_main = idx[(status & (mask_main))==mask_main]

    if type_ == 'LRG' or type_ == 'QSO':
        ran_tot = np.random.uniform(size = len(idx_main))
        idx_main = idx_main[(ran_tot<=downsampling[type_])]
        data = data[idx_main]
        data = Table(data)
        data['DESI_TARGET'] = desitar[type_]
        data['PRIORITY_INIT'] = priority[type_]
        data['PRIORITY'] = priority[type_]
        data['NUMOBS_MORE'] = numobs[type_]
        data['NUMOBS_INIT'] = numobs[type_]

    else:
        #abacus 2nd gen has this selection defined to split LOP/VLO
        datat = []
        mask_LOP = mask_secondgen(nz=1, foot='Y5', nz_lop=1)
        idx_LOP = idx[(status & (mask_LOP))==mask_LOP]
        idx_VLO = np.setdiff1d(idx_main, idx_LOP)

        ran_lop = np.random.uniform(size = len(idx_LOP))
        idx_LOP = idx_LOP[(ran_lop<=downsampling[type_])]
        ran_vlo = np.random.uniform(size = len(idx_VLO))
        idx_VLO = idx_VLO[(ran_vlo<=downsampling[type_])]

        data_lop = Table(data[idx_LOP])
        data_vlo = Table(data[idx_VLO])
        
        df_lop=data_lop.to_pandas()
        df_vlo=data_vlo.to_pandas()

        num_HIP_LOP = int(len(df_lop) * percentage_elg_hip)
        df_HIP_LOP = df_lop.sample(n=num_HIP_LOP)
        remaining_LOP = df_lop.drop(df_HIP_LOP.index)
        df_HIP_LOP.reset_index(drop=True, inplace=True)
        remaining_LOP.reset_index(drop=True, inplace=True)

        num_HIP_VLO = int(len(df_vlo) * percentage_elg_hip)
        df_HIP_VLO = df_vlo.sample(n=num_HIP_VLO)
        remaining_VLO = df_vlo.drop(df_HIP_VLO.index)
        df_HIP_VLO.reset_index(drop=True, inplace=True)
        remaining_VLO.reset_index(drop=True, inplace=True)

        remaining_LOP['PRIORITY_INIT'] = 3100
        remaining_LOP['PRIORITY'] = 3100
        remaining_LOP['DESI_TARGET'] = 2**5 + 2**1


        remaining_VLO['PRIORITY_INIT'] = 3000
        remaining_VLO['PRIORITY'] = 3000
        remaining_VLO['DESI_TARGET'] = 2**7 + 2**1

        df_HIP_LOP['PRIORITY_INIT'] = 3200
        df_HIP_LOP['PRIORITY'] = 3200
        df_HIP_LOP['DESI_TARGET'] = 2**6 + 2**1 + 2**5

        df_HIP_VLO['PRIORITY_INIT'] = 3200
        df_HIP_VLO['PRIORITY'] = 3200
        df_HIP_VLO['DESI_TARGET'] = 2**6 + 2**1 + 2**7

        remaining_LOP['NUMOBS_MORE'] = numobs[type_]
        remaining_LOP['NUMOBS_INIT'] = numobs[type_]
        remaining_VLO['NUMOBS_MORE'] = numobs[type_]
        remaining_VLO['NUMOBS_INIT'] = numobs[type_]
        df_HIP_LOP['NUMOBS_MORE'] = numobs[type_]
        df_HIP_LOP['NUMOBS_INIT'] = numobs[type_]
        df_HIP_VLO['NUMOBS_MORE'] = numobs[type_]
        df_HIP_VLO['NUMOBS_INIT'] = numobs[type_]

        datat.append(Table.from_pandas(remaining_LOP))
        datat.append(Table.from_pandas(remaining_VLO))
        datat.append(Table.from_pandas(df_HIP_LOP))
        datat.append(Table.from_pandas(df_HIP_VLO))
        data = vstack(datat)
        del datat

else:
    data['DESI_TARGET'] = desitar[type_[:3]]
    data['PRIORITY_INIT'] = priority[type_[:3]]
    data['PRIORITY'] = priority[type_[:3]]
    if type_ == 'ELG':
        if args.ELGsplit == 'y':
            sel_LOP = data[args.ELGtpcol] == 1 #assuming that the input mocks have set the LOP sample to 1 for this column
            data['DESI_TARGET'][~sel_LOP] = 2+2**7
            data['PRIORITY_INIT'][~sel_LOP] = 3000
            data['PRIORITY'][~sel_LOP] = 3000
            data['DESI_TARGET'][sel_LOP] = 2+2**5
            data['PRIORITY_INIT'][sel_LOP] = 3100
            data['PRIORITY'][sel_LOP] = 3100

        rans = rng.random(len(data))
        sel_HIP = rans < 0.1 #10% of ELG get promoted to HIP
        data['DESI_TARGET'][sel_HIP] += 2**6
        data['PRIORITY_INIT'][sel_HIP] = 3200
        data['PRIORITY'][sel_HIP] = 3200
        print('ELG priorities',str(np.unique(data['PRIORITY'],return_counts=True)))
    


obstype = type_
if type_ == 'ELG':
    obstype= 'ELGnotqso'

obsdata = fitsio.read(f'/global/cfs/cdirs/desi/survey/catalogs/{args.survey}/LSS/{args.specdata}/LSScats/{args.dataversion}/{obstype}_full_HPmapcut.dat.fits',columns=['DESI_TARGET'])
obsarea = fitsio.read_header(f'/global/cfs/cdirs/desi/survey/catalogs/{args.survey}/LSS/{args.specdata}/LSScats/{args.dataversion}/{obstype}_0_full_HPmapcut.ran.fits',ext=1)['NAXIS2']/2500

obsdens = len(obsdata)/obsarea

mdens = len(data)/area

print('The mock target density for '+type_+' is '+str(round(mdens,3)))
print('The data target density for '+type_+' is '+str(round(obsdens,3)))

if type_ == 'ELG' and args.ELGsplit == 'y':
    sel_mock = (data['DESI_TARGET'] & 2**7) > 0
    sel_obs = (obsdata['DESI_TARGET'] & 2**7) > 0
    obsdens = np.sum(sel_obs)/obsarea
    mdens = np.sum(sel_mock)/area
    print('The mock target density for ELG_VLO is '+str(round(mdens,3)))
    print('The data target density for ELG_VLO is '+str(round(obsdens,3)))
    sel_mock = (data['DESI_TARGET'] & 2**5) > 0
    sel_obs = (obsdata['DESI_TARGET'] & 2**5) > 0
    obsdens = np.sum(sel_obs)/obsarea
    mdens = np.sum(sel_mock)/area
    print('The mock target density for ELG_LOP is '+str(round(mdens,3)))
    print('The data target density for ELG_LOP is '+str(round(obsdens,3)))

#put in something to make plots as function of z
