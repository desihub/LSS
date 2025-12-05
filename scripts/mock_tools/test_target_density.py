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

def mknz(df, area, bs = 0.01, zmin = 0.01, zmax = 1.6, randens = 2500., zcol='Z'):
    from LSS.tabulated_cosmo import TabulatedDESI
    cosmo = TabulatedDESI()
    dis_dc = cosmo.comoving_radial_distance
    print('area is '+str(area))


    nbin = int((zmax-zmin)*(1+bs/10)/bs)
    zhist = np.histogram(df[zcol],bins=nbin,range=(zmin,zmax)) #,weights=wts)
    
    zreturn, nzreturn = [], []
    for i in range(0,nbin):
        zl = zhist[1][i]
        zh = zhist[1][i+1]
        zm = (zh+zl)/2.
        voli = area/(360.*360./np.pi)*4.*np.pi/3.*(dis_dc(zh)**3.-dis_dc(zl)**3.)
        nbarz =  zhist[0][i]/voli
        zreturn.append(zm)
        nzreturn.append(nbarz)
    return zreturn, nzreturn

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

###ranf = '/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v1.0/randoms/rands_intiles_DARK_0_withimagingmask_withz.fits'
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
'''
def return_north(ra, dec):
    #given a table that already includes RA,DEC, add PHOTSYS column denoting whether
    #the data is in the DECaLS ('S') or BASS/MzLS ('N') photometric region
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    c = SkyCoord(ra* u.deg, dec* u.deg,frame='icrs')
    gc = c.transform_to('galactic')
    sel_ngc = gc.b > 0
            
    seln = dec > 32.375

    sel = seln&sel_ngc
    return sel


randata = Table.read(ranf)
maskN = return_north(randata['RA'], randata['DEC'])

area_N = len(randata[maskN])/2500.

area_S = len(randata[~maskN])/2500.
print(area_N, area_S)
'''
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
elif args.input_mockfile[-2:] == 'h5':
    import h5py
    import hdf5plugin #need to be in the cosmodesi test environment, as of Sep 4th 25
    data = Table()
    with h5py.File(args.input_mockpath+args.input_mockfile) as fn:
        columns = fn.keys()
        for col in columns:
            data[col] = fn[col][:]
else:
    data = Table.read(args.input_mockpath+args.input_mockfile)
#    remove = np.random.choice(len(data), size=13000, replace=False)
#    keep = np.setdiff1d(np.arange(len(data)), remove)
#    data = data[keep]

#cut data to area in tiles file
print(len(data),' in all areas')
sel_tiles = is_point_in_desi(tiletab,data['RA'],data['DEC'])
data = data[sel_tiles]
print(len(data),' in tiles area')
print(data.columns)
for name, col in data.columns.items():
    print(f"{name}: {col.dtype}")
#print(data.columns.dtypes)


#data.write('HOMe_%s_zrsd_in_DR2.fits')
'''
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
    
'''

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


import plotext as plt
#do plot n(z)
if type_ == 'LRG':
    z, nz = np.loadtxt('/pscratch/sd/a/acarnero/codes/desi-cutsky-mock/nz_files/DA2/nz_lrg_da2_mocks.txt' , unpack = True)
    
    ztarget, nztarget = mknz(data, area, zmax=1.6, zcol=args.zrsdcol)
    plt.plot(z, nz, color='red')
    plt.plot(ztarget, nztarget, color='blue')
    np.savetxt('nz_lrg_mockComparing.txt', np.array([ztarget, nztarget]).T)
    plt.title(type_)
    plt.show()
    print("\nLegend:")
    print("  red  Reference")
    print("  blue Comparison")

if type_ == 'QSO':

    z, nz = np.loadtxt('/global/homes/a/acarnero/nz_qso_v5.txt', unpack = True)
    #z, nz = np.loadtxt('/pscratch/sd/a/acarnero/codes/desi-cutsky-mock/nz_files/DA2/nz_qso_da2_mocks.txt', unpack = True)
    ztarget, nztarget = mknz(data, area, zmax=4, zcol=args.zrsdcol)
    plt.plot(z, nz, color='red')
    plt.plot(ztarget, nztarget, color='blue')
    np.savetxt('nz_qso_mockComparing.txt', np.array([ztarget, nztarget]).T)
    plt.title(type_)
    plt.show()
    print("\nLegend:")
    print("  red  Reference")
    print("  blue Comparison")

def return_north(ra, dec):
    '''
    given a table that already includes RA,DEC, add PHOTSYS column denoting whether
    the data is in the DECaLS ('S') or BASS/MzLS ('N') photometric region
    '''
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    c = SkyCoord(ra* u.deg, dec* u.deg,frame='icrs')
    gc = c.transform_to('galactic')
    sel_ngc = gc.b > 0

    seln = dec > 32.375

    sel = seln&sel_ngc
    return sel


if type_ == 'ELG':
    divide = True
    if divide:
        masknorth = return_north(data['RA'], data['DEC'])
        datanorth = data[masknorth]
        datasouth = data[~masknorth]

        zN, nzN = np.loadtxt('/pscratch/sd/a/acarnero/codes/desi-cutsky-mock/nz_files/DA2/nz_elg_N_da2_mocks.txt', unpack = True, usecols=([0,1]))
        zS, nzS = np.loadtxt('/pscratch/sd/a/acarnero/codes/desi-cutsky-mock/nz_files/DA2/nz_elg_S_da2_mocks.txt', unpack = True, usecols=([0,1]))

#        maskS = (zS>=0.4)
#        zS = zS[maskS]
#        nzS = nzS[maskS]
        
#        maskN = (zN>=0.4)
#        zN = zN[maskN]
#        nzN = nzN[maskN]

        random = Table.read(ranf)
        masknorth = return_north(random['RA'], random['DEC'])
        randomnorth = random[masknorth]
        randomsouth = random[~masknorth]

        ztargetN, nztargetN = mknz(datanorth, len(randomnorth)/2500., zmax=2., zmin=0.,zcol=args.zrsdcol)
        ztargetS, nztargetS = mknz(datasouth, len(randomsouth)/2500., zmax=2., zmin=0.,zcol=args.zrsdcol)

        plt.plot(zS, nzS, color='red')
        plt.plot(ztargetS, nztargetS, color='blue')


        #plt.plot(zS,nzS/nztargetS)
        np.savetxt('nz_elgS_mockComparing.txt', np.array([ztargetS, nztargetS]).T)
        plt.title(type_ + ' SOUTH')
        plt.show()
        print("\nLegend SOUTH:")
        print("  red  Reference")
        print("  blue Comparison")
        plt.clear_figure()
        plt.plot(zN, nzN, color='red')
        plt.plot(ztargetN, nztargetN, color='blue')
        np.savetxt('nz_elgN_mockComparing.txt', np.array([ztargetS, nztargetS]).T)
        plt.title(type_ + ' NORTH')
        plt.show()
        print("\nLegend NORTH:")
        print("  red  Reference")
        print("  blue Comparison")

    else:
        print('')
        zN, nzN = np.loadtxt('/pscratch/sd/a/acarnero/codes/desi-cutsky-mock/nz_files/HighFidelity/nz_elg_highfidel_max.txt', unpack = True, usecols=([0,1]))
        ztargetS, nztargetS = mknz(data, len(ranf)/2500., zmax=2., zmin=0., zcol=args.zrsdcol)#, zmin=0.4)
        plt.plot(zN, nzN, color='red')
        plt.plot(ztargetS, nztargetS, color='blue')
        plt.title(type_ + ' COMBINED')
        plt.show()
        print("\nLegend:")
        print("  red  Reference")
        print("  blue Comparison")
