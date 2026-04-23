import numpy as np
from astropy.table import Table,join,vstack
import os
import sys
from matplotlib import pyplot as plt
#sys.path.append(os.environ['HOME']+'/LSS/py')#wherever you git cloned the LSS repo
import LSS.common_tools as common
import h5py
import hdf5plugin
import astropy.io.fits as pf
from calibrate_nz_prep import calibrate_nz_bgs
import LSS.common_tools as cm
from astropy.io import fits # Access to FITS (Flexible Image Transport System) files.

from desitarget.targetmask import obsconditions

def mknz(df, lenran,fout, bs=0.01, zmin=0.02, zmax=0.6, randens=2500.):

    # cd = distance(om,1-om)
    headers = []


        # should have originally had 2500/deg2 density, so can convert to area

    area = lenran/randens

    print('area is '+str(area))
    headers.append('#area is '+str(area)+'square degrees')

    nbin = int((zmax-zmin)*(1+bs/10)/bs)
    #wts = None # no weights, if relevant at all
    zhist = np.histogram(df['Z'], bins=nbin, range=(zmin, zmax))#, weights=wts)
    headers.append('#zmid zlow zhigh n(z) Nbin Vol_bin')
    zl = zhist[1][:-1]
    zh = zhist[1][1:]
    zm = (zh+zl)/2.
    voli = area/(360.*360./np.pi)*4.*np.pi/3. * np.diff(cm.dis_dc(zhist[1])**3.)
    nbarz = zhist[0]/voli
    np.savetxt(fout, np.column_stack([zm, zl, zh, nbarz, zhist[0], voli]), fmt='%.18g', header='\n'.join(headers), comments='')
    return True


tile = 'BRIGHT'

wd=Table.read('/pscratch/sd/d/desica/DA2/mocks/white_dwarfs_bgs/wd_sample.fits')
wd.remove_columns(['TARGETID', 'REF_EPOCH', 'PARALLAX', 'PMRA','PMDEC','NUMOBS', 'ZTILEID','Z_QN','IS_QSO_QN','DELTACHI2','TARGET_STATE','TIMESTAMP','VERSION'])
wd.rename_column('Z', 'RSDZ')
swd = len(wd)
wd['NZ'] = [0.0001]*swd
wd['RA'] = wd['RA'].astype(np.float64)
wd['DEC'] = wd['DEC'].astype(np.float64)
wd['TRUEZ'] = np.ones(swd)
wd['R_MAG_ABS'] = np.ones(swd)*-21
wd['BRICKNAME'] = np.full(swd, '000p0000')    #- required !?!
wd['WEIGHT'] = 1.
wd['ZWARN'] = np.zeros(swd, dtype='i8')+int(0)



files = ['holi_BGS_v4.82_GCcomb_clustering.dat.h5', 'holi_BGS-NONKP_v4.82_GCcomb_clustering.dat.h5']
path = '/global/cfs/cdirs/desi/mocks/cai/holi/webjax_v4.82'


seeds_good = []
seeds_missing = []
for i in range(1000):
    if i<50:
        continue
    seed = str(i).zfill(4)
    dest = os.path.join(path,'seed%s' % seed)
    if os.path.isfile(os.path.join(dest,files[0])) and os.path.isfile(os.path.join(dest,files[1])):
        seeds_good.append(i)



path_to_save_nz = '/pscratch/sd/d/desica/holi_bgs_nz_temp'

for i in seeds_good:
    seed = str(i).zfill(4)

    print(seed)
    dest = os.path.join(path,'seed%s' % seed)
    datos_nonkp = os.path.join(os.path.join(dest, files[1]))   #'/global/cfs/cdirs/desi/mocks/cai/holi/webjax_v4.80/seed0001/holi_BGS-NONKP_v4.80_GCcomb_clustering.dat.h5'
    mock = Table(h5py.File(datos_nonkp, 'r+'))
    mknz(mock, 35957774, os.path.join(path_to_save_nz, 'nz_holi_bgs_nonkp_extended_%s.txt' % seed), zmax=1.6)
    z,nz = np.loadtxt(os.path.join(path_to_save_nz, 'nz_holi_bgs_nonkp_extended_%s.txt' % seed), unpack=True, usecols=(0,3))
    the_bgs_any = calibrate_nz_bgs(mock, redshift_column = 'Z', tracer_type='BGS', survey='DA2', save_mock_nz = 'n', n_mean=[z,nz], nzfile='/global/cfs/cdirs/desi/survey/catalogs/DA2/LSS/loa-v1/LSScats/v2/BGS_ANY_full_HPmapcut_nz.txt')
    the_bgs_any_new = the_bgs_any[the_bgs_any['STATUS']==1]
    print(len(the_bgs_any_new), 'should be around 19623289')
    val=19623289-len(the_bgs_any_new)
    print(val)
    mknz(the_bgs_any_new, 35957774, os.path.join(path_to_save_nz, 'nz_holi_bgs_any_v4_seed%s.txt' % seed), zmax=1.6)
    the_bgs_any_new['DESI_TARGET'] = 2**60
    the_bgs_any_new['R_MAG_ABS'] = -21.
    z,nz = np.loadtxt(os.path.join(path_to_save_nz, 'nz_holi_bgs_any_v4_seed%s.txt' % seed), unpack=True,usecols=(0,3))
    the_bgs_bright = calibrate_nz_bgs(the_bgs_any_new, redshift_column = 'Z', tracer_type='BGS', survey='DA2', save_mock_nz = 'n', n_mean=[z,nz], nzfile='/global/cfs/cdirs/desi/survey/catalogs/DA2/LSS/loa-v1/LSScats/v2/BGS_BRIGHT_full_HPmapcut_nz.txt')

    the_bgs_bright_new = the_bgs_bright[the_bgs_bright['STATUS']==1]
    the_bgs_rest = the_bgs_bright[the_bgs_bright['STATUS']==0]
#print(len(the_bgs_bright_new), len(the_bgs_rest))
    mknz(the_bgs_bright_new, 35957774, os.path.join(path_to_save_nz, 'nz_holi_bgs_bright_v4_seed%s.txt' % seed), zmax=1.6)
    the_bgs_bright_new['BGS_TARGET'] = 2**1

    dowisebgs = 0.00585
    ran_rest = np.random.uniform(size = len(the_bgs_rest))
    wisebgs_mask = (ran_rest <= dowisebgs)

    the_bgs_faint_final = the_bgs_rest[~wisebgs_mask]
    the_bgs_wise_final = the_bgs_rest[wisebgs_mask]
    the_bgs_wise_final['BGS_TARGET'] = 2**2
    the_bgs_faint_final['BGS_TARGET'] = 2**0
    print(len(the_bgs_wise_final), 'size BGS WISE should be around 43742')
    val=43742-len(the_bgs_wise_final)
    print(val)
    the_bgs_faint_final['PRIORITY_INIT'] = 2000
    the_bgs_faint_final['PRIORITY'] = 2000
    the_bgs_wise_final['PRIORITY_INIT'] = 2000
    the_bgs_wise_final['PRIORITY'] = 2000

    PromoteFracBGSFaint=0.2
    ran_hip = np.random.uniform(size = len(the_bgs_faint_final))
    faint_hip_mask = (ran_hip <= PromoteFracBGSFaint)

    the_bgs_faint_final['BGS_TARGET'][faint_hip_mask] += 2**3   # for high-priority BGS faint
    the_bgs_faint_final['PRIORITY'][faint_hip_mask] = 2100   # for high-priority BGS faint
    the_bgs_faint_final['PRIORITY_INIT'][faint_hip_mask] = 2100

    targets_faint = vstack([the_bgs_faint_final,the_bgs_wise_final])
    targets_faint['R_MAG_ABS'] = -19.


    maska = (the_bgs_rest['Z'] > 0.5) & (the_bgs_rest['Z']<0.61)
    some = the_bgs_rest[maska]

    random_indices = np.random.choice(len(some), size=100000, replace=False)
    t_sampled = some[random_indices]

    datos_nonkp = os.path.join(dest, files[0])
    mock = Table(h5py.File(datos_nonkp, 'r+'))
    mocka = vstack([mock, t_sampled])
    mknz(mocka, 35957774, os.path.join(path_to_save_nz, 'nz_holi_bgs_cosmosample_seed%s.txt' % seed), zmax=0.61)
    zsim, nzsim = np.loadtxt(os.path.join(path_to_save_nz, 'nz_holi_bgs_cosmosample_seed%s.txt' % seed), unpack=True, usecols=(0,3))

    the_bgs_cosmo = calibrate_nz_bgs(mocka, redshift_column = 'Z', tracer_type='BGS', survey='DA2', save_mock_nz = 'n', n_mean=[zsim,nzsim], nzfile='/global/cfs/cdirs/desi/survey/catalogs/DA2/LSS/loa-v1/LSScats/v2/nonKP/BGS_BRIGHT-21.35_NGC_nz.txt')
    the_real_bgs_cosmo = the_bgs_cosmo[the_bgs_cosmo['STATUS']==1]
    mknz(the_real_bgs_cosmo, 35957774, os.path.join(path_to_save_nz, 'nz_holi_bgs_bright-21.35_v4_seed%s.txt' % seed), zmax=0.61)
    zsim,nzsim=np.loadtxt(os.path.join(path_to_save_nz, 'nz_holi_bgs_bright_v4_seed%s.txt' % seed), unpack=True,usecols=(0,3))
    temp_holi_bright = calibrate_nz_bgs(the_bgs_bright_new, redshift_column = 'Z', tracer_type='BGS', survey='DA2', save_mock_nz = 'n', n_mean=[zsim,nzsim], nzfile=os.path.join(path_to_save_nz, 'nz_holi_bgs_bright-21.35_v4_seed%s.txt' % seed))
    final_bright = temp_holi_bright[temp_holi_bright['STATUS'] == 0]
    the_real_bgs_cosmo['R_MAG_ABS'] = -22.
    final_bright['R_MAG_ABS'] = -21.




    datat = []
    the_real_bgs_cosmo['BGS_TARGET'] = 2**1
    datat.append(Table(final_bright))
    datat.append(Table(the_real_bgs_cosmo))

    targets_bright = vstack(datat)
    targets_bright['PRIORITY_INIT'] = 2100
    targets_bright['PRIORITY'] = 2100

    all_targets = vstack([targets_faint,targets_bright])
    all_targets['DESI_TARGET'] = 2**60
    all_targets['NUMOBS_MORE'] = 2
    all_targets['NUMOBS_INIT'] = 2
    all_targets['RA'] = all_targets['RA'].astype(np.float64)
    all_targets['DEC'] = all_targets['DEC'].astype(np.float64)

    all_targets.rename_column('Z_COSMO', 'TRUEZ')
    all_targets.rename_column('Z', 'RSDZ')
        

    n = len(all_targets)

    all_targets['MWS_TARGET'] = np.zeros(n, dtype='i8')
    all_targets['SUBPRIORITY'] = np.random.uniform(0, 1, n)
    all_targets['BRICKNAME'] = np.full(n, '000p0000')    #- required !?!
    all_targets['SCND_TARGET'] = np.zeros(n, dtype='i8')+int(0)
    all_targets['ZWARN'] = np.zeros(n, dtype='i8')+int(0)
    all_targets['WEIGHT'] = 1.
    tile = 'BRIGHT'
    all_targets['OBSCONDITIONS'] = obsconditions.mask(tile)
    all_targets.remove_column('STATUS')

    if len(all_targets.columns) != len(wd.columns):
        print('something went wrong with columns')
    
    cata = vstack([all_targets, wd])
    cata['TARGETID'] = (np.arange(1,len(cata)+1)+1e8).astype(int)
    out_file_name = os.path.join('/pscratch/sd/d/desica/DA2/mocks/holi_bgs', 'forFA%d.fits' %i)
    cata.write(out_file_name, overwrite=True)


    fits.setval(out_file_name, 'EXTNAME', value='TARGETS', ext=1)
    fits.setval(out_file_name, 'OBSCON', value=tile, ext=1)
    print('end seed',seed)



