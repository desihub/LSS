from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import fitsio
import gc
import desimodel.io
import desitarget.mtl
import desisim.quickcat
from astropy.io import fits
from astropy.table import Table, Column, vstack
import json
import shutil
import healpy
from desitarget.targetmask import desi_mask, obsconditions
from collections import Counter
import subprocess


def ra_dec_subset(data, ra_min=130, ra_max=190, dec_min=-5, dec_max=15):
    subset_ii = (data['RA']>ra_min) & (data['RA']<ra_max)
    subset_ii &= (data['DEC']>dec_min) & (data['DEC']<dec_max)
    return subset_ii
       


def accurate_assign_lya_qso(initial_mtl_file, pixweight_file):
    print("Finding targets that will be lya in the truth file")
    print('started reading targets')
    targets = Table.read(initial_mtl_file)
    print('finished reading targets')

    pixweight, header = fits.getdata(pixweight_file, 'PIXWEIGHTS', header=True)
    hpxnside = header['HPXNSIDE']

    theta_w, phi_w = healpy.pix2ang(hpxnside, pixweight['HPXPIXEL'], nest=True)

    hpxnside_sample = 64 # pixel area on which we will sample the lyaqso
    npix_sample = healpy.nside2npix(hpxnside_sample)
    pixnumber_sample = healpy.ang2pix(hpxnside_sample, theta_w, phi_w, nest=True)

    subpixels = {} # store the pixels at the resolution in the input target catalog that are included within the pixels at the new resolution.
    for i in range(npix_sample):
        ii_sample = pixnumber_sample==i
        subpixels[i] =  healpy.ang2pix(hpxnside, theta_w[ii_sample], phi_w[ii_sample], nest=True)
    
    
    # redefine the covered area for the new pixels from the higher resolution FRACAREA
    covered_area = np.ones(npix_sample)
    for i in range(npix_sample):
        sum_weight = np.sum(pixweight['FRACAREA'][subpixels[i]])
        if sum_weight>0.0:
            covered_area[i] = np.sum(pixweight['FRACAREA'][subpixels[i]]**2)/np.sum(pixweight['FRACAREA'][subpixels[i]])
        else:
            covered_area[i] = 0.0
    print('finished computing covered area')
    theta_s, phi_s = healpy.pix2ang(hpxnside_sample, np.arange(npix_sample), nest=True)

    pixelarea_sample = healpy.pixelfunc.nside2pixarea(hpxnside_sample, degrees=True)
    n_lya_qso_in_pixel = np.int_(covered_area * 50 * pixelarea_sample)

    # compute angular coordinates from the targets
    targets_phi = np.deg2rad(targets['RA'])
    targets_theta = np.deg2rad(90.0-targets['DEC'])

    # find the pixnumber to which the target belongs (in the hpxnside_sample resolution)
    pixnumber_targets = healpy.ang2pix(hpxnside_sample, targets_theta, targets_phi, nest=True)
    print('finished computed pixnumber for all targets')
    
    # what target are QSOs?
    is_qso = (targets['DESI_TARGET'] & desi_mask.QSO)!=0

    # list of unique pixels covered by the targets
    pixnumber_target_list = list(set(pixnumber_targets)) # list of pixelsIDs covered by the targets in the new resolution

    n_qso_per_pixel_targets = np.zeros(len(pixnumber_target_list))
    n_pixnumber_target_list = len(pixnumber_target_list)
    print('started counting pixnumber target list', n_pixnumber_target_list)
    for i in range(n_pixnumber_target_list):
        ii_targets = is_qso & (pixnumber_targets==pixnumber_target_list[i])
        n_qso_per_pixel_targets[i] = np.count_nonzero(ii_targets)
    print('finished counting pixnumber target list')
    
    n_lya_desired_pixel_targets = np.random.poisson(n_lya_qso_in_pixel[pixnumber_target_list])
    print('finished random poisson')
    
    # Generate the boolean array to determine whether a target is a lyaqso or not
    n_targets = len(targets)
    is_lya_qso = np.repeat(False, n_targets)
    target_ids = np.arange(n_targets)
    print('started looping over pixnumber_target_list')
    for i in range(len(pixnumber_target_list)):
        ii_targets = is_qso & (pixnumber_targets==pixnumber_target_list[i])
        n_qso_in_pixel = np.count_nonzero(ii_targets)
        n_lya_desired = n_lya_desired_pixel_targets[i]
        if n_lya_desired >= n_qso_in_pixel:
            is_lya_qso[ii_targets] = True
        else:
            #print(len(target_ids[ii_targets]), n_lya_desired)
            ii_lya_qso = np.random.choice(target_ids[ii_targets], n_lya_desired, replace=False)
            is_lya_qso[ii_lya_qso] = True
    return is_lya_qso

def make_global_DR8_sky(output_path="./"):
    # Create output directory
    os.makedirs(output_path, exist_ok=True)
    global_DR8_sky_file = os.path.join(output_path, "global_DR8_sky.fits")

    if os.path.exists(global_DR8_sky_file):
        print("File {} already exist".format(global_DR8_sky_file))
        return global_DR8_sky_file
    
    print('Preparing file {}'.format(global_DR8_sky_file))

    columns = ['TARGETID', 'DESI_TARGET', 'MWS_TARGET', 'BGS_TARGET', 'SUBPRIORITY', 'NUMOBS_INIT', 'PRIORITY_INIT', 'RA', 'DEC', 'HPXPIXEL', 'BRICKNAME', 'OBSCONDITIONS']
    
    # List all the fits files to read
    path_to_targets = '/global/cfs/projectdirs/desi/target/catalogs/dr8/0.39.0/skies/'
    target_files = glob.glob(os.path.join(path_to_targets, "skies-*.fits"))
    print('sky files to read:', len(target_files))
    target_files.sort()
    
    # Read the first file, only the columns that are useful for MTL
    data = fitsio.FITS(target_files[0], 'r')
    target_data = data[1].read(columns=columns)
    data.close()
    
    # Read all the other files
    for i, i_name in enumerate(target_files[1:]): 
        data = fitsio.FITS(i_name, 'r')
        tmp_data = data[1].read(columns=columns)
        target_data = np.hstack((target_data, tmp_data))
        data.close()
        print('reading file', i, len(target_files), len(tmp_data))

    #target_data = Table(target_data)

    print('Started writing file {}'.format(global_DR8_sky_file))
    #target_data.write(global_DR8_sky_file, overwrite=True)
    outfd = fitsio.FITS(global_DR8_sky_file, "rw")
    outfd.write(None, header=None, extname="PRIMARY")
    outfd.write(target_data, header=None, extname="TARGETS")
    outfd.close()

    print('Finished writing file {}'.format(global_DR8_sky_file))
    del target_data
    return global_DR8_sky_file

def random_assign_lya_qso(targets, fraction=0.25):
    """
    Assign a fraction of all QSOs to by a Lya QSO.
    """
        
    n_targets = len(targets)

    # what target are QSOs?
    is_qso = (targets['DESI_TARGET'] & desi_mask.QSO)!=0
    n_qso = np.count_nonzero(is_qso)
    
    is_lya_qso = np.repeat(False, n_targets)

    target_ids = np.arange(n_targets)

    n_lya_qso = int(fraction * n_qso)
    ii_lya_qso = np.random.choice(target_ids[is_qso], n_lya_qso, replace=False)
    is_lya_qso[ii_lya_qso] = True
    print('Number of total QSOs: {}'.format(n_qso))
    print('Numer of lya QSOs: {}'.format(n_lya_qso))
    return is_lya_qso

def make_global_DR8_mtl(output_path='./', program='dark'):
    os.makedirs(output_path, exist_ok=True)
    
    global_DR8_mtl_file = os.path.join(output_path, 'global_DR8_mtl_{}.fits'.format(program))
    if os.path.exists(global_DR8_mtl_file):
        print("File {} already exists".format(global_DR8_mtl_file))
        return global_DR8_mtl_file
    
    print('Preparing file {}'.format(global_DR8_mtl_file))
    # List all the fits files to read
    path_to_targets = '/global/cfs/projectdirs/desi/target/catalogs/dr8/0.39.0/targets/main/resolve/'+program+'/'
    target_files = glob.glob(os.path.join(path_to_targets, "targets*fits"))
    print('target files to read:', len(target_files))
    target_files.sort()
    
    columns = ['TARGETID', 'DESI_TARGET', 'MWS_TARGET', 'BGS_TARGET', 'SUBPRIORITY', 'NUMOBS_INIT', 'PRIORITY_INIT', 'RA', 'DEC', 'HPXPIXEL', 'BRICKNAME', 'FLUX_R', 'MW_TRANSMISSION_R']

    data = fitsio.FITS(target_files[0], 'r')
    target_data = data[1].read(columns=columns)
    data.close()
    for i, i_name in enumerate(target_files[1:]):
        data = fitsio.FITS(i_name, 'r')
        tmp_data = data[1].read(columns=columns)
        target_data = np.hstack((target_data, tmp_data))
        data.close()
        print('reading file', i, len(target_files), len(tmp_data))
    
    if program=='dark':
        full_mtl = desitarget.mtl.make_mtl(target_data, 'DARK|GRAY')
    if program=='bright':
        full_mtl = desitarget.mtl.make_mtl(target_data, 'BRIGHT')

    print('Started writing file {}'.format(global_DR8_mtl_file))
    full_mtl.write(global_DR8_mtl_file, overwrite=True)
    print('Finished writing file {}'.format(global_DR8_mtl_file))

    del full_mtl
    return global_DR8_mtl_file

def make_global_DR8_truth(global_DR8_mtl_file, output_path='./', program='dark'):
    import desitarget.mock.mockmaker as mb
    from desitarget.targetmask import desi_mask, bgs_mask, mws_mask
    
    os.makedirs(output_path, exist_ok=True)
    global_DR8_truth_file = os.path.join(output_path, 'global_DR8_truth_{}.fits'.format(program))
    if os.path.exists(global_DR8_truth_file):
        print("File {} already exists".format(global_DR8_truth_file))
        return global_DR8_truth_file
    
    print('Started reading file {}'.format(global_DR8_mtl_file))
    targets = Table.read(global_DR8_mtl_file)
    print('Finished reading file {}'.format(global_DR8_mtl_file))

    # Find what targets will be associated to lya targets
    is_lya_qso = random_assign_lya_qso(targets)
    
    # Initialized truth Table
    colnames = list(targets.dtype.names)
    print(colnames)
    nobj = len(targets)
    truth = mb.empty_truth_table(nobj=nobj)[0]
    print(truth.keys())

    for k in colnames:
        if k in truth.keys():
            print(k)
            truth[k][:] = targets[k][:]

    nothing = '          '
    truth['TEMPLATESUBTYPE'] = np.repeat(nothing, nobj)

    masks = ['MWS_ANY', 'BGS_ANY', 'STD_FAINT', 'STD_BRIGHT','ELG', 'LRG', 'QSO', ]
    dict_truespectype = {'BGS_ANY':'GALAXY', 'ELG':'GALAXY', 'LRG':'GALAXY', 'QSO':'QSO', 
                    'MWS_ANY':'STAR', 'STD_FAINT':'STAR', 'STD_BRIGHT':'STAR'}
    dict_truetemplatetype = {'BGS_ANY':'BGS', 'ELG':'ELG', 'LRG':'LRG', 'QSO':'QSO', 
                        'MWS_ANY':'STAR', 'STD_FAINT':'STAR', 'STD_BRIGHT':'STAR'}
    dict_truez = {'BGS_ANY':0.2, 'ELG':1.5, 'LRG':0.7, 'QSO':2.0, 
                        'MWS_ANY':0.0, 'STD_FAINT':0.0, 'STD_BRIGHT':0.0}

    for m in masks:
        istype = (targets['DESI_TARGET'] & desi_mask.mask(m))!=0
        print(m, np.count_nonzero(istype))
        truth['TRUESPECTYPE'][istype] = np.repeat(dict_truespectype[m], np.count_nonzero(istype))
        truth['TEMPLATETYPE'][istype] = np.repeat(dict_truetemplatetype[m], np.count_nonzero(istype))
        truth['MOCKID'][istype] = targets['TARGETID'][istype]
        truth['TRUEZ'][istype] = dict_truez[m]
        
    truth['TRUEZ'][is_lya_qso] = 3.0
    truth['RA'] = targets['RA']
    truth['DEC'] = targets['DEC']

    # Check that all targets have been assigned to a class
    iii = truth['MOCKID']==0
    assert np.count_nonzero(iii)==0
    
    del targets
    
    print('Started writing to file {}'.format(global_DR8_truth_file))
    truth.write(global_DR8_truth_file, overwrite=True)
    print('Finished writing to file {}'.format(global_DR8_truth_file))
    
    del truth
    return global_DR8_truth_file
    
def prepare_tile_batches(surveysim_file, output_path='./', program='dark', start_day=0, end_day=365, batch_cadence=7, 
                        select_subset_sky=False, ra_min=130, ra_max=190, dec_min=-5, dec_max=15,use_last_date=False):
    
    os.makedirs(output_path, exist_ok=True)

    all_exposures = Table.read(surveysim_file, hdu=1)
    all_tiledata = Table.read(surveysim_file, hdu=2)
    all_tiles = desimodel.io.load_tiles()

    all_exposures['MJD_OFFSET'] = all_exposures['MJD'] - all_exposures['MJD'].min()

    if program=='dark':
        ii = (all_tiles['PROGRAM']=='DARK') | (all_tiles['PROGRAM']=='GRAY')
        all_tiles = all_tiles[ii]
    elif program=='bright':
        ii = (all_tiles['PROGRAM']=='BRIGHT')
        all_tiles = all_tiles[ii]

    if select_subset_sky:
        ii = ra_dec_subset(all_tiles, ra_min=ra_min, ra_max=ra_max, dec_min=dec_min, dec_max=dec_max)
        tiles = all_tiles[ii]
    else:
        tiles = all_tiles.copy()
    tiles = Table(tiles)
    
    ii = np.isin(all_exposures['TILEID'], tiles['TILEID'])    
    exposures = all_exposures[ii]
    a, b = np.unique(exposures['TILEID'], return_index=True)
    unique_tiles = exposures['TILEID'][np.sort(b)]
    unique_dates = exposures['MJD_OFFSET'][np.sort(b)]
    
    if use_last_date:
        for tile in unique_tiles:
            wt = all_exposures['TILEID'] == tile
            md = np.max(all_exposures['MJD_OFFSET'][wt])
            wt = unique_tiles == tile
            od = unique_dates[wt]
            unique_dates[wt] = md
            print(tile,md,od)
        

    i_day  = start_day 
    batch_id = int(start_day/batch_cadence)
    if start_day/batch_cadence - batch_id != 0:
        print('mismatch between starting day and initial batch id')
        print(start_day/batch_cadence,batch_id)
    while (i_day + batch_cadence) < end_day:
        min_day  = i_day 
        max_day  = min_day + batch_cadence
        ii = (unique_dates>=min_day) & (unique_dates<max_day)
        tiles_in_batch = unique_tiles[ii]
        
        # count the batch only if it has tiles in it
        if len(tiles_in_batch):
            print('batch_{:04d} {}'.format(batch_id, len(tiles_in_batch)), max_day)
            jj = np.isin(tiles['TILEID'], tiles_in_batch)
            batch_filename = os.path.join(output_path, 'batch_{:04d}_{}.fits'.format(batch_id, program))
            tiles[jj].write(batch_filename, overwrite=True)
            batch_id += 1
            
        i_day += batch_cadence

    return batch_id
 

def make_patch_file(data_filename, ra_min=130, ra_max=190, dec_min=-5, dec_max=15):
    patch_filename = data_filename.replace("global", "patch")
    if os.path.exists(patch_filename):
        print("File {} already exists".format(patch_filename))
        return patch_filename
    print('Creating file {}'.format(patch_filename))
    
    print('started reading targets')
    data = Table.read(data_filename)
    print('finished reading targets')
    
    print('started patch selection')
    ii = ra_dec_subset(data, ra_min=ra_min, ra_max=ra_max, dec_min=dec_min, dec_max=dec_max)
    print('finished patch selection')
    
    print('Started writing to file {}'.format(patch_filename))
    data[ii].write(patch_filename, overwrite=True)
    print('Finished writing to file {}'.format(patch_filename))
    
    del data
    
    return patch_filename 

    
def run_strategy(initial_mtl_file, truth_file, sky_file, output_path="./", batch_path="./", program='dark',sbatch=0,mxbatch=100):
    os.makedirs(output_path, exist_ok=True)
    targets_path='{}/targets'.format(output_path)
    zcat_path = '{}/zcat'.format(output_path)
    os.makedirs(targets_path, exist_ok=True)
    os.makedirs(zcat_path, exist_ok=True)
    os.makedirs('{}/fiberassign'.format(output_path), exist_ok=True)

    batch_files = glob.glob(batch_path+"/batch_*_"+program+".fits")
    batch_files.sort()
    
    # Read targets and truth
    print('reading truth file')
    truth = Table.read(truth_file)
    print('done readying truth file')
    #obsconditions
    obsconditions = None
    if program=='dark':
        obsconditions = 'DARK|GRAY'
    if program=='bright':
        obsconditions = 'BRIGHT'
    
    n_batch = len(batch_files)
    if mxbatch < n_batch:
        n_batch = mxbatch
    for i_batch in range(sbatch,n_batch):
        print()
        print("Batch {}".format(i_batch))
        fiberassign_path = '{}/fiberassign/{:04d}'.format(output_path, i_batch)
        os.makedirs(fiberassign_path, exist_ok=True)

        footprint = batch_files[i_batch]
        mtl_filename = os.path.join(targets_path, '{:04d}_mtl.fits'.format(i_batch))
        new_mtl_filename = os.path.join(targets_path, '{:04d}_mtl.fits'.format(i_batch+1))

        zcat_filename = os.path.join(zcat_path, '{:04d}_zcat.fits'.format(i_batch))
        old_zcat_filename = os.path.join(zcat_path, '{:04d}_zcat.fits'.format(i_batch-1))
        
        
        if i_batch == 0:
            shutil.copyfile(initial_mtl_file, mtl_filename)
        print(footprint)
        
        fba_run = 'fba_run --targets {} --sky {} --footprint {}  --dir {} --rundate 2020-01-01T00:00:00 --overwrite'.format(
            mtl_filename, sky_file, footprint, fiberassign_path)
        print(fba_run)
        os.system(fba_run)
    
        # Gather fiberassign files
        fba_files = np.sort(glob.glob(os.path.join(fiberassign_path,"fba-*.fits")))
        
        # read the current mtl file
        targets = Table.read(mtl_filename)

        # Compute zcat
        if i_batch==0:
            zcat = desisim.quickcat.quickcat(fba_files, targets, truth, fassignhdu='FASSIGN', perfect=True)
        else:
            old_zcat = Table.read(old_zcat_filename)
            zcat = desisim.quickcat.quickcat(fba_files, targets, truth, fassignhdu='FASSIGN', zcat=old_zcat, perfect=True)  
            del old_zcat      
    
        zcat.write(zcat_filename, overwrite=True)
        mtl = desitarget.mtl.make_mtl(targets, obsconditions, zcat=zcat)
        
        del targets
        del zcat
        #gc.collect()
        print('writing mtl ')

        mtl.write(new_mtl_filename, overwrite=True)
        print('wrote mtl')
        del mtl
        gc.collect()
        
    
    return True    
        