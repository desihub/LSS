from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import fitsio
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



def ra_dec_subset(data, ra_min=130, ra_max=180, dec_min=-10, dec_max=40):
    subset_ii = (data['RA']>ra_min) & (data['RA']<ra_max)
    subset_ii &= (data['DEC']>dec_min) & (data['DEC']<dec_max)
    return subset_ii

def assign_lya_qso(initial_mtl_file, pixweight_file):
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

def write_initial_mtl_file(initial_mtl_file, ra_min=130, ra_max=180, dec_min=-10, dec_max=40):
    
    target_files = glob.glob(os.path.join(path_to_targets, "targets*fits"))
    print('target files to read:', len(target_files))
    target_files.sort()
    
    data = fitsio.FITS(target_files[0], 'r')
    target_data = data[1].read(columns=['TARGETID', 'DESI_TARGET', 'MWS_TARGET', 'BGS_TARGET', 'SUBPRIORITY', 'NUMOBS_INIT', 'PRIORITY_INIT', 'RA', 'DEC', 'HPXPIXEL', 'BRICKNAME'])
    data.close()
    for i, i_name in enumerate(target_files[1:]):
        data = fitsio.FITS(i_name, 'r')
        tmp_data = data[1].read(columns=['TARGETID', 'DESI_TARGET', 'MWS_TARGET', 'BGS_TARGET', 'SUBPRIORITY', 'NUMOBS_INIT', 'PRIORITY_INIT', 'RA', 'DEC', 'HPXPIXEL', 'BRICKNAME'])
        target_data = np.hstack((target_data, tmp_data))
        data.close()
        print('reading file', i, len(target_files), len(tmp_data))
    full_mtl = desitarget.mtl.make_mtl(target_data, 'DARK|GRAY')

    ii_mtl_dark = (full_mtl['OBSCONDITIONS'] & obsconditions.DARK)!=0
    ii_mtl_gray = (full_mtl['OBSCONDITIONS'] & obsconditions.GRAY)!=0
    ii_north = (full_mtl['RA']>85) & (full_mtl['RA']<300) & (full_mtl['DEC']>-15)

    print("Writing nothern cap")
    mtl_file = "targets/dr8_mtl_dark_gray_NGC.fits"
    full_mtl[(ii_mtl_dark | ii_mtl_gray) & ii_north].write(outdir+mtl_file, overwrite=True)
    
    print("Writing subset in the northern cap")
    mtl_data = Table.read(mtl_file)
    subset_ii = ra_dec_subset(mtl_data,ra_min,ra_max,dec_min,dec_max)
    mtl_data[subset_ii].write(outdir+initial_mtl_file, overwrite=True)

def write_initial_mtl_file_AJR(initial_mtl_file, ra_min=130, ra_max=180, dec_min=-10, dec_max=40,outdir=''):
    
    data = fitsio.FITS(full_target_data, 'r')
    tmp_data = data[1].read(columns=['TARGETID', 'DESI_TARGET', 'MWS_TARGET', 'BGS_TARGET', 'SUBPRIORITY', 'NUMOBS_INIT', 'PRIORITY_INIT', 'RA', 'DEC', 'HPXPIXEL', 'BRICKNAME'])
    full_mtl = desitarget.mtl.make_mtl(tmp_data, 'DARK|GRAY')

    ii_mtl_dark = (full_mtl['OBSCONDITIONS'] & obsconditions.DARK)!=0
    ii_mtl_gray = (full_mtl['OBSCONDITIONS'] & obsconditions.GRAY)!=0
    ii_north = (full_mtl['RA']>85) & (full_mtl['RA']<300) & (full_mtl['DEC']>-15)

    print("Writing nothern cap")
    mtl_file = "targets/dr8_mtl_dark_gray_NGC.fits"
    full_mtl[(ii_mtl_dark | ii_mtl_gray) & ii_north].write(outdir+mtl_file, overwrite=True)
    
    print("Writing subset in the northern cap")
    mtl_data = Table.read(outdir+mtl_file)
    
    subset_ii = ra_dec_subset(mtl_data,ra_min,ra_max,dec_min,dec_max)
    mtl_data[subset_ii].write(outdir+initial_mtl_file, overwrite=True)
    
    print('NOTHING done yet in SGC')

def write_initial_mtl_files(cap='NGC',dr='dr8', ra_min=130, ra_max=180, dec_min=-10, dec_max=40,outdir='',full_target_data='',obscon='DARK|GRAY',sky_data_file=''):
    
    data = fitsio.FITS(full_target_data, 'r')
    tmp_data = data[1].read(columns=['TARGETID', 'DESI_TARGET', 'MWS_TARGET', 'BGS_TARGET', 'SUBPRIORITY', 'NUMOBS_INIT', 'PRIORITY_INIT', 'RA', 'DEC', 'HPXPIXEL', 'BRICKNAME'])
    full_mtl = desitarget.mtl.make_mtl(tmp_data, obscon)

    ii_mtl_dark = (full_mtl['OBSCONDITIONS'] & obsconditions.DARK)!=0
    ii_mtl_gray = (full_mtl['OBSCONDITIONS'] & obsconditions.GRAY)!=0
    ii_mtl_bright = (full_mtl['OBSCONDITIONS'] & obsconditions.BRIGHT)!=0
    ii_north = (full_mtl['RA']>85) & (full_mtl['RA']<300) & (full_mtl['DEC']>-15)
    
    if obscon == 'DARK|GRAY':
        str_obscon = 'dark_gray'

    
    mtl_file = "targets/"+dr+"_mtl_"+str_obscon+"_"+cap+".fits"

    initial_mtl_file = "targets/subset_"+dr+"_mtl_"+str_obscon+"_"+cap+".fits"
        
    initial_std_file = "targets/subset_"+dr+"_std.fits"
    
    
    if obscon == 'DARK|GRAY' and cap == 'NGC':
        sel_mtl = full_mtl[(ii_mtl_dark | ii_mtl_gray) & ii_north]
        print("Writing subset "+obscon+" in the northern cap")
    else:
        print('THIS SELECTION WAS NOT DEFINED YET!')    
    sel_mtl.write(outdir+mtl_file, overwrite=True)
    
    print('cutting to ra,dec selection',ra_min,ra_max,dec_min,dec_max)
    mtl_data = Table.read(outdir+mtl_file)
    subset_ii = ra_dec_subset(mtl_data,ra_min,ra_max,dec_min,dec_max)
    mtl_data[subset_ii].write(outdir+initial_mtl_file, overwrite=True)
    
    del mtl_data
    
    print('doing standards')
    mtl_data = Table.read(outdir+initial_mtl_file)
    std_mask = desi_mask.STD_FAINT | desi_mask.STD_WD | desi_mask.STD_BRIGHT
    print('STDMASK', std_mask)
    std_ii = (mtl_data['DESI_TARGET'] & std_mask)!=0
    print(len(std_ii), np.count_nonzero(std_ii))
    mtl_data[std_ii].write(outdir+initial_std_file, overwrite=True)
    
    del mtl_data
    
    return initial_mtl_file
    

def write_initial_std_file(initial_mtl_file, initial_std_file):
    mtl_data = Table.read(outdir+initial_mtl_file)
    std_mask = desi_mask.STD_FAINT | desi_mask.STD_WD | desi_mask.STD_BRIGHT
    print('STDMASK', std_mask)
    std_ii = (mtl_data['DESI_TARGET'] & std_mask)!=0
    print(len(std_ii), np.count_nonzero(std_ii))
    mtl_data[std_ii].write(outdir+initial_std_file, overwrite=True)

def write_initial_sky_file(initial_sky_file,sky_data_file='', ra_min=130, ra_max=180, dec_min=-10, dec_max=40,outdir=''):
    print('reading sky')
    #sky_data = fitsio.read(sky_data_file)
    #subset_ii = ra_dec_subset(sky_data,ra_min,ra_max,dec_min,dec_max)
    #sd = sky_data[subset_ii]
    #print('writing sky')
    #outfd = fitsio.FITS(outdir+initial_sky_file, "rw")
    #outfd.write(None, header=None, extname="PRIMARY")
    #outfd.write(sd, header=None, extname="TARGETS")
    #outfd.close()
    #print('done writing sky')
    fd = fitsio.FITS(sky_data_file)
    fdata = np.array(fd[1].read())#(columns=keep_columns))

    inside = np.where(
        np.logical_and(
            np.logical_and((fdata["RA"] > ra_min), (fdata["RA"] < ra_max)),
            np.logical_and((fdata["DEC"] > dec_min), (fdata["DEC"] < dec_max))
        )
    )[0]
    sky_mtl = fdata[inside]


    fd.close()
    del fd

    # Sanity check that these are all sky, supp_sky, or bad_sky

    print("{} input targets in sky file".format(len(sky_mtl)))

    sky_sky_rows = np.where(
        np.bitwise_and(sky_mtl["DESI_TARGET"], desi_mask["SKY"].mask)
    )[0]

    print("  {} SKY targets".format(len(sky_sky_rows)))

    sky_suppsky_rows = np.where(
        np.bitwise_and(sky_mtl["DESI_TARGET"], desi_mask["SUPP_SKY"].mask)
    )[0]

    print("  {} SUPP_SKY targets".format(len(sky_suppsky_rows)))

    sky_badsky_rows = np.where(
        np.bitwise_and(sky_mtl["DESI_TARGET"], desi_mask["BAD_SKY"].mask)
    )[0]

    print("  {} BAD_SKY targets".format(len(sky_badsky_rows)))

#     sky_mask = 0
#     sky_mask |= desi_mask["SKY"].mask
#     sky_mask |= desi_mask["SUPP_SKY"].mask
#     sky_mask |= desi_mask["BAD_SKY"].mask


    if os.path.isfile(initial_sky_file):
        os.remove(initial_sky_file)
    with fitsio.FITS(initial_sky_file, "rw") as fd:
        fd.write(sky_mtl)


def write_initial_truth_file(initial_truth_file,initial_mtl_file='',outdir='',pixweight_file=''):
    import desitarget.mock.mockmaker as mb
    from desitarget.targetmask import desi_mask, bgs_mask, mws_mask
    

    is_lya_qso = assign_lya_qso(outdir+initial_mtl_file, pixweight_file)
    
    targets = Table.read(outdir+initial_mtl_file)
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

    # Check that all targets have been assigned to a class
    iii = truth['MOCKID']==0
    assert np.count_nonzero(iii)==0
    
    print('writing truth')
    truth.write(outdir+initial_truth_file, overwrite=True)
    print('done truth')
    
def prepare_tiles():
    tiles = Table(desimodel.io.load_tiles())

    ii_tiles = tiles['PROGRAM'] != 'BRIGHT'
    ii_subset = ra_dec_subset(tiles) 

    
    tilefile = 'footprint/subset_tiles.fits'
    tiles[ii_tiles&ii_subset].write(outdir+tilefile, overwrite='True')
    tiles = Table.read(outdir+tilefile)

    ii_gray = tiles['PROGRAM']=='GRAY'
    ii_dark_0 = (tiles['PROGRAM']=='DARK') & (tiles['PASS']==0)
    ii_dark_1 = (tiles['PROGRAM']=='DARK') & (tiles['PASS']==1)
    ii_dark_2 = (tiles['PROGRAM']=='DARK') & (tiles['PASS']==2)
    ii_dark_3 = (tiles['PROGRAM']=='DARK') & (tiles['PASS']==3)

    footprint = dict()
    footprint['gray'] = tiles[ii_gray]
    footprint['dark0'] = tiles[ii_dark_0]
    footprint['dark1'] = tiles[ii_dark_1]
    footprint['dark2'] = tiles[ii_dark_2]
    footprint['dark3'] = tiles[ii_dark_3]

    footprint['gray'].write(outdir+'footprint/subset_gray.fits', overwrite=True)
    footprint['dark0'].write(outdir+'footprint/subset_dark0.fits', overwrite=True)
    footprint['dark1'].write(outdir+'footprint/subset_dark1.fits', overwrite=True)
    vstack([footprint['dark2'], footprint['dark3']]).write(outdir+'footprint/subset_dark2_dark3.fits', overwrite=True)
    vstack([footprint['dark1'], footprint['dark2'], footprint['dark3']]).write(outdir+'footprint/subset_dark1_dark2_dark3.fits', overwrite=True)
    vstack([footprint['dark0'], footprint['dark1'], footprint['dark2'], footprint['dark3']]).write(outdir+'footprint/subset_dark0_dark1_dark2_dark3.fits', overwrite=True)
    vstack([footprint['gray'], footprint['dark0'], footprint['dark1'], footprint['dark2'], footprint['dark3']]).write(outdir+'footprint/subset_gray_dark0_dark1_dark2_dark3.fits', overwrite=True)

def create_multi_footprint(sim_path, footprint_path, cadence=28,outdir='',ramin=0,ramax=10,decmin=0,decmax=10):
    
    # load exposures and tiles
    #exposures = Table.read(os.path.join(sim_path,'exposures.fits'), hdu=1)
    exposures = fitsio.read(os.path.join(sim_path,'exposures.fits'), ext=1)
    tiles = desimodel.io.load_tiles() #seems slow
    #tfn  = os.getenv('DESIMODEL')+'/data/footprint/desi-tiles.fits'
    #tiles = fitsio.read(tfn)
    
    # select tiles to be dark+gray in a special region of the sky
    ii_subset = ra_dec_subset(tiles,ramin,ramax,decmin,decmax) 
    tiles = tiles[ii_subset]
    not_bright = tiles['PROGRAM']!='BRIGHT'
    dark_gray_tiles = tiles[not_bright]
    
    # create a "month" of "cadence" days.
    #exposures['MONTH'] = np.int_((exposures['MJD']-exposures['MJD'].min())/cadence)
    emonth = np.int_((exposures['MJD']-exposures['MJD'].min())/cadence)
   
    all_tiles_in_month = {}
    #month_id = list(set(exposures['MONTH']))
    month_id = list(set(emonth))
    month_id.sort()
    print(month_id)
    subsetnames = []
    for month in month_id:
        # gather all tiles available in a month from the surveysim file
        #all_tiles_in_month[month] = list(set(exposures['TILEID'][exposures['MONTH']==month]))
        all_tiles_in_month[month] = list(set(exposures['TILEID'][emonth==month]))
        # check that the available tiles are in the subset of tiles we are interested in
        ii = np.in1d(dark_gray_tiles['TILEID'], all_tiles_in_month[month])
        n_tiles = np.count_nonzero(ii)
        print(month, len(month_id), n_tiles)
        # write those tiles in a single gile
        if n_tiles > 0:
            table_tiles = Table(dark_gray_tiles[ii])
            subsetname = '{:02d}'.format(month)
            tilefile = os.path.join(footprint_path, 'subset_{}.fits'.format(subsetname))
            subsetnames.append(subsetname)
            print('writing to '+outdir+tilefile)
            table_tiles.write(outdir+tilefile, overwrite=True)
            print('wrote')
    return subsetnames
    
def consolidate_favail(fba_files):
    # getting all the targetids of the assigned fibers
    print('reading individual fiberassign files')
    favail = list()
    for i_tile, tile_file in enumerate(fba_files):
        if i_tile%50 ==0:
            print(i_tile)
        id_favail, header = fits.getdata(tile_file, 'FAVAIL', header=True)
        favail.extend(id_favail['TARGETID'])
    return list(set(favail))
    
def run_strategy(footprint_names, pass_names, obsconditions, strategy, initial_mtl_file, initial_sky_file, initial_std_file, 
                 fiberassign_script='fiberassign_legacy', legacy=None,outdir='',initial_truth_file=''):
    for i_pass in range(len(footprint_names)-1):
    
        footprint_name = footprint_names[i_pass]
        old_pass_name = pass_names[i_pass-1]
        pass_name = pass_names[i_pass]
        new_pass_name = pass_names[i_pass+1]
    
        os.makedirs(outdir+'{}/fiberassign_{}'.format(strategy, pass_name), exist_ok=True)
        os.makedirs(outdir+'{}/targets'.format(strategy), exist_ok=True)
        os.makedirs(outdir+'{}/zcat'.format(strategy), exist_ok=True)

    
        assign_footprint_filename = outdir+'footprint/subset_{}.fits'.format(footprint_name)
        zcat_footprint_filename = outdir+'footprint/subset_{}.fits'.format(pass_name)
        fiberassign_dir = outdir+'{}/fiberassign_{}/'.format(strategy, pass_name)
        mtl_filename = outdir+'{}/targets/{}_subset_dr8_mtl_dark_gray_NGC.fits'.format(strategy, pass_name)
        new_mtl_filename = outdir+'{}/targets/{}_subset_dr8_mtl_dark_gray_NGC.fits'.format(strategy, new_pass_name)
        old_zcat_filename = outdir+'{}/zcat/{}_zcat.fits'.format(strategy, old_pass_name)
        zcat_filename = outdir+'{}/zcat/{}_zcat.fits'.format(strategy, pass_name)
        sky_filename = outdir+'targets/subset_dr8_sky.fits'
    
        if i_pass == 0:
            shutil.copyfile(outdir+initial_mtl_file, mtl_filename)
        
    
        # Run fiberassign
 
        if legacy==True:
            cmd = '{} --mtl {} --sky {} --std {}'.format(fiberassign_script, mtl_filename, initial_sky_file, initial_std_file)
            cmd += ' --footprint {} --outdir {} --overwrite '.format(assign_footprint_filename, fiberassign_dir)
            cmd += ' --fibstatusfile fiberstatus.ecsv --starmask 60129542144'
        if legacy==False:
            cmd = 'fiberassign --mtl {} --sky {} '.format(mtl_filename,sky_filename)
            cmd +=' --footprint {} --outdir {} --overwrite'.format(assign_footprint_filename, fiberassign_dir)
            
        print(cmd)
        os.system(cmd)
    
        # Gather fiberassign files
        fba_files = np.sort(glob.glob(os.path.join(fiberassign_dir,"fiberassign*.fits")))

        # remove tilefiles that are not in the list of tiles to build zcat
        footprint = Table.read(zcat_footprint_filename)
        to_keep = []
        for i_file, fba_file in enumerate(fba_files):
            fibassign, header = fits.getdata(fba_file, header=True)
            tileid = int(header['TILEID'])
            if tileid in footprint['TILEID']:
                print(tileid, 'in list', outdir+zcat_footprint_filename)
                print('keeping {}'.format(fba_file))
                to_keep.append(i_file)
            else:
                print('renaming {}'.format(fba_file))
                fiberassign_file = fba_file.replace('fiberassign-', 'fba_')
                renamed_file = fiberassign_file.replace('.fits', '_unused.fits')
                print('renaming', fba_file, renamed_file)
                os.rename(fba_file, renamed_file)
            
        fba_files = fba_files[to_keep]
        print('Files to keep', len(fba_files))
    
        # Read targets and truth
        targets = Table.read(mtl_filename)
        truth = Table.read(outdir+initial_truth_file)
    
        # Compute zcat
        if i_pass==0:
            zcat = desisim.quickcat.quickcat(fba_files, targets, truth, perfect=True)
        else:
            old_zcat = Table.read(old_zcat_filename)
            zcat = desisim.quickcat.quickcat(fba_files, targets, truth, zcat=old_zcat, perfect=True)        
    
        zcat.write(zcat_filename, overwrite=True)
        mtl = desitarget.mtl.make_mtl(targets, obsconditions[i_pass], zcat=zcat)
        mtl.write(new_mtl_filename, overwrite=True)

 
