'''
Some functions to make it easy to generate mtl files, etc., for use in fiberassign tests

'''

import os
import sys
from collections import OrderedDict
import shutil
from datetime import datetime
import glob
import re

import numpy as np
from numpy.lib.recfunctions import append_fields

import matplotlib.pyplot as plt


from scipy.spatial import KDTree

import fitsio

from desimodel.io import findfile as dm_findfile
from desimodel.io import load_tiles as dm_load_tiles

from desitarget.targetmask import desi_mask, obsconditions

from desitarget.mtl import make_mtl

from fiberassign.targets import (
    Targets,
    TargetTree,
    TargetsAvailable,
    LocationsAvailable,
    load_target_table,
    default_target_masks,
    TARGET_TYPE_SCIENCE, 
    TARGET_TYPE_SKY,
    TARGET_TYPE_SUPPSKY,
    TARGET_TYPE_STANDARD
)


def mktarfile(target_ra_min=0,target_ra_max=360,target_dec_min=-90,target_dec_max=90,dr ='dr8',tarver = '0.39.0',outdir='/global/cscratch1/sd/ajross/fiberassigntest/fiducialtargets/temp/',prog='dark'):

    # First select the science targets

    input_dir = "/global/cfs/projectdirs/desi/target/catalogs/"+dr+"/"+tarver+"/targets/main/resolve/"+prog

    print('working with target data in '+input_dir+' IS THAT CORRECT???')

    sky_dir = "/global/cfs/projectdirs/desi/target/catalogs/"+dr+"/"+tarver+"/skies"

    print('working with skies data in '+sky_dir+' IS THAT CORRECT???')


    input_files = glob.glob(os.path.join(input_dir, "*.fits"))

    target_data = []

    for file in input_files:
        print("Working on {}".format(os.path.basename(file)), flush=True)
        fd = fitsio.FITS(file, "r")
        fdata = fd[1].read()
        inside = np.where(
            np.logical_and(
                np.logical_and((fdata["RA"] > target_ra_min), (fdata["RA"] < target_ra_max)),
                np.logical_and((fdata["DEC"] > target_dec_min), (fdata["DEC"] < target_dec_max))
            )
        )[0]
        target_data.append(fdata[inside])
        fd.close()

    target_data = np.concatenate(target_data)

    out_file = outdir+"target_science_sample.fits" 
    if os.path.isfile(out_file):
        os.remove(out_file)

    fd = fitsio.FITS(out_file, "rw")
    fd.write(None, header=None, extname="PRIMARY")
    fd.write(target_data, header=None, extname="TARGETS")
    fd.close()
    
    del target_data

    # Now select the sky targets

    print("Working on sky...", flush=True)

    input_files = glob.glob(os.path.join(sky_dir, "*.fits"))

    sky_data = []

    for file in input_files:
        print("Working on {}".format(os.path.basename(file)), flush=True)
        fd = fitsio.FITS(file, "r")
        fdata = fd[1].read()
        inside = np.where(
            np.logical_and(
                np.logical_and((fdata["RA"] > target_ra_min), (fdata["RA"] < target_ra_max)),
                np.logical_and((fdata["DEC"] > target_dec_min), (fdata["DEC"] < target_dec_max))
            )
        )[0]
        sky_data.append(fdata[inside])
        fd.close()

    sky_data = np.concatenate(sky_data)


    out_file = outdir+"target_sky_sample.fits" 
    if os.path.isfile(out_file):
        os.remove(out_file)


    outfd = fitsio.FITS(out_file, "rw")
    outfd.write(None, header=None, extname="PRIMARY")
    outfd.write(sky_data, header=None, extname="TARGETS")
    outfd.close()

    fd.close()



def mkmtl(obscon="DARK|GRAY",target_ra_min=0,target_ra_max=360,target_dec_min=-90,target_dec_max=90,outdir='/global/cscratch1/sd/ajross/fiberassigntest/fiducialtargets/temp/',target_sample='/project/projectdirs/desi/users/ajross/dr8tar/target_science_sample.fits'):
    '''
    initially copied from https://github.com/desihub/tutorials/blob/master/FiberAssignAlgorithms_Part2.ipynb
    
    '''
    
    science_file = outdir + 'mtl_science.fits'
    std_file = outdir + 'mtl_std.fits'

    
    # Load the raw science / standard target sample and prune columns

    keep_columns = [
        'TARGETID', 
        'RA', 
        'DEC',
        'RA_IVAR',
        'DEC_IVAR',
        'PMRA',
        'PMDEC',
        'PMRA_IVAR',
        'PMDEC_IVAR',
        'DESI_TARGET', 
        'BGS_TARGET', 
        'MWS_TARGET', 
        'SUBPRIORITY', 
        'BRICKNAME',
        'BRICKID',
        'BRICK_OBJID',
        'PRIORITY_INIT', 
        'NUMOBS_INIT'
    ]

    fd = fitsio.FITS(target_sample)
    fdata = fd[1].read(columns=keep_columns)

    inside = np.where(
        np.logical_and(
            np.logical_and((fdata["RA"] > target_ra_min), (fdata["RA"] < target_ra_max)),
            np.logical_and((fdata["DEC"] > target_dec_min), (fdata["DEC"] < target_dec_max))
        )
    )[0]
    targets_raw = fdata[inside]


    # Get the default target masks for this target file

    (filesurvey, 
     filecol, 
     def_sciencemask, 
     def_stdmask, 
     def_skymask, 
     def_suppskymask,
     def_safemask, 
     def_excludemask) = default_target_masks(targets_raw)

    print("Detected targets for survey '{}', using bitfield column '{}'".format(filesurvey, filecol))

    # Force our science and std masks to a more restrictive set.  Only keep ELG, LRG and QSO targets.
    # Cut any targets with multiple of those set.

    science_mask = 0
    science_mask |= desi_mask["LRG"].mask
    science_mask |= desi_mask["ELG"].mask
    science_mask |= desi_mask["QSO"].mask

    std_mask = 0
    std_mask |= desi_mask["STD_FAINT"].mask
    std_mask |= desi_mask["STD_WD"].mask
    std_mask |= desi_mask["STD_BRIGHT"].mask
    
    elg_rows = np.where(
        np.logical_and(
            np.logical_and(
                np.logical_and(
                    np.bitwise_and(targets_raw["DESI_TARGET"], desi_mask["ELG"].mask),
                    np.logical_not(
                        np.bitwise_and(targets_raw["DESI_TARGET"], std_mask)
                    )
                ),
                np.logical_not(
                    np.bitwise_and(targets_raw["DESI_TARGET"], desi_mask["QSO"].mask)
                )
            ),
            np.logical_not(
                np.bitwise_and(targets_raw["DESI_TARGET"], desi_mask["LRG"].mask)
            )
        )
    )[0]

    qso_rows = np.where(
        np.logical_and(
            np.logical_and(
                np.logical_and(
                    np.bitwise_and(targets_raw["DESI_TARGET"], desi_mask["QSO"].mask),
                    np.logical_not(
                        np.bitwise_and(targets_raw["DESI_TARGET"], std_mask)
                    )
                ),
                np.logical_not(
                    np.bitwise_and(targets_raw["DESI_TARGET"], desi_mask["ELG"].mask)
                )
            ),
            np.logical_not(
                np.bitwise_and(targets_raw["DESI_TARGET"], desi_mask["LRG"].mask)
            )
        )
    )[0]

    lrg_rows = np.where(
        np.logical_and(
            np.logical_and(
                np.logical_and(
                    np.bitwise_and(targets_raw["DESI_TARGET"], desi_mask["LRG"].mask),
                    np.logical_not(
                        np.bitwise_and(targets_raw["DESI_TARGET"], std_mask)
                    )
                ),
                np.logical_not(
                    np.bitwise_and(targets_raw["DESI_TARGET"], desi_mask["QSO"].mask)
                )
            ),
            np.logical_not(
                np.bitwise_and(targets_raw["DESI_TARGET"], desi_mask["ELG"].mask)
            )
        )
    )[0]
    
    n_elg = len(elg_rows)
    n_qso = len(qso_rows)
    n_lrg = len(lrg_rows)

    science_rows = np.concatenate([elg_rows, qso_rows, lrg_rows])

    std_rows = np.where(
        np.logical_and(
            np.bitwise_and(targets_raw["DESI_TARGET"], std_mask),
            np.logical_not(
                np.bitwise_and(targets_raw["DESI_TARGET"], science_mask)
            )
        )
    )[0]

    print(
        "Using {} science and {} standards from input catalog".format(
            len(science_rows),
            len(std_rows)
        )
    )

    # Split out the science and standard targets, although this is actually not necessary for passing
    # to fiberassign.

    science_targets = np.array(targets_raw[science_rows])

    std_targets = np.array(targets_raw[std_rows])

    # Close the input fits file so it doesn't take up extra memory
    del targets_raw
    fd.close()
    del fd

    # We have concatenated the 3 target types in the new table, so now the rows are
    # different:
    elg_rows = np.arange(n_elg, dtype=np.int64)
    qso_rows = np.arange(n_qso, dtype=np.int64) + n_elg
    lrg_rows = np.arange(n_lrg, dtype=np.int64) + n_elg + n_qso

    # Make the MTLs

    science_mtl = make_mtl(science_targets, "DARK|GRAY").as_array()
    if len(science_mtl) != len(science_targets):
        print("WARNING:  science MTL has {} rows, input has {}".format(len(science_mtl), len(science_targets)))

    std_mtl = make_mtl(std_targets, "DARK|GRAY").as_array()
    if len(std_mtl) != len(std_targets):
        print("WARNING:  standards MTL has {} rows, input has {}".format(len(std_mtl), len(std_targets)))

    # Delete the large intermediate arrays
    
    del science_targets
    del std_targets
    
    # Write MTLs

    if os.path.isfile(science_file):
        os.remove(science_file)
    with fitsio.FITS(science_file, "rw") as fd:
        fd.write(science_mtl)

    if os.path.isfile(std_file):
        os.remove(std_file)
    with fitsio.FITS(std_file, "rw") as fd:
        fd.write(std_mtl)    

    print("{} science targets".format(len(science_mtl)))
    print("    {} ELG targets".format(len(elg_rows)))
    print("    {} QSO targets".format(len(qso_rows)))
    print("    {} LRG targets".format(len(lrg_rows)))
    print("{} std targets".format(len(std_mtl)))

    # We'll be loading later science MTLs as we go through the survey, so delete that now.
    # the standards are constant so we'll keep those in memory.

    del science_mtl
    
def mkmtl_sky(target_ra_min=0,target_ra_max=360,target_dec_min=-90,target_dec_max=90,outdir='/global/cscratch1/sd/ajross/fiberassigntest/fiducialtargets/temp/',target_sample='/project/projectdirs/desi/users/ajross/dr8tar/target_sky_sample.fits'):
	# Now create the skies file.

	std_file = outdir + 'mtl_sky.fits'
	
	keep_columns = [
		'TARGETID', 
		'RA', 
		'DEC', 
		'DESI_TARGET', 
		'BGS_TARGET', 
		'MWS_TARGET', 
		'SUBPRIORITY', 
		'BRICKNAME',
		'BRICKID',
		'BRICK_OBJID',
		'APFLUX_G',
		'APFLUX_R',
		'APFLUX_Z',
		'APFLUX_IVAR_G',
		'APFLUX_IVAR_R',
		'APFLUX_IVAR_Z',
		'OBSCONDITIONS'
	]


    fdata = np.array(fd[1].read(columns=keep_columns))

    inside = np.where(
        np.logical_and(
            np.logical_and((fdata["RA"] > target_ra_min), (fdata["RA"] < target_ra_max)),
            np.logical_and((fdata["DEC"] > target_dec_min), (fdata["DEC"] < target_dec_max))
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

	sky_mask = 0
	sky_mask |= desi_mask["SKY"].mask
	sky_mask |= desi_mask["SUPP_SKY"].mask
	sky_mask |= desi_mask["BAD_SKY"].mask

	sky_unknown_rows = np.where(
		np.logical_not(
			np.bitwise_and(sky_mtl["DESI_TAR 
		)
    )[0]

	print("  {} targets are not one of the 3 recognized types".format(len(sky_unknown_rows)))

	if os.path.isfile(sky_file):
		os.remove(sky_file)
	with fitsio.FITS(sky_file, "rw") as fd:
		fd.write(sky_mtl)   