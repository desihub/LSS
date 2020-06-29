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

# (Change this whenever target samples update)



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


def mkmtl_dt(obscon="DARK|GRAY",target_ra_min=0,target_ra_max=360,target_dec_min=-90,target_dec_max=90,outf='/global/cscratch1/sd/ajross/fiberassigntest/fiducialtargets/temp/mtl.fits',infl='/project/projectdirs/desi/users/ajross/dr8tar/target_science_sample.fits'):

    fd = fitsio.FITS(infl, "r")
    fdata = fd[1].read()
    inside = np.where(
        np.logical_and(
            np.logical_and((fdata["RA"] > target_ra_min), (fdata["RA"] < target_ra_max)),
            np.logical_and((fdata["DEC"] > target_dec_min), (fdata["DEC"] < target_dec_max))
        )
    )[0]
    target_data = fdata[inside]

    mtl = make_mtl(target_data,obscon)

    fd = fitsio.FITS(outf, "rw")
    fd.write(None, header=None, extname="PRIMARY")
    fd.write(mtl, header=None, extname="TARGETS")
    fd.close()


def mkmtl_tk():
    '''
    initially copied from https://github.com/desihub/tutorials/blob/master/FiberAssignAlgorithms_Part2.ipynb
    should probably use desitarget function instead?
    '''
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
    targets_raw = fd[1].read(columns=keep_columns)

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