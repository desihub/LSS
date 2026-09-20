# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
LSS.sky_maps
============

Routines for building weight maps from randoms, etc., for systematics
"""
import os
import re
import fitsio
import numpy as np
from time import time
import healpy as hp
from glob import glob

from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.wcs import WCS

from desitarget.io import read_targets_header, find_target_files, read_targets_in_box
from desitarget.geomask import match, nside2nside
from desitarget.gaiamatch import get_gaia_dir, gaia_psflike
from desitarget.randoms import dr8_quantities_at_positions_in_a_brick, get_dust
from desitarget.targets import resolve
from desitarget import __version__ as desitarget_version
from desitarget.internal import sharedmem

# ADM the DESI default logger.
from desiutil.log import get_logger
from desiutil import brick

# ADM general bitmasks for LSS
from LSS.masks import lrg_mask, elg_mask, skymap_mask

# ADM initialize the DESI default logger.
log = get_logger()

# ADM start the clock.
start = time()

mapdt = [
    ('MAPNAME', 'O'), ('SUBDIR', 'O'), ('FILENAME', 'O'), ('NSIDE', 'i'),
    ('MAPTYPE', 'O'), ('COLNAME', 'O'), ('MASKCHECK', 'U3'), ('NESTED', '?'),
    ('GALACTIC', '?')
]
# ADM update with new maps in this global array on a new line.
# ADM 'NESTED' is True for nested HEALPix and False for the ring scheme.
# ADM `GALACTIC' is True for a map in Galactic coords, False for RA/Dec.
# -------
# ADM 'MASKCHECK' is a 3-character string containing the logical comparison used to define a value that IS masked;
# ADM for example if MASKCHECK is "> 1" then all values in the file that are greater than 1 ARE masked.
# MMM if Maskvalue is anything other than '' degraded pixel values may have contributions from masked areas.
# MMM if Maskvalue does not exist (no formally masked pixels), it defaults to '' for PIXMAP
# MMM WARNING: not all columns have the same dtype and range
# MMM WARNING: Need to decide on column with MaskValue
# MMM Provisional notes:
#     Halpha -   Temperature and error: any value f4, mask: u1 1 (ok) or above (masked)
#     Calibration - no column names as it is not a table but an image fits file; mask data == 0 same file.
#     EBVdustGaia - column names Recon_Mean, Recon_Variance, Recon_VarianceCorr f8 ; no masks
#     EBV_SGF     - column names: ebv f4, mask: status i4 (masked area if status > 0)
#     kappa       - values: index, real, imag;  mask :
maparray = np.array([
    ('HALPHA',                  'Halpha',                    'Halpha_fwhm06_0512.fits',  512,  'PIXMAP', 'TEMPERATURE',        '',  True, True),
    ('HALPHA_ERROR',            'Halpha',              'Halpha_error_fwhm06_0512.fits',  512,  'PIXMAP', 'ERROR',              '',  True, True),
    ('HALPHA_MASK',             'Halpha',               'Halpha_mask_fwhm06_0512.fits',  512, 'PIXMASK', 'MASK',            '> 1',  True, True),
    ('CALIB_G',            'calibration',                      'decam-ps1-0128-g.fits',  128,  'PIXMAP', 'NONE-IMAGE',         '', False, False),
    ('CALIB_R',            'calibration',                      'decam-ps1-0128-r.fits',  128,  'PIXMAP', 'NONE-IMAGE',         '', False, False),
    ('CALIB_Z',            'calibration',                      'decam-ps1-0128-z.fits',  128,  'PIXMAP', 'NONE-IMAGE',         '', False, False),
    ('CALIB_G_MASK',       'calibration',                      'decam-ps1-0128-g.fits',  128, 'PIXMASK', 'NONE-IMAGE',      '==0', False, False),
    ('CALIB_R_MASK',       'calibration',                      'decam-ps1-0128-r.fits',  128, 'PIXMASK', 'NONE-IMAGE',      '==0', False, False),
    ('CALIB_Z_MASK',       'calibration',                      'decam-ps1-0128-z.fits',  128, 'PIXMASK', 'NONE-IMAGE',      '==0', False, False),
    ('EBV_CHIANG_SFDCORR',         'EBV',    'Chiang23_SFD_corrected_hp2048_nest.fits', 2048,  'PIXMAP', 'T',                  '',  True, True),
    ('EBV_CHIANG_LSS_MASK',        'EBV',             'Chiang23_mask_hp2048_nest.fits', 2048, 'PIXMASK', 'T',               '==0',  True, True),
    ('EBV_CHIANG_LSS_ERROR',       'EBV',                    'Chiang23_lss_error.fits', 2048,  'PIXMAP', 'T',                  '',  True, True),
    ('EBV_CHIANG_STD',             'EBV',                      'Chiang23_std_ebv.fits', 2048,  'PIXMAP', 'T',                  '',  True, True),
    ('EBV_CHIANG_LSS_INTENSITY',   'EBV',                'Chiang23_lss_intensity.fits', 2048,  'PIXMAP', 'T',                  '',  True, True),
    ('EBV_MPF_MEAN_FW15',          'EBV',                 'recon_fw15_final_mult.fits', 2048,  'PIXMAP', 'Recon_Mean',         '', False, True),
    ('EBV_MPF_MEAN_ZPTCORR_FW15',  'EBV',                 'recon_fw15_final_mult.fits', 2048,  'PIXMAP', 'Recon_Mean_ZptCorr', '', False, True),
    ('EBV_MPF_VAR_FW15',           'EBV',                 'recon_fw15_final_mult.fits', 2048,  'PIXMAP', 'Recon_Variance',     '', False, True),
    ('EBV_MPF_VARCORR_FW15',       'EBV',                 'recon_fw15_final_mult.fits', 2048,  'PIXMAP', 'Recon_VarianceCorr', '', False, True),
    ('EBV_MPF_MEAN_FW6P1',         'EBV',                'recon_fw6-1_final_mult.fits', 2048,  'PIXMAP', 'Recon_Mean',         '', False, True),
    ('EBV_MPF_MEAN_ZPTCORR_FW6P1', 'EBV',                'recon_fw6-1_final_mult.fits', 2048,  'PIXMAP', 'Recon_Mean_ZptCorr', '', False, True),
    ('EBV_MPF_VAR_FW6P1',          'EBV',                'recon_fw6-1_final_mult.fits', 2048,  'PIXMAP', 'Recon_Variance',     '', False, True),
    ('EBV_MPF_VARCORR_FW6P1',      'EBV',                'recon_fw6-1_final_mult.fits', 2048,  'PIXMAP', 'Recon_VarianceCorr', '', False, True),
    ('EBV_SGF14',                  'EBV',                        'ps1-ebv-4.5kpc.fits',  512,  'PIXMAP', 'ebv',                '', False, True),
    ('EBV_SGF14_MASK',             'EBV',                        'ps1-ebv-4.5kpc.fits',  512, 'PIXMASK', 'status',          '< 0', False, True),
    ('EBV_ZHOU_DESI_GR',           'EBV',                'Zhou23_v1_desi_ebv_256.fits',  256,  'PIXMAP', 'EBV_DESI_GR',        '', False, False),
    ('EBV_ZHOU_DESI_RZ',           'EBV',                'Zhou23_v1_desi_ebv_256.fits',  256,  'PIXMAP', 'EBV_DESI_RZ',        '', False, False),
    ('BETA_ML',                    'EBV', 'COM_CompMap_dust-commander_0256_R2.00.fits',  256,  'PIXMAP', 'BETA_ML',            '',  True, True),
    ('BETA_MEAN',                  'EBV', 'COM_CompMap_dust-commander_0256_R2.00.fits',  256,  'PIXMAP', 'BETA_MEAN',          '',  True, True),
    ('BETA_RMS',                   'EBV', 'COM_CompMap_dust-commander_0256_R2.00.fits',  256,  'PIXMAP', 'BETA_RMS',           '',  True, True),
    ('HI',                         'NHI',                            'NHI_HPX.fits.gz', 1024,  'PIXMAP', 'NHI',                '', False, True),
    ('KAPPA_PLANCK',             'kappa',                               'dat_klm.fits', 2048,  'ALMMAP', 'NONE-3col',          '', False, True),
    ('KAPPA_PLANCK_MASK',        'kappa',                               'mask.fits.gz', 2048, 'PIXMASK', 'I',               '==0', False, True),
    ('FRACAREA',        'pixweight-dark',                      'pixweight-1-dark.fits',  256,  'PIXMAP', 'FRACAREA',           '',  True, False),
    ('STARDENS',              'stardens',                              'stardens.fits',  512,  'PIXMAP', 'STARDENS',           '',  True, False),
    ('ELG',             'pixweight-dark',                      'pixweight-1-dark.fits',  256,  'PIXMAP', 'ELG',                '',  True, False),
    ('LRG',             'pixweight-dark',                      'pixweight-1-dark.fits',  256,  'PIXMAP', 'LRG',                '',  True, False),
    ('QSO',             'pixweight-dark',                      'pixweight-1-dark.fits',  256,  'PIXMAP', 'QSO',                '',  True, False),
    ('BGS_ANY',       'pixweight-bright',                    'pixweight-1-bright.fits',  256,  'PIXMAP', 'BGS_ANY',            '',  True, False)
    ], dtype=mapdt)


def sanity_check_map_array():
    """Convenience function to check the format of the map_array global.
    """
    log.info("Running sanity checks on maparray...")

    for skymap in maparray:

        mapname = skymap['MAPNAME']

        # ADM check named files/directories exist (at least at NERSC).
        lssmapdir = os.getenv("LSS_MAP_DIR")
        if lssmapdir is not None:
            fn = os.path.join(lssmapdir, skymap["SUBDIR"], skymap["FILENAME"])
            if not os.path.exists(fn):
                msg = "{} does not exist".format(fn)
                raise_myerror(msg)

        # MMM check nside is an integer.
        if not isinstance(skymap["NSIDE"].tolist(), int):
            msg = "NSIDE is not an integer in {}"
            log.critical(msg.format(mapname))
            raise ValueError(msg.format(mapname))

        # MMM perform a sanity check on options or maptype.
        if skymap['MAPTYPE'] not in ['PIXMAP', 'PIXMASK', 'ALMMAP']:
            msg = "There is NO acceptable value for MAPTYPE"
            log.critical(msg.format(mapname))
            raise ValueError(msg.format(mapname))

        # ADM check the conditionals in the MASKCHECK column.
        if skymap["MAPTYPE"] == "PIXMASK":
            parse_mask_check(np.empty(2), skymap["MASKCHECK"], check=True)
            # ADM check the name for the skymap mask makes sense.
            if not skymap["MAPNAME"][-5:] == "_MASK":
                msg = "Mask-maps need MAPNAMEs ending in _MASK; {} does not!"
                log.critical(msg.format(skymap["MAPNAME"]))
                raise ValueError(msg.format(skymap["MAPNAME"]))

        # ADM check somebody didn't include two maps with the same name.
        pixmap = maparray[maparray["MAPNAME"] == mapname]
        if len(pixmap) != 1:
            if len(pixmap) > 1:
                msg = "There are TWO maps in maparray that have MAPNAME={}!"
            # ADM check there's an entry in maparray for the passed map name.
            elif len(pixmap) < 1:
                msg = "There are NO maps in maparray that have MAPNAME={}!"
            log.critical(msg.format(mapname))
            raise ValueError(msg.format(mapname))

        # ADM an ALMMAP should never have COLNAME NONE_IMAGE, or we won't
        # ADM know how to read the map.
        if (skymap["MAPTYPE"] == "ALMMAP") & (skymap["COLNAME"] == "NONE-IMAGE"):
            msg = "COLNAME can't be 'NONE-IMAGE' if MAPTYPE is 'ALMMAP' (see {})"
            log.critical(msg.format(skymap["MAPNAME"]))
            raise ValueError(msg.format(skymap["MAPNAME"]))

    log.info("...maparray seems to be correctly formatted")

    return


def read_main_survey_targets(obscon):
    """Read in DESI Main Survey targets.

    Parameters
    ----------
    obscon : :class:`list`
        Pass "dark" to read in dark-time targets (e.g. ELGs, LRGs, QSOs)
        or "bright" to read in bright-time targets (BGS).

    Returns
    -------
    :class:`~numpy.ndarray`
        A numpy structured array of all DESI Main Survey targets for the
        passed observing condition (`obscon`) that contains columns "RA",
       "DEC", "DESI_TARGET" and "BGS_TARGET". Useful for limiting targets
       to specific target classes.

    Notes
    -----
    - If the environment variable $TARG_DIR is set, then that is used
      as the root directory to find the target files. Otherwise, the
      filenames are hardcoded with their locations at NERSC.
    """
    # ADM hard-code the root target directory, unless $TARG_DIR is set.
    targdir = os.environ.get('TARG_DIR')
    if targdir is None:
        targdir = "/global/cfs/cdirs/desi/target/catalogs"

    # ADM check whether the passed observing condition was valid.
    if obscon not in ["bright", "dark"]:
        msg = 'Allowed obscons are "bright" and "dark". You passed "{}"!'
        raise_myerror(msg.format(obscon))

    # ADM use desitarget I/O code to find the appropriate filename for
    # ADM Main Survey targets and the given observing conditions.
    filename = find_target_files(targdir, dr="9", flavor="targets",
                                 survey="main", obscon=obscon)

    # ADM find_target_files() defaults to the currenct version of
    # ADM desitarget, so we need to replace this with the version
    # ADM used for the main survey.
    filename = filename.replace(desitarget_version, "1.1.1")

    # ADM use desitarget I/O code to read in targets.
    targets = read_targets_in_box(filename, quick=True,
                                  columns=["RA", "DEC", "DESI_TARGET", "BGS_TARGET"])

    log.info("Read {} {}-time targets in {:.1f}s".format(
        len(targets), obscon, time()-start))

    return targets


def get_lss_map_dir(lssmapdir=None):
    """Function to grab the $LSS_MAP_DIR environment variable.

    Parameters
    ----------
    lssmapdir : :class:`str`, optional, defaults to $LSS_MAP_DIR
        If `lssmapdir` is passed, it's returned from this function. If it
        is not passed, the $LSS_MAP_DIR environment variable is returned.

    Returns
    -------
    :class:`str`
        If `lssmapdir` is passed, it is returned from this function. If
        it is not passed, the directory stored in the $LSS_MAP_DIR
        environment variable is returned.

    Notes
    -----
    - At NERSC, $LSS_MAP_DIR is typically:
      /global/cfs/cdirs/desi/survey/catalogs/external_input_maps
    """
    if lssmapdir is None:
        lssmapdir = os.environ.get('LSS_MAP_DIR')
        # ADM check that the $LSS_MAP_DIR environment variable is set.
        if lssmapdir is None:
            msg = "Pass lssmapdir or set the $LSS_MAP_DIR environment variable!"
            log.critical(msg)
            raise ValueError(msg)

    return lssmapdir


def get_lss_dir(lssmapdir=None, survey="main"):
    """Derive the LSS directory from the $LSS_MAP_DIR environment variable.

    Parameters
    ----------
    lssmapdir : :class:`str`, optional, defaults to $LSS_MAP_DIR
        The location of the directory that host the LSS maps. If `None`
        is passed then the $LSS_MAP_DIR environment variable is used.
    survey : :class:`str`, optional, defaults to "main"
        A survey phase. Either upper-case or lower-case can be passed.

    Returns
    -------
    :class:`str`
        The base directory for LSS products.

    Notes
    -----
    - At NERSC, the LSS directory is typically:
      /global/cfs/cdirs/desi/survey/catalogs/$survey/LSS
    """
    lssmapdir = get_lss_map_dir(lssmapdir=lssmapdir)

    # ADM convert between upper/lower case as needed.
    if "s" in survey:
        surv = survey.upper()
    else:
        surv = survey.lower()

    if surv not in ["SV1", "SV2", "SV3", "main"]:
        msg = "Survey: {} not recognized".format(survey)
        log.critical(msg)
        raise ValueError(msg)

    return os.path.join(os.path.dirname(lssmapdir), surv, "LSS")


def get_ls_brickmask_dir(tracer, lssmapdir=None):
    """Derive Legacy Surveys brickmasks directory from $LSS_MAP_DIR.

    Parameters
    ----------
    tracer : :class:`str`
        The name of a tracer that has masks in the brickmasks directory.
        Options are "ELG" and "LRG".
    lssmapdir : :class:`str`, optional, defaults to $LSS_MAP_DIR
        The location of the directory that host the LSS maps. If `None`
        is passed then the $LSS_MAP_DIR environment variable is used.

    Returns
    -------
    :class:`str`
        The base directory for the brickmask products for the `tracer`.

    Notes
    -----
    - At NERSC, the brickmask directory is typically:
      /global/cfs/cdirs/desi/survey/catalogs/brickmasks
    """
    lssmapdir = get_lss_map_dir(lssmapdir=lssmapdir)

    # ADM make sure the tracer names are upper-case
    tracer = tracer.upper()

    if tracer not in ["LRG", "ELG"]:
        msg = "Tracer: {} not recognized".format(tracer)
        log.critical(msg)
        raise ValueError(msg)

    version = "v1"
    if tracer == "LRG":
        version = "v1.1"

    basedir = os.path.dirname(lssmapdir)

    return os.path.join(basedir, "brickmasks", tracer, version)


def write_atomically(filename, data, extname=None, header=None):
    """Write a FITS file in an atomic fashion.

    Parameters
    ----------
    filename : :class:`str`
        The output file.
    data : :class:`~numpy.ndarray`
        The numpy structured array of data to write.
    extname, header optional
        Passed through to `fitsio.write()`. `header` can be either
        a FITShdr object or a dictionary.

    Returns
    -------
    Nothing, but writes the `data` to the `filename`.

    Notes
    -----
    - Always OVERWRITES existing files!
    - Always makes the `filename` directory if it doesn't exist.
    - By "in an atomic fashion" it is meant that files that died
      mid-write will be appended by ".tmp".
    """
    # ADM make the necessary directory if it doesn't exist.
    os.makedirs(os.path.dirname(filename), exist_ok=True)

    # ADM write the file atomically by making a .tmp file and moving it.
    fitsio.write(filename+'.tmp', data, extname=extname, header=header,
                 clobber=True)
    os.rename(filename+'.tmp', filename)

    # ADM note we successfully completed.
    log.info("Wrote {} objects to {}...t={:.1f}s".format(
        len(data), filename, time()-start))

    return


def ls_bitmask_one_brick(brickname, ra, dec, photsys, tracer, mxdir=None):
    """Slow look up of Legacy Surveys bitmask information for one brick.

    Parameters
    ----------
    brickname : :class:`str`
        The name of a Legacy Surveys brick, e.g. '0001m015'
    ra : :class:`float`
        The Right Ascension of a location (degrees).
    dec : :class:`float
        The Declination of a matching location (degrees).
    photsys : :class:`float`
        Pass "N" for locations in a northern Legacy Surveys brick, or
        "S" for locations in a southern brick.
    tracer : :class:`str`
        The name of a tracer that has masks in the brickmasks directory.
        Options are "ELG" and "LRG".
    mxdir : :class:`str`
        The full path to the brickmasks directory. Pass ``None`` to guess
        this from the $LSS_MAP_DIR environment variable.

    Returns
    -------
    :class:`~numpy.ndarray`
        The lrg_mask or elg_mask bitmasks at these locations for this
        tracer type. The meanings of the bitmasks can be loaded via
        from LSS.masks import lrg_mask, elg_mask

    Notes
    -----
    - Adapted from:
      https://github.com/desihub/LSS/blob/main/scripts/readwrite_pixel_bitmask.py
    - I (ADM) ended up not taking this approach, as it's so slow. But, I kept
      the code in case it's useful for later development.
    """
    # ADM assume a default location if mxdir was not passed.
    if mxdir is None:
        mxdir = get_brickmask_dir(tracer)

    # ADM convert "N" to "north" and "S" to "south".
    if photsys == "N":
        field = "north"
    elif photsys == "S":
        field = "south"
    else:
        msg = "photsys should be 'N' or 'S', not {}".format(photsys)
        log.critical(msg)
        raise ValueError(msg)

    fn = os.path.join(mxdir, '{}/coadd/{}/{}/{}-{}mask.fits.gz'.format(
        field, brickname[:3], brickname, brickname, tracer.lower()))

    img = fitsio.read(fn)

    img, hdr = fitsio.read(fn, header=True)
    w = WCS(hdr)

    coadd_x, coadd_y = w.wcs_world2pix(ra, dec, 0)
    coadd_x, coadd_y = np.round(coadd_x).astype(int), np.round(coadd_y).astype(int)

    mx = img[coadd_y, coadd_x]

    return mx


def ls_bitmask_for_randoms(randoms, ident, lssmapdir=None, survey="main"):
    """Assign Legacy-Surveys-based tracer bitmask information to randoms.

    Parameters
    ----------
    randoms : :class:`~numpy.ndarray`
        Random catalog, as made by, e.g. :func:`read_randoms()`.
    ident : :class:`~numpy.ndarray`
        Structured array with one column "IDENT" that has the same number
        of rows as the input `randoms`. Can also be generated by
        :func:`read_randoms()`. Used to associate `randoms` to the
        corresponding random catalogs on disk.
    lssmapdir : :class:`str`, optional, defaults to $LSS_MAP_DIR
        The location of the directory that host the LSS maps. If `None`
        is passed then the $LSS_MAP_DIR environment variable is used.
        Used to look up the LSS directory.
    survey : :class:`str`, optional, defaults to "main"
        A survey phase. Either upper-case or lower-case can be passed.

    Returns
    -------
    :class:`~numpy.ndarray`
        An array of integers populated for all Legacy Surveys tracers
        that have bitmasks (typically ELGs and LRGs). The populated
        values can be loaded via: from LSS.masks import skymap_mask

    Notes
    -----
    - Performs a cross-match on TARGETID to existing, pre-constructed
      per-random catalogs of bitmasks in the LSS directory. If these
      catalogs don't exist they'll need made (triggering an OSerror).
    """
    # ADM find the location of the pre-constructed bitmask catalogs.
    lss_dir = get_lss_dir(lssmapdir=lssmapdir, survey=survey)

    # ADM set up an output array. I use "i8" here because (as of the time
    # ADM of writing) fitsio does not support I/O for "u8" (uint64).
    done = np.zeros(len(randoms), dtype="i8")

    # ADM need to process once per different ident/random catalog.
    sidents = sorted(list(set(ident["IDENT"])))
    for sident in sidents:
        log.info("Working on random catalog {}...t={:.1f}s".format(
            sident, time()-start))
        # ADM loop through each tracer and populate a bitmask transposed
        # ADM from the LS representation to the skymap representation.
        for i, tracer in enumerate(["elg", "lrg"]):
            # ADM the file that contains the Legacy Surveys info.
            fn = "randoms-{}{}imask.fits".format(sident, tracer)
            fn = os.path.join(lss_dir, fn)
            lsmx = fitsio.read(fn)
            # ADM first-time-through, set up an array to hold the bitmask
            # ADM transposed from the LS to skymap representation.
            if i == 0:
                trans = np.zeros(len(lsmx), dtype="i8")
                targetids = lsmx["TARGETID"]
            else:
                # ADM a check that all of the files of tracers are
                # ADM ordered consistently on TARGETID.
                if np.any(lsmx["TARGETID"] != targetids):
                    msg = ("Files of tracers (e.g. {}) not ordered"
                           " consistently on TARGETID".format(fn))
                    log.critical(msg)
                    raise ValueError(msg)
            # ADM the name of the relevant column/bitmask.
            mxname = "{}_mask".format(tracer)
            log.info("Transposing bits from {} to skymap_mask..t={:.1f}s".format(
                mxname, time()-start))
            # ADM the bitmask itself.
            mx = globals()[mxname]
            # ADM loop through the relevant bitmask names and transpose
            # ADM from the LS bitmask to skymap bitmask representation.
            for nom in mx.names():
                altnom = "{}_{}".format(tracer.upper(), nom)
                # ADM the mask name will be the same in LS and skymap...
                if nom in skymap_mask.names():
                    snom = nom
                # ADM ...or it will be preceeded by the tracer name...
                elif altnom in skymap_mask.names():
                    snom = altnom
                # ADM ...or something went wrong when naming masks.
                else:
                    msg = "Can't transpose {} in {} to {} or {} in {}".format(
                        nom, mxname, nom, altnom, "skymap_mask")
                    log.critical(msg)
                    raise ValueError(msg)
                log.info("Transposing {} to {}".format(nom, snom))
                # ADM find instances of bit in the Legacy Surveys mask.
                is_bit_set = lsmx[mxname] & mx[nom] != 0
                # ADM set instances of this bit using the skymap mask.
                # ADM I cast as 'i8', here, because (as of the time of
                # ADM writing) fitsio does not support I/O for "u8".
                trans |= np.array(is_bit_set * skymap_mask[snom], dtype='i8')
        # ADM now match on TARGETID and assign mask bits.
        rii, tii = match(randoms["TARGETID"], targetids)
        done[rii] = trans[tii]

    return done


def wrap_pixmap(randoms, targets, nside=512, gaialoc=None):
    """HEALPix map from randoms (wrapper on desitarget.randoms.pixmap)

    Parameters
    ----------
    randoms : :class:`~numpy.ndarray` or `str`
        Filename of random catalog or catalog itself. Catalogs must have
        columns 'RA', 'DEC', 'EBV', 'PSFDEPTH_W1/W2/G/R/Z', 'NOBS_G/R/Z'
        'GALDEPTH_G/R/Z', 'PSFSIZE_G/R/Z', 'MASKBITS' and have been
        generated at the same density. If `randoms` is a list files will
        be concatenated in list-order.
    targets : :class:`~numpy.ndarray` or `str`
        Corresponding (same Data Release as `randoms`) file of targets,
        or name of a directory containing HEALPixel-split targets that
        can be read by :func:`desitarget.io.read_targets_in_box()`.
    nside : :class:`int`, optional, defaults to nside=512
        Resolution (HEALPix nside) at which to build the (NESTED) map.
    gaialoc : :class:`str`, optional, defaults to ``None``
        Name of a FITS file that already contains a column "STARDENS",
        which is simply read in. If ``None``, the stellar density is
        constructed from files in $GAIA_DIR.

    Returns
    -------
    :class:`~numpy.ndarray`
        An array of useful information that includes
            - HPXPIXEL: HEALPixel integers at the passed `nside`.
            - FRACAREA: Fraction of pixel with at least one observation
                        in any band. Made with :func:`pixweight()`.
            - STARDENS: The stellar density in a pixel from Gaia. Made
                        with :func:`stellar_density()`.
            - EBV: E(B-V) in pixel from the SFD dust map, from the
                   median of EBV values in the passed `randoms`.
            - PSFDEPTH_G, R, Z: PSF depth in the pixel, from the median
                                of PSFDEPTH values in `randoms`.
            - GALDEPTH_G, R, Z: Galaxy depth in the pixel, from the
                                median of GALDEPTH values in `randoms`.
            - PSFDEPTH_W1, W2: (AB PSF) depth in the pixel, from the
                               median of values in the passed `randoms`.
            - PSFSIZE_G, R, Z: Weighted average PSF FWHM, in arcsec, in
                               the pixel, from the median of PSFSIZE
                               values in the passed random catalog.
            - FRACAREA_X: Fraction of pixel with at least one observation
                          in any band with MASKBITS==X (bitwise OR, so,
                          e.g. if X=7 then fraction for 2^0 | 2^1 | 2^2).
            - One column for every bit that is returned by
              :func:`desitarget.QA._load_targdens()`. Each column
              contains the target density in the pixel.
    :class:`str`
        Survey to which `targets` corresponds, e.g., 'main', 'svX', etc.

    Notes
    -----
    - If `gaialoc` is ``None`` then $GAIA_DIR must be set.
    - Docstring mostly stolen from :func:`desitarget.randoms.pixmap()`.
    """
    # ADM desitarget function to wrap.
    from desitarget.randoms import pixmap

    return pixmap(randoms, targets, dens, nside=nside, gaialoc=gaialoc)


def write_pixmap(randoms, targets, hdr=None, nside=512, gaialoc=None,
                 outdir=None):
    """Write pixmap made by :func:`wrap_pixmap()`

    Parameters
    ----------
    randoms : :class:`~numpy.ndarray`
        Random catalog.
    targets : :class:`str`
        Corresponding (same Data Release as `randoms`) file of targets,
        or name of a directory containing HEALPixel-split targets that
        can be read by :func:`desitarget.io.read_targets_in_box()`. The
        file (or all files in the directory) must contain "OBSCON" in the
        header so the code can determine if we're working with dark-time
        or bright-time targets.
    hdr : :class:`dict` or `FITSHDR`
        Header to write to the pixweight file.
    nside : :class:`int`, optional, defaults to nside=512
        Resolution (HEALPix nside) at which to build the (NESTED) map.
    gaialoc : :class:`str`, optional, defaults to ``None``
        Name of a FITS file that already contains a column "STARDENS",
        which is simply read in. If ``None``, the stellar density is
        constructed from files in $GAIA_DIR.
    outdir : :class:`str`, optional, defaults to ``None``
        Name of output directory to which to write pixel map. If ``None``
        then default to the $LSS_MAP_DIR environment variable.

    Returns
    -------
    Nothing, but writes a pixel map to:
        `outdir`/pixweight_maps_all/pixweight-<obscon>.fits
    if the keyword SEED isn't provided in `hdr`, or, if it is:
        `outdir`/pixweight_maps_all/pixweight-<seed>-<obscon>.fits
    """
    # ADM read in the observing program from the target file header.
    obscon = read_targets_header(targets)["OBSCON"].lower()

    # ADM construct the output file.
    outfile = "pixweight-{}.fits".format(obscon)
    if hdr is not None:
        if "SEED" in hdr:
            outfile = outfile.replace("eight-", "eight-{}-".format(hdr["SEED"]))
    else:
        # ADM if no header was passed, we need to construct one.
        hdr = fitsio.FITSHDR()

    lss_map_dir = get_lss_map_dir(outdir)
    outfile = os.path.join(lss_map_dir, "pixweight_maps_all", outfile)

    # ADM augment the output header.
    hdr['GAIALOC'] = gaialoc
    hdr['HPXNSIDE'] = nside
    hdr['HPXNEST'] = True

    pixmap, survey = wrap_pixmap(randoms, targets, nside=nside, gaialoc=gaialoc)

    hdr["SURVEY"] = survey

    # ADM write out the map.
    write_atomically(outfile, pixmap, extname='PIXWEIGHTS', header=hdr)
    log.info('wrote map of HEALPixel weights to {}...t={:.1f}s'.format(
        outfile, time()-start))


def ident_for_randoms(nrandoms, filename):
    """Get the unique identifier string for each row in a random catalog.

    Parameters
    ----------
    nrandoms : :class:`int`
        Number of rows in a random catalog.
    filename : :class:`str`
        The filename of a random catalog.

    Returns
    -------
    :class:`~numpy.ndarray`
        Structured array with one column "IDENT" that is `nrandoms` long.

    Notes
    -----
    - Randoms are typically demarcated by the phrase randoms-ISEED-ISPLIT
      in the `filename`. The ISEED-ISPLIT populates the column "IDENT".
    """
    # ADM set up the output array.
    dt = [('IDENT', '<U4')]
    done = np.zeros(nrandoms, dtype=dt)

    # ADM extract the part of the filename after "randoms-".
    ender = os.path.basename(filename).split("randoms-")[-1]
    # ADM extract the regex that looks like ISEED-ISPLIT.
    done['IDENT'] = re.findall("[0-9]{1,2}-[0-9]{1,2}", ender)[0]

    return done


def read_randoms(infiles, test=False):
    """Read a random catalog to use for constructing sky maps.

    Parameters
    ----------
    infiles : :class:`str` or `list`
        Filename of random catalog or list of filenames. Files must have
        columns 'RA', 'DEC', 'EBV', 'PSFDEPTH_W1/W2/G/R/Z', 'NOBS_G/R/Z',
        'GALDEPTH_G/R/Z', 'PSFSIZE_G/R/Z', 'MASKBITS' and have been
        generated at the same density. If `randoms` is a list files will
        be concatenated in list-order.
    test : :class:`bool`, optional, defaults to ``False``
        If ``True`` then only read the first 100,000 entries in each
        random catalog. Useful for testing the code.

    Returns
    -------
    :class:`~numpy.ndarray`
        The random catalog read or concatenated from `infiles`.
    :class:`FITSHDR`
        The header of the FINAL file read from `infiles`. If `infiles`
        is a list then the DENSITY keyword in the header is returned as
        the SUM of the DENSITY in each file header.
    :class:`~numpy.ndarray`
        Structured array with one column "IDENT" that has the same number
        of rows as the output random catalog.

    Notes
    -----
    - The header of each filename in `infiles` must include the keyword
      "DENSITY" to establish the density used to make the random catalog.
    - If a list of filenames is passed, then the associated catalogs must
      all have been generated at the same density.
    - Randoms are typically demarcated by the phrase randoms-ISEED-ISPLIT
      in the filename. The ISEED-ISPLIT is what is returned as the array
      with column IDENT (to help track provenance).
    """
    # ADM if we're testing, only read in a subset of randoms.
    rows = None
    if test:
        rows = np.arange(100000)

    # ADM if a filename was passed for the random catalog, read it in...
    if isinstance(infiles, str):
        log.info("Reading in random catalog...t = {:.1f}s".format(time()-start))
        # ADM also need to know the density of randoms in the catalog.
        randoms, hdr = fitsio.read(infiles, rows=rows, header=True)
        # ADM add the IDENTity of this random catalog.
        ident = ident_for_randoms(len(randoms), infiles)
    # ADM ...otherwise if a list was passed, concatenate the randoms in
    # ADM the list and check they were generated at the same density.
    elif isinstance(infiles, list):
        randomsall = []
        densall = []
        identall = []
        for fn in infiles:
            log.info("Reading random catalog {}...t = {:.1f}s".format(
                fn, time()-start))
            randoms, hdr = fitsio.read(fn, rows=rows, header=True)
            # ADM add the IDENTity of this random catalog.
            ident = ident_for_randoms(len(randoms), fn)
            # ADM concatenate the random catalogs.
            randomsall.append(randoms)
            identall.append(ident)
            densall.append(hdr["DENSITY"])
            # ADM check all of the densities are the same.
            if not len(set(densall)) == 1:
                msg = "Densities in random catalogs do not match."
                log.critical(msg)
                for r, d in zip(randoms, densall):
                    log.info("{}: {}".format(r, d))
                raise ValueError(msg)
        # ADM concatenate randoms and store density.
        randoms = np.concatenate(randomsall)
        ident = np.concatenate(identall)
        hdr["DENSITY"] = np.sum(densall)
    else:
        msg = "randoms must be passed as either a list or a string!"
        log.critical(msg)
        raise ValueError(msg)

    log.info("Read {} total randoms at density {}... t = {:.1f}s".format(
        len(randoms), hdr["DENSITY"], time()-start))

    return randoms, hdr, ident


def rancat_name_to_mask_name(rancatname, lssmapdir=None):
    """Convert a random catalog name to the corresponding mask filename.

    Parameters
    ----------
    rancatname : :class:`str`
        Name of, or full path to, a random catalog.
    lssmapdir : :class:`str`, optional, defaults to $LSS_MAP_DIR
        Location of the directory that hosts all of the sky maps. If
       `lssmapdir` is ``None`` (or not passed), $LSS_MAP_DIR is used.

    Returns
    -------
    :class:`str`
        The full path to the corresponding mask filename in the
        lssmapdir directory.
    """
    outfn = os.path.basename(rancatname).replace(".fits", "-skymapmask.fits")

    # ADM formally grab $LSS_MAP_DIR in case lssmapdir=None was passed.
    lssmapdir = get_lss_map_dir(lssmapdir=lssmapdir)

    return os.path.join(lssmapdir, "maskvalues", outfn)


def rancat_name_to_map_name(rancatname, lssmapdir=None):
    """Convert random catalog name to corresponding map values filename.

    Parameters
    ----------
    rancatname : :class:`str`
        Name of, or full path to, a random catalog.
    lssmapdir : :class:`str`, optional, defaults to $LSS_MAP_DIR
        Location of the directory that hosts all of the sky maps. If
       `lssmapdir` is ``None`` (or not passed), $LSS_MAP_DIR is used.

    Returns
    -------
    :class:`str`
        The full path to the corresponding map values filename in the
        lssmapdir directory.
    """
    outfn = os.path.basename(rancatname).replace(".fits", "-skymapvalues.fits")

    # ADM formally grab $LSS_MAP_DIR in case lssmapdir=None was passed.
    lssmapdir = get_lss_map_dir(lssmapdir=lssmapdir)

    return os.path.join(lssmapdir, "mapvalues", outfn)


def rancat_names_to_pixweight_name(rancatlist, lssmapdir=None):
    """Convert random catalog name to corresponding pixweight filename.

    Parameters
    ----------
    rancatlist : :class:`list`
        List of strings of names of, or full paths to, random catalogs.
    lssmapdir : :class:`str`, optional, defaults to $LSS_MAP_DIR
        Location of the directory that hosts all of the sky maps. If
       `lssmapdir` is ``None`` (or not passed), $LSS_MAP_DIR is used.

    Returns
    -------
    :class:`str`
        The full path to the corresponding pixweight filename in the
        pixweight_maps_all directory, which is expected to exist one
        directory below the lssmapdir directory.

    Notes
    -----
    - Assumes a standard form for the names of the random catalogs.
    - Output format for filename resembles pixweight-SEEDS-ITERS
      where SEEDS are the random catalog seeds and ITERS are the
      random catalog iterations.
    """
    # ADM recover the idents for all random catalogs in the input list.
    idents = [ident_for_randoms(1, fn)[0][0] for fn in rancatlist]

    # ADM split into the seed and iteration for each random catalog.
    seed, it = np.array([ident.split('-') for ident in idents], dtype='int').T

    # ADM combine the seeds and iterations into a single filename..
    seedstr = "".join(np.array(sorted(list(set(seed))), dtype="str"))
    itstr = "".join(np.array(sorted(list(set(it))), dtype="str"))
    outfn = "randoms-pixweight-{}-{}.fits".format(seedstr, itstr)

    # ADM formally grab $LSS_MAP_DIR in case lssmapdir=None was passed.
    lssmapdir = get_lss_map_dir(lssmapdir=lssmapdir)

    return os.path.join(os.path.dirname(lssmapdir), "pixweight_maps_all", outfn)


def parse_mask_check(mxdata, maskcheck, check=False):
    """Turn a MASKCHECK string into a conditional and apply it.

    Parameters
    ----------
    mxdata : :class:`~numpy.ndarray`
        Array containing values to be masked.
    maskcheck : :class:`str`
        Conditional to apply to array. Must be a 3-character string where
        the first 2 characters are a logical comparison and the final
        character is an integer, such as "> 1" or "<=1".
    check : :class:`bool`
        If ``True`` then just check is `maskcheck` is an allowed string
        and return (Nothing is returned).

    Returns
    -------
    :class:`~numpy.ndarray`
        Boolean array which is ``True`` for values where the `maskcheck`
        conditional IS met.

    Notes
    -----
    - This is a tad unwieldy but is likely much safer than allowing an
      eval() or exec() function.
    """
    allowed = [">=", "> ", "<=", "< ", "!=", "=="]
    # ADM check maskcheck is constructed properly.
    if maskcheck[:2] not in allowed:
        msg = "First two characters of MASKCHECK must be one of {} (|{}| passed)"
        log.critical(msg.format(allowed, maskcheck))
        raise ValueError(msg.format(allowed, maskcheck))

    if len(maskcheck) != 3:
        msg = "MASKCHECK must be a 3-character string (|{}| passed)"
        log.critical(msg.format(maskcheck))
        raise ValueError(msg.format(maskcheck))

    try:
        checknum = int(maskcheck[2])
    except ValueError:
        msg = "MASKCHECK must end in an integer (|{}| passed)"
        log.critical(msg.format(maskcheck))
        raise ValueError(msg.format(maskcheck))

    if check:
        return

    if maskcheck[:2] == ">=":
        return mxdata >= checknum
    elif maskcheck[:2] == "> ":
        return mxdata > checknum
    elif maskcheck[:2] == "<=":
        return mxdata <= checknum
    elif maskcheck[:2] == "< ":
        return mxdata < checknum
    elif maskcheck[:2] == "!=":
        return mxdata != checknum
    elif maskcheck[:2] == "==":
        return mxdata == checknum
    # ADM redundant, but worth keeping for future code development.
    else:
        msg = "Conditional in allowed list not included in elif statements!"
        log.critical(msg)
        raise ValueError(msg)


def read_sky_map(mapname, lssmapdir=None):
    """A generic function to read a sky map, regardless of map format.

    Parameters
    ----------
    mapname : :class:`str`
        Name of a map that appears in the `maparray` global array, above.
    lssmapdir : :class:`str`, optional, defaults to $LSS_MAP_DIR
        Location of the directory that hosts all of the sky maps. If
       `lssmapdir` is ``None`` (or not passed), $LSS_MAP_DIR is used.

    Returns
    -------
    :class:`~numpy.ndarray`
        The data read from the map.
    """
    # ADM formally grab $LSS_MAP_DIR in case lssmapdir=None was passed.
    lssmapdir = get_lss_map_dir(lssmapdir=lssmapdir)

    # ADM extract the relevant map from the name.
    try:
        pixmap = maparray[maparray["MAPNAME"] == mapname][0]
    except IndexError:
        msg = "{} is not a named map in the maparray".format(mapname)
        log.critical(msg)
        raise ValueError(msg)

    # ADM construct the filename for, and read, the relevant map.
    fn = os.path.join(lssmapdir, pixmap["SUBDIR"], pixmap["FILENAME"])

    # MMM obtain nside for the relevant map.
    nsidemap = pixmap['NSIDE']

    # ADM try a few generic ways to read all types of maps.
    # ADM some maps are 1-D and have no column names.
    if pixmap["COLNAME"] == "NONE-IMAGE":
        mapdata = fitsio.read(fn)
    # ADM some maps are ALM maps.
    elif pixmap["MAPTYPE"] == "ALMMAP":
        # MMM WARNING - Hardwired values of ellmin, ellmax
        ellmin = 3
        ellmax = 2048
        alms = hp.read_alm(fn)
        mapdata = get_map_from_alms(alms, ellmin, ellmax,
                                    nside_out=nsidemap, nside_in=nsidemap)

    # ADM these are the MAPTYPE cases of either PIXMAP or PIXMASK.
    else:
        mapdata = fitsio.read(fn, 1, columns=pixmap["COLNAME"])
        # ADM if we're dealing with a 2-D map, use hp.read_map.
        if len(mapdata.shape) > 1:
            colnames = fitsio.read(fn, rows=0).dtype.names
            w = np.where([pixmap["COLNAME"] in i for i in colnames])
            # MMM test what passes through this piece of code
            # mapdata = hp.read_map(fn, field=w[0][0])
            # print("HELLO2", len(mapdata), pixmap["NSIDE"])
            # msg = "Non-regular field, specified column name ({}) HERE for (2-D) map {}?"
            # log.critical(msg.format(pixmap["COLNAME"], mapname))
            # ADM guard against a common incorrect-column-name error.
            if len(w) == 0:
                msg = "is the specified column name ({}) wrong for (2-D) map {}?"
                log.critical(msg.format(pixmap["COLNAME"], mapname))
                raise ValueError(msg.format(pixmap["COLNAME"], mapname))
            else:
                mapdata = hp.read_map(fn, field=w[0][0], nest=pixmap["NESTED"])

    return mapdata


def make_stardens(nside=512, gaiadir=None, dr="dr2", outdir=None, write=True):
    """Make a stellar density map using Gaia.

    Parameters
    ----------
    nside : :class:`int`, optional, defaults to nside=512
        Resolution (HEALPixel NESTED nside) at which to build the map.
    gaiadir : :class:`str`, optional, defaults to $GAIA_DIR
        Location of the directory that hosts HEALPixel-split Gaia files.
        See the Notes, below. Must be passed if $GAIA_DIR is not set (or
        if $GAIA_DIR is ``None``).
    dr : :class:`str`, optional, defaults to "dr2"
        If `gaiadir` is NOT passed, `dr` is used to determine which Gaia
        Data Release to use at NERSC. `dr` also sets criteria for a Gaia
        point-source (via :func:`desitarget.gaiamatch.gaia_psflike()`).
    outdir : :class:`str`, optional, defaults to $LSS_MAP_DIR/stardens
        Location of the directory to write output files. Must be passed
        if $LSS_MAP_DIR is not set (or if $LSS_MAP_DIR is ``None``).
    write : :class:`bool`, optional, defaults to ``True``
        If ``True`` then also write the output to file.

    Notes
    -----
    - Uses Gaia to generate HEALPixel map of stellar density. If the
      parameter `gaiadir` is not passed then the environment variable
      $GAIA_DIR must be set. At NERSC, $GAIA_DIR typically points to
      /global/cfs/cdirs/desi/target/gaia_dr3 or
      /global/cfs/cdirs/desi/target/gaia_dr2.
    - Mostly stolen from :func:`desitarget.randoms.stellar_density()`.
    """
    # ADM If gaiadir was not passed, then check that the GAIA_DIR is set
    # ADM and retrieve it.
    if gaiadir is None:
        gaiadir = get_gaia_dir(dr=dr)

    # ADM default to an output directory of lssmapdir/stardens.
    if outdir is None:
        outdir = get_lss_map_dir()
        outdir = os.path.join(outdir, "stardens")
        log.info("Setting output directory to {}".format(outdir))

    # ADM check that all the needed directories are set.
    msg = "{} must be passed or {} must be set!"
    if outdir is None:
        raise_myerror(msg.format("outdir", "$LSS_MAP_DIR"))

    # ADM retrieve the HEALPixel-ized Gaia sub-directory.
    hpdir = os.path.join(gaiadir, 'healpix')

    if not os.path.exists(hpdir):
        msg = "The Gaia HEALPixel directory is set to {} which doesn't exist. "
        msg += "Is $GAIA_DIR set correctly? Or did you pass a bad gaiadir?"
        raise_myerror(msg.format(hpdir))
    else:
        log.info("Gaia HEALPixel directory is set to {}".format(hpdir))

    # ADM the gaia_psflike function is only set for "edr3," which should
    # ADM have the same criteria as "dr3". Switch to "edr3", if needed.
    psfdr = dr
    if psfdr == "dr3":
        psfdr = "edr3"
    log.info("Using point-source criteria for {}".format(psfdr))

    # ADM the number of pixels and the pixel area at nside.
    npix = hp.nside2npix(nside)
    pixarea = hp.nside2pixarea(nside, degrees=True)

    # ADM an output array of all possible HEALPixels at nside.
    pixout = np.zeros(npix, dtype='int32')

    # ADM find all of the Gaia files.
    filenames = sorted(glob(os.path.join(hpdir, '*fits')))

    # ADM read in each file, restricting to the criteria for point
    # ADM sources and storing in a HEALPixel map at resolution nside.
    nfiles = len(filenames)
    t0 = time()
    for nfile, filename in enumerate(filenames):
        if nfile % 1000 == 0 and nfile > 0:
            elapsed = time() - t0
            rate = nfile / elapsed
            log.info("{}/{} files; {:.1f} files/sec; {:.1f} total mins elapsed"
                     .format(nfile, nfiles, rate, elapsed/60.))

        # ADM save memory, speed up by only reading a subset of columns.
        gobjs = fitsio.read(
            filename,
            columns=['RA', 'DEC', 'PHOT_G_MEAN_MAG', 'ASTROMETRIC_EXCESS_NOISE']
        )

        # ADM restrict to subset of point sources.
        ra, dec = gobjs["RA"], gobjs["DEC"]
        gmag = gobjs["PHOT_G_MEAN_MAG"]
        aen = gobjs["ASTROMETRIC_EXCESS_NOISE"]
        pointlike = gaia_psflike(aen, gmag, dr=psfdr)

        # ADM calculate the HEALPixels for the point sources.
        theta, phi = np.radians(90-dec[pointlike]), np.radians(ra[pointlike])
        pixnums = hp.ang2pix(nside, theta, phi, nest=True)

        # ADM return the counts in each pixel number...
        pixnum, pixcnt = np.unique(pixnums, return_counts=True)
        # ADM...and populate the output array with the counts.
        pixout[pixnum] += pixcnt

    # ADM calculate the stellar density.
    sd = pixout/pixarea

    # ADM set up an output structure with the STARDENS column.
    npix = hp.nside2npix(nside)
    done = np.zeros(npix, dtype=[('STARDENS', '>f4')])
    done["STARDENS"] = sd

    # ADM write the results to file.
    if write:
        hdr = fitsio.FITSHDR()
        hdr["GAIADIR"] = gaiadir
        hdr["NSIDE"] = nside
        outfn = os.path.join(outdir, "stardens.fits")
        write_atomically(outfn, done, extname='STARDENS', header=hdr)

    return done


def generate_mask(rancatname, lssmapdir=None, outdir=None, write=True):
    """Generate a file of mask values and TARGETID for a random catalog.

    Parameters
    ----------
    rancatname : :class:`str`
        Full path to a random catalog.
    lssmapdir : :class:`str`, optional, defaults to $LSS_MAP_DIR
        Location of the directory that hosts all of the sky maps. If
       `lssmapdir` is ``None`` (or not passed), $LSS_MAP_DIR is used.
    outdir : :class:`str`, optional, defaults to $LSS_MAP_DIR
        Location of the directory to write output files.
    write : :class:`bool`, optional, defaults to ``True``
        If ``True`` then also write the output to file.

    Returns
    -------
    :class:`~numpy.ndarray`
        An array that contains TARGETID and SKYMAP_MASK columns. This is
        also written to lssmapdir/maskvalues/rancatname-skymapmask.fits
        if `write` is passed as ``True``.
    """
    # ADM formally grab $LSS_MAP_DIR in case lssmapdir=None was passed.
    lssmapdir = get_lss_map_dir(lssmapdir=lssmapdir)

    # ADM default to an output directory of lssmapdir.
    if outdir is None:
        outdir = lssmapdir

    # ADM read the random catalog.
    randoms, hdr, ident = read_randoms(rancatname)

    # ADM store the Galactic coordinates for the randoms.
    c = SkyCoord(randoms["RA"]*u.degree, randoms["DEC"]*u.degree)
    lgal, bgal = c.galactic.l.value, c.galactic.b.value

    # ADM set up an output array. I use "i8" here because (as of the time
    # ADM of writing) fitsio does not support I/O for "u8" (uint64).
    dt = [('SKYMAP_MASK', 'i8'), ('TARGETID', '>i8')]
    done = np.zeros(len(randoms), dtype=dt)
    done["TARGETID"] = randoms["TARGETID"]

    # ADM grab the output filename.
    outfn = rancat_name_to_mask_name(rancatname, lssmapdir=outdir)

    # ADM first generate the bits from the LS masks.
    outmx = ls_bitmask_for_randoms(randoms, ident, lssmapdir=lssmapdir)

    # ADM limit to just the maps that correspond to masks...
    mxarray = maparray[maparray["MAPTYPE"] == "PIXMASK"]
    # ADM ... and loop through them.
    for mx in mxarray:
        log.info("Working on mask {}...t={:.1f}s".format(
            mx["MAPNAME"], time()-start))

        mxdata = read_sky_map(mx["MAPNAME"], lssmapdir=lssmapdir)

        # ADM construct a True/False version of this mask
        # ADM and store it in the array "ismasked".
        ismasked = parse_mask_check(mxdata, mx["MASKCHECK"])

        # ADM look up the nside of the mask-map.
        nsidemx = mx["NSIDE"]

        # ADM the coordinates to use for this mask-map.
        c1, c2 = randoms["RA"], randoms["DEC"]
        # ADM if needed, use Galactic coordinates.
        if mx["GALACTIC"]:
            log.info("Using Galactic coordinates for {}".format(mx["MAPNAME"]))
            c1, c2 = lgal, bgal

        # ADM determine whether each of the randoms is masked for the
        # ADM mask-map scheme (i.e. nested or ring).
        theta, phi = np.radians(90-c2), np.radians(c1)
        pixnums = hp.ang2pix(nsidemx, theta, phi, nest=mx["NESTED"])
        randmx = ismasked[pixnums]

        # ADM now we know the mask name is sensible, add the bit-mask for
        # ADM randoms that need masked (randoms with ismasked==True).
        mxnom = mx["MAPNAME"][:-5].upper()
        # ADM I cast as 'i8', here, because (as of the time of
        # ADM writing) fitsio does not support I/O for "u8" (uint64).
        outmx |= np.array(randmx * skymap_mask[mxnom], dtype='i8')

    # ADM now we've looped over all masks, construct the final array...
    done["SKYMAP_MASK"] = outmx
    # ADM ...and write it to file.
    if write:
        write_atomically(outfn, done, extname='PIXMASK', header=hdr)

    return done


def raise_myerror(msg):
    """Convenience function to raise a ValueError with a message"""
    log.critical(msg)
    raise ValueError(msg)


# MMM test create pixweight files.
def aux_test_mask():
    """Convenience function for testing create_pixweight_file()"""
    testd = '/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve'
    testfn = 'randoms-1-0.fits'
    randomcatlist = [os.path.join(testd, testfn), os.path.join(testd, testfn)]
    fieldslist = ['GALDEPTH_G', 'HALPHA_ERROR', 'APFLUX_IVAR_R',
                  'WISEMASK_W2', 'CALIB_Z']
    masklist = [131072, ['MASKBITS', 'ARTIFACTS', 'ELG_GAIA', 'LRG_UNWISE',
                         'EBV_SGF14'], 131072, ['KAPPA_PLANCK'], 4063232]
    outfn = '/global/u1/m/manera/pixweight.fits'
    nside_out = 512
    # lssmapdir='/global/cfs/cdirs/desi/survey/catalogs/external_input_maps'
    create_pixweight_file(
        randomcatlist, fieldslist, masklist, nside_out=nside_out,
        lssmapdir=None, outfn=outfn, write=True)


def create_pixweight_file(randomcatlist, fieldslist, masklist, nside_out=512,
                          lssmapdir=None, outfn=None, write=True, reg=None):
    """
    Creates a pixweight file from randoms filtered by bitmasks.

    Parameters
    ----------
    randomcatlist : :class:`list`
        List of strings representing (full paths to) random catalogs.
    fieldslist : :class:`list`
        List of fields/columns to process.
    masklist : :class:`list`
        Masks associated with `fieldslist` fields/columns. Entries must
        be either an integer or a list of mask names (strings), e.g.:
        [131072, ['MASKBITS', 'ELG_GAIA'], ['KAPPA_PLANCK'], 4063232]
    nside_out : :class:`int`, optional, defaults to 512
        Resolution (HEALPix nside) at which to build the output (NESTED)
        pixweight map.
    lssmapdir : :class:`str`, optional, defaults to $LSS_MAP_DIR
        Location of the directory that hosts all of the sky maps. If
        `lssmapdir` is ``None`` (or not passed), $LSS_MAP_DIR is used.
    outfn : :class:`str`, optional, defaults to ``None``
        Output filename. If not passed, the output from
        :func:`rancat_names_to_pixweight_name()` is used.
    write : :class:`bool`, optional, defaults to ``True``
        If ``True`` then also write the output to file.
    reg : :class:`str`, optional, defaults to ``None``
        If 'N' or 'S' are chosen and PHOTSYS is in the randoms, the
        randoms are cut to PHOTSYS==reg

    Returns
    -------
    :class:`~numpy.ndarray`
        Pixweight array of the requested masked fields. This is also
        written to file if `write`=``True``.
    """
    # MMM formally grab $LSS_MAP_DIR in case lssmapdir=None was passed.
    lssmapdir = get_lss_map_dir(lssmapdir=lssmapdir)

    #  ---- format checks -----
    # MMM check inputs are lists of correct length.
    for listy, word in zip([randomcatlist, fieldslist, masklist],
                           ["file(s)", "fields", "mask(s)"]):
        if not isinstance(listy, list):
            raise_myerror("the input {} is not a list".format(word))
    if len(fieldslist) != len(masklist):
        raise_myerror("number of masks and input fields do not match")

    # MMM check passed catalog names and fields are strings.
    for nom in randomcatlist + fieldslist:
        if not isinstance(nom, str):
            msg = "file and field names must be strings ({} is not)".format(nom)
            raise_myerror(msg)

    # MMM Check there are no repeated field names
    repeats = set([x for x in fieldslist if fieldslist.count(x) > 1])
    if repeats != set():
        msg = "Please don't use repeated field names in field list. \n \
        If you do need this feature contact the developers. \n \
        You have repeated {} ".format(repeats)
        raise raise_myerror(msg)

    # MMM Determine output filename.
    if write and not outfn:
        outfn = rancat_names_to_pixweight_name(rancatlist, lssmapdir=lssmapdir)
        log.warning("output filename not passed, defaulting to {}".format(outfn))

    # ------------------
    # MMM create bitmasklist from (and check) masklist.
    try:
        bitmasklist = [skymap_mask.mask("|".join(i)) if isinstance(i, list)
                       else int(i) for i in masklist]
    except (ValueError, TypeError, KeyError):
        msg = "input maskbits list should comprise integers or lists of strings "
        msg += "(and mask names must be strings), e.g.:\n"
        msg += "[131072, ['MASKBITS', 'ELG_GAIA'], ['CALIB_R'], 4063232]"
        raise_myerror(msg)

    # MMM---------  create header for later ------
    # MMM document fields
    hdr = {field: bitmask for field, bitmask in zip(fieldslist, bitmasklist)}
    # ADM document the input random catalogs...
    hdr["INFILES"] = randomcatlist
    # ADM and the directory from which we read the LSS maps.
    hdr["LSSMAPDIR"] = lssmapdir

    # ------ get columns/dtypes for pixweight files

    # Check which set of files to use
    # ADM need chxhdr if I want to check random catalogs generated at same density.
    randomcat = randomcatlist[0]
    stdfield, chxhdr = fitsio.read(randomcat, rows=[0], header=True)
    maskcol = ['SKYMAP_MASK']

    if 'SKYMAP_MASK' in stdfield.dtype.names:

        # MMM ra, dec, all fields and masks should be in the same file
        # MMM for now won't check if they come from the same density ***
        randomswithallfields = True
        skyfield = np.array([], dtype=[])  # dtype is needed

    else:

        # MMM Reading just one line to get names of columns and header.
        skymapvaluescat = rancat_name_to_map_name(randomcat, lssmapdir=lssmapdir)
        skymapmaskcat = rancat_name_to_mask_name(randomcat, lssmapdir=lssmapdir)

        skyfield = fitsio.read(skymapvaluescat, rows=[0])

    # MMM check if there are no foreign or misspelled items in fieldlist.
    foreign = [fieldslist[i] for i, x in enumerate(fieldslist) if x
               not in list(stdfield.dtype.names) + list(skyfield.dtype.names)]
    if foreign:
        msg = "You have some wrong or misspelled items in the field list\n \
        They are {} \n".format(foreign)
        raise raise_myerror(msg)

    # MMM select unique columns by matching to field list
    # MMM (standard field, sky_image, and mask).
    stdfcol = list(set(fieldslist).intersection(stdfield.dtype.names))
    skyfcol = list(set(fieldslist).intersection(skyfield.dtype.names))

    # MMM sanity check on ra dec.
    if not {"RA", "DEC"}.issubset(set(stdfield.dtype.names)):
        raise_myerror("RA or DEC field not found in randoms")

    # ------------- pixweight counts and creation ----------
    # MMM create healpix rec arrays for output pixweight table.
    npix = hp.nside2npix(nside_out)
    counts = np.zeros(npix, dtype=[(field, '>i4') for field in fieldslist])
    wcounts = np.zeros(npix, dtype=[(field, '>f4') for field in fieldslist])

    # ADM useful to cast lists as arrays to facilitate boolean indexing.
    fieldsarray, bitmaskarray = np.array(fieldslist), np.array(bitmasklist)

    # MMM loop over sets of files.
    for randomcat in randomcatlist:
        #read from dvs_ro
        randomcat = randomcat.replace('global','dvs_ro')
        # MMM log file we are reading.
        log.info("Reading in random catalog {} and associated files...t = {:.1f}s"
                 .format(randomcat, time()-start))

        # MMM names of the sky-map field and mask values, if needed
        skymapvaluescat = rancat_name_to_map_name(randomcat, lssmapdir=lssmapdir)
        skymapmaskcat = rancat_name_to_mask_name(randomcat, lssmapdir=lssmapdir)

        # MMM read RA DEC and SKYMAP_MASK for each random.
        # ADM read ALL needed columns from randomcat here as a speed-up.
        ranvalues, ranhdr = fitsio.read(randomcat, columns=stdfcol+['RA', 'DEC'],
                                        header=True)

        log.info("Read in random catalog {} and associated files...t = {:.1f}s"
                 .format(randomcat, time()-start))

        # MMM read field values; only if need be.
        if skyfcol:
            skymapvalues = fitsio.read(skymapvaluescat, columns=skyfcol)
        else:
            skymapvalues = []

        regsel = np.ones(len(ranvalues),dtype='bool')
        if reg is not None:
            if skymapvalues != []:
                print('for region selection, all fields must be in the randoms, code exiting!!!')
                print('fields that were loaded with randoms are '+str(ranvalues.dtype.names))
                return 'ERROR'
            regcol = fitsio.read(randomcat,columns=['PHOTSYS'])
            regsel = regcol['PHOTSYS'] == reg
            oldlength = len(ranvalues)
            ranvalues = ranvalues[regsel]
            log.info('cut randoms to selected region, kept '+str(len(ranvalues))+' out of '+str(oldlength)) 
        if not randomswithallfields:
            skymapmask = fitsio.read(skymapmaskcat, columns=maskcol)
        else:
            skymapmask = np.zeros(len(ranvalues), dtype=[('SKYMAP_MASK', 'i8')])
            skymapmask["SKYMAP_MASK"] = fitsio.read(randomcat,
                                                    columns=['SKYMAP_MASK'])[regsel]

        # ADM check all random catalogs were generated at same density.
        # MMM I can only do this if not reading from user made randoms
        # MMM Also don't check targetids match (they should by construction).
        if not randomswithallfields:
            if ranhdr["DENSITY"] != chxhdr["DENSITY"]:
                raise_myerror("Random catalogs {} and {} made at different densities"
                              .format(randomcat, randomcatlist[0]))

        # MMM find nested HEALPixel in the passed nside for each random.
        theta, phi = np.radians(90-ranvalues['DEC']), np.radians(ranvalues['RA'])
        randpixnums = hp.ang2pix(nside_out, theta, phi, nest=True)

        # MMM if all bitmasks are same, no need to set mask every time.
        # MMM mask-in (i.e., list selected) randoms.
        need2setmask = True
        if bitmasklist.count(bitmasklist[0]) == len(bitmasklist):
            bitmask = bitmasklist[0]
            need2setmask = False
            maskin = (skymapmask['SKYMAP_MASK'] & bitmask) == 0
            # uniq, ii, cnt = np.unique(randpixnums[maskin], return_inverse=True,
            #                          return_counts=True)

        ############################
        # MMM ----- read all fields at once ----
        log.info("Determining counts for {}...t = {:.1f}s".format(
            randomcat, time()-start))
        for col, values in zip([stdfcol, skyfcol], [ranvalues, skymapvalues]):
            if len(col) > 0:
                # ADM limit to just the fields/bitmasks corresponding to col.
                jj = np.array([fld in col for fld in fieldslist])
                for field, bitmask in zip(fieldsarray[jj], bitmaskarray[jj]):
                    if need2setmask:
                        maskin = (skymapmask['SKYMAP_MASK'] & bitmask) == 0
                    #    uniq, ii, cnt = np.unique(
                    #        randpixnums[maskin], return_inverse=True,
                    #        return_counts=True)
                    masknan = values[field]*0 == 0
                    maskhpun = values[field] != hp.UNSEEN
                    uniq, ii, cnt = np.unique(
                        randpixnums[maskin & masknan & maskhpun],
                        return_inverse=True, return_counts=True)
                    wcnt = np.bincount(ii, values[field][maskin & masknan & maskhpun])
                    counts[field][uniq] += cnt
                    wcounts[field][uniq] += wcnt

    ##########################
    # MMM compute weighted means.
    # MMM healpix unseen pixel value is -1.6375e+30.
    for field in fieldslist:

        log.info("raw numbers {} {} {} {}".format(
            field, np.sum(wcounts[field]), np.sum(counts[field]),
            np.sum(wcounts[field])/np.sum(counts[field]))
        )

        ii = counts[field] > 0
        wcounts[field][ii] = wcounts[field][ii]/counts[field][ii]
        wcounts[counts[field] == 0][field] = hp.UNSEEN

        log.info("final numbers {} {} {}".format(
            field, np.mean(wcounts[ii][field]), np.mean(counts[ii][field])))

        # for ii in range(0,len(counts[field])):
        #    if counts[ii][field] > 0:
        #        wcounts[ii][field] = wcounts[ii][field]/counts[ii][field]
        #    else:
        #        wcounts[ii][field] =hp.UNSEEN
        #
        # ii = counts[field] > 0
        # below was not actually dividing by counts for some reason
        # print(field,np.sum(wcounts[ii][field]),np.sum(counts[ii][field]),np.sum(wcounts[ii][field])/np.sum(counts[ii][field]),np.mean(wcounts[ii][field]/counts[ii][field]))

    # MMM Write atomically (sanity check done before).
    if write:
        write_atomically(outfn, wcounts, extname='PIXWEIGHT', header=hdr)

    return wcounts

def create_pixweight_file_allinone(randomcatlist, fieldslist, masklist, nside_out=512,
                          lssmapdir=None, outfn=None, write=True, regl=None):
    """
    Creates a pixweight file from randoms filtered by bitmasks.
    This concatenates all of the random info before putting onto maps and expects a list for regl if it is not None
    Parameters
    ----------
    randomcatlist : :class:`list`
        List of strings representing (full paths to) random catalogs.
    fieldslist : :class:`list`
        List of fields/columns to process.
    masklist : :class:`list`
        Masks associated with `fieldslist` fields/columns. Entries must
        be either an integer or a list of mask names (strings), e.g.:
        [131072, ['MASKBITS', 'ELG_GAIA'], ['KAPPA_PLANCK'], 4063232]
    nside_out : :class:`int`, optional, defaults to 512
        Resolution (HEALPix nside) at which to build the output (NESTED)
        pixweight map.
    lssmapdir : :class:`str`, optional, defaults to $LSS_MAP_DIR
        Location of the directory that hosts all of the sky maps. If
        `lssmapdir` is ``None`` (or not passed), $LSS_MAP_DIR is used.
    outfn : :class:`str`, optional, defaults to ``None``
        Output filename. If not passed, the output from
        :func:`rancat_names_to_pixweight_name()` is used.
    write : :class:`bool`, optional, defaults to ``True``
        If ``True`` then also write the output to file.
    regl : :class:`str`, optional, defaults to ``None``
        Use regl=['N','S'] to make separate N and S maps

    Returns
    -------
    :class:`~numpy.ndarray`
        Pixweight array of the requested masked fields. This is also
        written to file if `write`=``True``.
    """
    # MMM formally grab $LSS_MAP_DIR in case lssmapdir=None was passed.
    lssmapdir = get_lss_map_dir(lssmapdir=lssmapdir)

    #  ---- format checks -----
    # MMM check inputs are lists of correct length.
    for listy, word in zip([randomcatlist, fieldslist, masklist],
                           ["file(s)", "fields", "mask(s)"]):
        if not isinstance(listy, list):
            raise_myerror("the input {} is not a list".format(word))
    if len(fieldslist) != len(masklist):
        raise_myerror("number of masks and input fields do not match")

    # MMM check passed catalog names and fields are strings.
    for nom in randomcatlist + fieldslist:
        if not isinstance(nom, str):
            msg = "file and field names must be strings ({} is not)".format(nom)
            raise_myerror(msg)

    # MMM Check there are no repeated field names
    repeats = set([x for x in fieldslist if fieldslist.count(x) > 1])
    if repeats != set():
        msg = "Please don't use repeated field names in field list. \n \
        If you do need this feature contact the developers. \n \
        You have repeated {} ".format(repeats)
        raise raise_myerror(msg)

    # MMM Determine output filename.
    if write and not outfn:
        outfn = rancat_names_to_pixweight_name(rancatlist, lssmapdir=lssmapdir)
        log.warning("output filename not passed, defaulting to {}".format(outfn))

    # ------------------
    # MMM create bitmasklist from (and check) masklist.
    try:
        bitmasklist = [skymap_mask.mask("|".join(i)) if isinstance(i, list)
                       else int(i) for i in masklist]
    except (ValueError, TypeError, KeyError):
        msg = "input maskbits list should comprise integers or lists of strings "
        msg += "(and mask names must be strings), e.g.:\n"
        msg += "[131072, ['MASKBITS', 'ELG_GAIA'], ['CALIB_R'], 4063232]"
        raise_myerror(msg)

    # MMM---------  create header for later ------
    # MMM document fields
    hdr = {field: bitmask for field, bitmask in zip(fieldslist, bitmasklist)}
    # ADM document the input random catalogs...
    hdr["INFILES"] = randomcatlist
    # ADM and the directory from which we read the LSS maps.
    hdr["LSSMAPDIR"] = lssmapdir

    # ------ get columns/dtypes for pixweight files

    # Check which set of files to use
    # ADM need chxhdr if I want to check random catalogs generated at same density.
    randomcat = randomcatlist[0]
    stdfield, chxhdr = fitsio.read(randomcat, rows=[0], header=True)
    maskcol = ['SKYMAP_MASK']

    if 'SKYMAP_MASK' in stdfield.dtype.names:

        # MMM ra, dec, all fields and masks should be in the same file
        # MMM for now won't check if they come from the same density ***
        randomswithallfields = True
        skyfield = np.array([], dtype=[])  # dtype is needed

    else:

        # MMM Reading just one line to get names of columns and header.
        skymapvaluescat = rancat_name_to_map_name(randomcat, lssmapdir=lssmapdir)
        skymapmaskcat = rancat_name_to_mask_name(randomcat, lssmapdir=lssmapdir)

        skyfield = fitsio.read(skymapvaluescat, rows=[0])

    # MMM check if there are no foreign or misspelled items in fieldlist.
    foreign = [fieldslist[i] for i, x in enumerate(fieldslist) if x
               not in list(stdfield.dtype.names) + list(skyfield.dtype.names)]
    if foreign:
        msg = "You have some wrong or misspelled items in the field list\n \
        They are {} \n".format(foreign)
        raise raise_myerror(msg)

    # MMM select unique columns by matching to field list
    # MMM (standard field, sky_image, and mask).
    stdfcol = list(set(fieldslist).intersection(stdfield.dtype.names))
    skyfcol = list(set(fieldslist).intersection(skyfield.dtype.names))

    # MMM sanity check on ra dec.
    if not {"RA", "DEC"}.issubset(set(stdfield.dtype.names)):
        raise_myerror("RA or DEC field not found in randoms")

    # ------------- pixweight counts and creation ----------
    # MMM create healpix rec arrays for output pixweight table.
    npix = hp.nside2npix(nside_out)
    counts = np.zeros(npix, dtype=[(field, '>i4') for field in fieldslist])
    wcounts = np.zeros(npix, dtype=[(field, '>f4') for field in fieldslist])

    # ADM useful to cast lists as arrays to facilitate boolean indexing.
    fieldsarray, bitmaskarray = np.array(fieldslist), np.array(bitmasklist)

    # MMM loop over sets of files.
    ranl = []
    for randomcat in randomcatlist:
        #read from dvs_ro
        randomcat = randomcat.replace('global','dvs_ro')
        # MMM log file we are reading.
        log.info("Reading in random catalog {} and associated files...t = {:.1f}s"
                 .format(randomcat, time()-start))

        # MMM names of the sky-map field and mask values, if needed
        skymapvaluescat = rancat_name_to_map_name(randomcat, lssmapdir=lssmapdir)

        # MMM read RA DEC and SKYMAP_MASK for each random.
        # ADM read ALL needed columns from randomcat here as a speed-up.
        ranvalues, ranhdr = fitsio.read(randomcat, columns=stdfcol+['RA', 'DEC','PHOTSYS','SKYMAP_MASK'],
                                        header=True)
        ranl.append(ranvalues)
        log.info("Read in random catalog {} and associated files...t = {:.1f}s"
                 .format(randomcat, time()-start))
        del ranvalues
    ranvalues = np.concatenate(ranl)
    # MMM find nested HEALPixel in the passed nside for each random.
    theta, phi = np.radians(90-ranvalues['DEC']), np.radians(ranvalues['RA'])
    randpixnums = hp.ang2pix(nside_out, theta, phi, nest=True)
    # MMM if all bitmasks are same, no need to set mask every time.
    # MMM mask-in (i.e., list selected) randoms.
    need2setmask = True
    if bitmasklist.count(bitmasklist[0]) == len(bitmasklist):
        bitmask = bitmasklist[0]
        need2setmask = False
        maskin = (ranvalues['SKYMAP_MASK'] & bitmask) == 0
            # uniq, ii, cnt = np.unique(randpixnums[maskin], return_inverse=True,
            #                          return_counts=True)

        ############################
       
    if regl is None:
        regl = [1]

        #regsel = np.ones(len(ranvalues),dtype='bool')
    for reg in regl:
        if reg != 1:    
            regsel = ranvalues['PHOTSYS'] == reg
            outfn_tot = outfn+reg+'.fits'  
        else:
            regsel = np.ones(len(ranvalues),dtype='bool')
            outfn_tot = outfn+'.fits'
        log.info('cut randoms to selected region, kept '+str(len(ranvalues[regsel]))+' out of '+str(len(ranvalues)))
        # MMM ----- read all fields at once ----
        log.info("Determining counts for {}...t = {:.1f}s".format(
            randomcat, time()-start))
        #for col, values in zip(stdfcol, ranvalues[regsel]):
        #    if len(col) > 0:
                # ADM limit to just the fields/bitmasks corresponding to col.
        #jj = np.array([for fld in fieldslist])
        for field, bitmask in zip(fieldsarray, bitmaskarray):
            if need2setmask:
                maskin = (ranvalues['SKYMAP_MASK'] & bitmask) == 0
            #    uniq, ii, cnt = np.unique(
            #        randpixnums[maskin], return_inverse=True,
            #        return_counts=True)
            masknan = ranvalues[regsel][field]*0 == 0
            maskhpun = ranvalues[regsel][field] != hp.UNSEEN
            uniq, ii, cnt = np.unique(
                randpixnums[regsel][maskin[regsel] & masknan & maskhpun],
                return_inverse=True, return_counts=True)
            wcnt = np.bincount(ii, ranvalues[regsel][field][maskin[regsel] & masknan & maskhpun])
            counts[field][uniq] += cnt
            wcounts[field][uniq] += wcnt

        ##########################
        # MMM compute weighted means.
        # MMM healpix unseen pixel value is -1.6375e+30.
        for field in fieldslist:
    
            log.info("raw numbers {} {} {} {}".format(
                field, np.sum(wcounts[field]), np.sum(counts[field]),
                np.sum(wcounts[field])/np.sum(counts[field]))
            )
    
            ii = counts[field] > 0
            wcounts[field][ii] = wcounts[field][ii]/counts[field][ii]
            wcounts[counts[field] == 0][field] = hp.UNSEEN
    
            log.info("final numbers {} {} {}".format(
                field, np.mean(wcounts[ii][field]), np.mean(counts[ii][field])))
    
            # for ii in range(0,len(counts[field])):
            #    if counts[ii][field] > 0:
            #        wcounts[ii][field] = wcounts[ii][field]/counts[ii][field]
            #    else:
            #        wcounts[ii][field] =hp.UNSEEN
            #
            # ii = counts[field] > 0
            # below was not actually dividing by counts for some reason
            # print(field,np.sum(wcounts[ii][field]),np.sum(counts[ii][field]),np.sum(wcounts[ii][field])/np.sum(counts[ii][field]),np.mean(wcounts[ii][field]/counts[ii][field]))
    
        # MMM Write atomically (sanity check done before).
        if write:
            write_atomically(outfn_tot, wcounts, extname='PIXWEIGHT', header=hdr)

    return True



def generate_map_values(rancatname, lssmapdir=None, outdir=None, write=True):
    """Generate a file of map values and TARGETID for a random catalog.

    Parameters
    ----------
    rancatname : :class:`str`
        Full path to a random catalog.
    lssmapdir : :class:`str`, optional, defaults to $LSS_MAP_DIR
        Location of the directory that hosts all of the sky maps. If
       `lssmapdir` is ``None`` (or not passed), $LSS_MAP_DIR is used.
    outdir : :class:`str`, optional, defaults to $LSS_MAP_DIR
        Location of the directory to write output files.
    write : :class:`bool`, optional, defaults to ``True``
        If ``True`` then also write the output to file.

    Returns
    -------
    :class:`~numpy.ndarray`
        An array that contains TARGETID and map-value columns. This is
        also written to lssmapdir/mapvalues/rancatname-skymapvalues.fits
        if `write` is passed as ``True``.
    """
    # ADM formally grab $LSS_MAP_DIR in case lssmapdir=None was passed.
    lssmapdir = get_lss_map_dir(lssmapdir=lssmapdir)

    # ADM default to an output directory of lssmapdir.
    if outdir is None:
        outdir = lssmapdir

    # ADM read the random catalog.
    randoms, hdr, ident = read_randoms(rancatname)

    # ADM store the Galactic coordinates for the randoms.
    c = SkyCoord(randoms["RA"]*u.degree, randoms["DEC"]*u.degree)
    lgal, bgal = c.galactic.l.value, c.galactic.b.value

    # ADM grab the output filename.
    outfn = rancat_name_to_map_name(rancatname, lssmapdir=outdir)

    # MMM limit to maps that are not masks.
    maps = maparray[(maparray["MAPTYPE"] == "PIXMAP") |
                    (maparray["MAPTYPE"] == "ALMMAP")]

    # ADM set up an initial output array. We'll modify the dtypes later.
    dt = [('TARGETID', '>i8')]
    dt += [(mapname, '>f4') for mapname in maps["MAPNAME"]]
    done = np.zeros(len(randoms), dtype=dt)
    done["TARGETID"] = randoms["TARGETID"]

    # ADM loop through the maps to find the values.
    for pixmap in maps:
        mapname = pixmap["MAPNAME"]
        log.info("Working on map {}...t={:.1f}s".format(mapname, time()-start))

        mapdata = read_sky_map(mapname, lssmapdir=lssmapdir)

        # ADM the coordinates to use for this map.
        c1, c2 = randoms["RA"], randoms["DEC"]
        if pixmap["GALACTIC"]:
            log.info("Using Galactic coordinates for {} map".format(mapname))
            c1, c2 = lgal, bgal

        # MMM obtain nside for the relevant map.
        nsidemap = pixmap['NSIDE']

        # ADM determine the map values for each of the randoms in the
        # ADM map scheme (i.e. nested or ring).
        theta, phi = np.radians(90-c2), np.radians(c1)
        pixnums = hp.ang2pix(nsidemap, theta, phi, nest=pixmap["NESTED"])

        # ADM alter the dtype of the output for this map, if needed.
        if done[mapname].dtype != mapdata.dtype:
            dtmod = [(nom, dtyp) if nom != mapname else (nom, mapdata.dtype.str)
                     for nom, dtyp in dt]
            done = done.astype(dtmod)
        # ADM add the map values for the randoms.
        done[mapname] = mapdata[pixnums]

    # ADM now we've looped over all maps, write the final array to file.
    if write:
        write_atomically(outfn, done, extname='PIXMAP', header=hdr)

    return done


# MMM map from alms.
def get_map_from_alms(alms, ellmin, ellmax, nside_in=512, nside_out=512):
    '''
    Create Healpix map fom healpix alms
    1) transform alm to map, with nside_res, only with ell <= ellmax
    2) if nside_out < nside_in, warns and degrades de map to nside_map
    Inputs:
       alms : array alm
       nside_out  : nside output map
       nside_in : nside of input and output map
       ellmin, ellmax : min and max ell-values of alm
    '''
    # MMM fill with zeros outside the ell range.
    # MMM might be optimized/pythonized.
    r = []
    for i in range(len(alms)):
        if i < ellmin:
            r.append(0.0)
        elif i > ellmax:
            r.append(0.0)
        else:
            r.append(1.0)
    alm2 = hp.almxfl(alms, r)  # multiply kap_al by r.
    ptest = hp.alm2map(alm2, nside_in)
    if(nside_out < nside_in):
        ptest = hp.ud_grade(ptest, nside_out)
    return ptest


def sample_map(mapname, randoms, lssmapdir=None, nside=512):
    """Sample a systematics map.

    Parameters
    ----------
    mapname : :class:`str`
        Name of a map that appears in the `maparray` global array, above.
    randoms : :class:`~numpy.ndarray`
        Random catalog, as made by, e.g. :func:`read_randoms()`.
    lssmapdir : :class:`str`, optional, defaults to $LSS_MAP_DIR
        Location of the directory that hosts all of the sky maps. If
       `lssmapdir` is ``None`` (or not passed), $LSS_MAP_DIR is used.
    nside : :class:`int`, optional, defaults to nside=512
        Resolution (HEALPix nside) at which to build the (NESTED) map.

    Returns
    -------
    :class:`~numpy.ndarray`
        Single-column array of the map values for the randoms in a NESTED
        HEALPixel map at the given nside. The name of the column in the
        output array is `mapname` in upper-case letters.
    """
    # ADM limit to just the map we are working with.
    pixmap = maparray[maparray["MAPNAME"] == mapname]

    if len(pixmap) != 1:
        # ADM check somebody didn't include two maps with the same name.
        if len(pixmap) > 1:
            msg = "There are TWO maps in maparray that have MAPNAME={}!"
        # ADM check there's an entry in maparray for the passed map name.
        elif len(pixmap) < 1:
            msg = "There are NO maps in maparray that have MAPNAME={}!"
        log.critical(msg.format(mapname))
        raise ValueError(msg.format(mapname))

    # ADM now we know for sure we have a 1-D map, we can enforce that.
    pixmap = pixmap[0]

    # ADM construct the filename for, and read, the relevant map.
    lssmapdir = get_lss_map_dir(lssmapdir=lssmapdir)
    fn = os.path.join(lssmapdir, pixmap["SUBDIR"], pixmap["FILENAME"])
    mapdata = fitsio.read(fn, columns=pixmap["COLNAME"])

    # MMM obtain nside directly (as opposed of using hp.npix2nside(len(mapdata))
    nsidemap = pixmap['NSIDE']

    # MMM obtain map through which to pass the randoms
    # MMM WARNING - Hardwired values of ellmin, ellmax
    if (pixmap['MAPTYPE'] == 'ALMMAP'):
        ellmin = 3
        ellmax = 2048
        lssmapdir = get_lss_map_dir(lssmapdir)
        fn = os.path.join(lssmapdir, pixmap["SUBDIR"], pixmap["FILENAME"])
        alms = hp.read_alm(fn)
        mapdata = get_map_from_alms(alms, ellmin, ellmax,
                                    nside_out=nsidemap, nside_in=nsidemap)

    # ADM if needed, convert the randoms to Galactic coordinates.
    c1, c2 = randoms["RA"], randoms["DEC"]
    if pixmap["GALACTIC"]:
        log.info("Using Galactic coordinates for {} map".format(mapname))
        c = SkyCoord(c1*u.degree, c2*u.degree)
        c1, c2 = c.galactic.l.value, c.galactic.b.value

    # ADM determine the map values for each of the randoms in the
    # ADM map scheme (i.e. nested or ring).
    theta, phi = np.radians(90-c2), np.radians(c1)
    pixnums = hp.ang2pix(nsidemap, theta, phi, nest=pixmap["NESTED"])
    randmapvals = mapdata[pixnums]

    # ADM find the nested HEALPixel in the passed nside for each random.
    theta, phi = np.radians(90-randoms["DEC"]), np.radians(randoms["RA"])
    randpixnums = hp.ang2pix(nside, theta, phi, nest=True)

    # ADM determine the mean in each HEALPixel, weighted by the randoms.
    uniq, ii, cnt = np.unique(randpixnums, return_inverse=1, return_counts=1)
    randmeans = np.bincount(ii, randmapvals)/cnt

    # ADM set up the output array.
    dt = [(mapname.upper(), mapdata.dtype.type)]
    npix = hp.nside2npix(nside)
    done = np.zeros(npix, dtype=[(mapname.upper(), mapdata.dtype.type)])
    # ADM The method to find the means will skip any missing pixels, so
    # ADM populate on uniq indices to retain the missing pixels as zeros.
    done[uniq] = randmeans

    return done


def pure_healpix_map_filename(outdir, nsideproc, hpxproc):
    """Standardize name to write output maps made by pure_healpix_map()

    Parameters
    ----------
    outdir : :class:`str`
        The root output directory, e.g. $SCRATCH.
    nsideproc : :class:`int`
        Resolution (nested HEALPix nside) at which map was made.
    hpxproc : :class:`int`
        Pixel number in which maps was made.

    Returns
    -------
    :class:`str`
        A standardized filename built from the input parameters.
    """
    outfile = f"skymap-nside-{nsideproc}-pixel-{hpxproc}.fits"

    return os.path.join(outdir, outfile)


def get_quantities_in_a_brick(ras, decs, brickname, drdir, dustdir=None):
    """Per-band quantities at locations in a Legacy Surveys brick.

    Parameters
    ----------
    ras : :class:`~numpy.ndarray`
        Array of Right Ascensions (degrees).
    decs : :class:`~numpy.ndarray`
        Array of Declinations (degrees).
    brickname : :class:`~numpy.array`
        Name of a brick in which to look up the RA/Dec locations, e.g.,
        '1351p320'. For any RA/Dec locations not in the brick, values
        of zero will be returned for all quntities.
    drdir : :class:`str`
        The root directory pointing to a Legacy Surveys Data Release
        e.g. /global/cfs/cdirs/cosmo/data/legacysurvey/dr9.
    dustdir : :class:`str`, optional, defaults to $DUST_DIR+'/maps'
        The root directory pointing to SFD dust maps. If not
        sent the code will try to use $DUST_DIR+'/maps' before failing.

    Returns
    -------
    :class:`~numpy.ndarray`
        A structured array with columns of Legacy Surveys quantities at
        the passed RA/Dec locations looked up in the passed brick.

    Notes
    -----
    - Based on :func:`desitarget.randoms.get_quantities_in_a_brick()`.
      More information is included in the docstring for that function.
    """
    # ADM only intended to work on one brick, so die for larger arrays.
    if not isinstance(brickname, str):
        log.fatal("Only one brick can be passed at a time!")
        raise ValueError

    # ADM retrieve the dictionary of quantities at each location.
    # ADM aprad=0 is a speed-up to skip calculating aperture fluxes.
    qdict = dr8_quantities_at_positions_in_a_brick(ras, decs, brickname, drdir,
                                                   aprad=0.)

    # ADM catch where a coadd directory is completely missing.
    if len(qdict) > 0:
        # ADM if 2 different camera combinations overlapped a brick
        # ADM then we also need to duplicate the ras, decs.
        if len(qdict['photsys']) == 2*len(ras):
            ras = np.concatenate([ras, ras])
            decs = np.concatenate([decs, decs])

    # ADM the structured array to output.
    dt = [('RELEASE', '>i2'), ('BRICKNAME', 'S8'), ('RA', '>f8'), ('DEC', 'f8'),
          ('NOBS_G', 'i2'), ('NOBS_R', 'i2'), ('NOBS_Z', 'i2'),
          ('GALDEPTH_G', 'f4'), ('GALDEPTH_R', 'f4'), ('GALDEPTH_Z', 'f4'),
          ('PSFDEPTH_G', 'f4'), ('PSFDEPTH_R', 'f4'), ('PSFDEPTH_Z', 'f4'),
          ('PSFDEPTH_W1', 'f4'), ('PSFDEPTH_W2', 'f4'),
          ('PSFSIZE_G', 'f4'), ('PSFSIZE_R', 'f4'), ('PSFSIZE_Z', 'f4'),
          ('EBV', 'f4'), ('PHOTSYS', '|S1')]
    qinfo = np.zeros(len(ras), dtype=dt)

    # ADM retrieve the E(B-V) values for each random point.
    ebv = get_dust(ras, decs, dustdir=dustdir)

    # ADM catch the case where a coadd directory was missing.
    if len(qdict) > 0:
        # ADM store each quantity of interest in the structured array
        # ADM remembering that the dictionary keys are lower-case text.
        cols = qdict.keys()
        for col in cols:
            if col.upper() in qinfo.dtype.names:
                qinfo[col.upper()] = qdict[col]

    # ADM add the RAs/Decs, SFD dust values and brick name.
    qinfo["RA"] = ras
    qinfo["DEC"] = decs
    qinfo["EBV"] = ebv
    qinfo["BRICKNAME"] = brickname

    return qinfo


def pure_healpix_map(drdir, nsideproc, hpxproc, nside=8192, numproc=1):
    """Build a skymap at HEALPixel centers (in nested scheme).

    Parameters
    ----------
    drdir : :class:`str`
        The root directory pointing to a Legacy Surveys Data Release
        e.g. /global/cfs/cdirs/cosmo/data/legacysurvey/dr9.
    nsideproc : :class:`int`
        Resolution (nested HEALPix nside) at which to process the map.
        To facilitate parallelization, results will only be calculated
        for HEALPixel `hpxproc` at nside `nsideproc`.
    hpxproc : :class:`int`
        Pixel number (nested HEALPixel) for which to process the map.
        To facilitate parallelization, results will only be calculated
        for HEALPixel `hpxproc` at nside `nsideproc`.
    nside : :class:`int`, optional, defaults to nside=8192
        Resolution (HEALPix nside) at which to build the (nested) map.
    numproc : :class:`int`, optional, defaults to 1
        The number of processes over which to parallelize.

    Returns
    -------
    :class:`~numpy.ndarray`
        A numpy structured array of (some) values from the Legacy Survey
        stacks and the `maparray` at the start of this module. Values
        are looked up at the HEALPixel centers at the passed `nside` and
        returned within the boundaries of the (nested) HEALPixel that
        corresponds to `nsideproc` and `hpxproc`.

    Notes
    -----
    - `nsideproc` cannot be larger than `nside` as it doesn't make sense
      to parallelize by grouping larger pixels into smaller pixels.
    """
    start = time()
    log.info(f"Starting...t={time()-start:.1f}s")

    # ADM see Notes.
    if nsideproc > nside:
        msg = f"nsideproc (={nsideproc}) cannot be larger than nside (={nside})"
        raise ValueError(msg)
        log.critical(msg)

    # ADM recover all the nested HEALPixel centers at the passed nside.
    # ADM first calculate all the HEALPixel numbers at the passed nside.
    npix = hp.nside2npix(nside)
    hpx = np.arange(npix)
    # ADM now calculate the HEALPixel numbers at the processing nside.
    fac = (nside//nsideproc)**2
    largerhpx = hpx//fac
    # ADM reduce to just HEALPixels within the processing nside.
    ii = largerhpx == hpxproc
    hpx = hpx[ii]
    # ADM finally, recover the locations of the centers.
    ras, decs = hp.pix2ang(nside, hpx, nest=True, lonlat=True)
    log.info(f"Found {len(hpx)} HEALPixels at nside={nside} within processing"
             f" pixel={hpxproc} at nside={nsideproc}...t={time()-start:.1f}s")

    # ADM build the Legacy Surveys bricks object.
    bricks = brick.Bricks(bricksize=0.25)

    # ADM recover all the brick names at the HEALPixel centers.
    allbricknames = bricks.brickname(ras, decs)

    # ADM group the locations (dic values) by brick (dic keys).
    dbrick = {bn: [] for bn in allbricknames}
    for bn, ra, dec in zip(allbricknames, ras, decs):
        dbrick[bn].append([ra, dec])

    # ADM just the unique brick names.
    bricknames = list(dbrick.keys())
    nbricks = len(bricknames)
    log.info(f"{len(hpx)} HEALPixels are spread over {nbricks}"
             f" unique bricks...t={time()-start:.1f}s")

    # ADM calculate quantities in a brick, parallelizing by brick names.
    # ADM the critical function to run on every brick.
    def _get_quantities(brickname):
        """wrap dr8_quantities_at_positions_a_brick() for a brick name"""
        # ADM extract the brick locations from the brick dictionary.
        ras, decs = np.array(dbrick[brickname]).T
        randoms = get_quantities_in_a_brick(ras, decs, brickname, drdir)
        return randoms

    # ADM this is just to count bricks in _update_status.
    nbrick = np.zeros((), dtype='i8')
    t0 = time()
    # ADM write a total of 25 output messages during processing.
    interval = nbricks // 25
    # ADM catch the case of very small numbers of bricks.
    if interval == 0:
        interval = 1

    def _update_status(result):
        """wrapper function for the critical reduction operation,
           that occurs on the main parallel process"""
        if nbrick % interval == 0:
            elapsed = time() - t0
            rate = (nbrick + 1) / elapsed
            log.info(f"Processed {nbrick+1}/{nbricks} bricks; {rate:.1f} "
                     f"bricks/sec; {elapsed/60.:.1f} total mins elapsed")
            # ADM warn the user if code might exceed 4 hours.
            if nbrick > 0 and nbricks/rate > 4*3600.:
                msg = "May take > 4 hours to run, so fail on interactive nodes."
                log.warning(msg)

        nbrick[...] += 1    # this is an in-place modification.
        return result

    # - Parallel process input files.
    if numproc > 1:
        pool = sharedmem.MapReduce(np=numproc)
        with pool:
            qinfo = pool.map(_get_quantities, bricknames, reduce=_update_status)
    else:
        qinfo = list()
        for brickname in bricknames:
            qinfo.append(_update_status(_get_quantities(brickname)))

    qinfo = np.concatenate(qinfo)

    # ADM resolve north/south bricks, first removing bricks that are
    # ADM outside of the imaging footprint.
    inside = qinfo["RELEASE"] != 0
    resolved = np.concatenate([qinfo[~inside], resolve(qinfo[inside])])

    # ADM build the output array.
    donedt = [('HPXNUM', '>i8')] + resolved.dtype.descr
    done = np.zeros(len(resolved), dtype=donedt)
    for col in resolved.dtype.names:
        done[col] = resolved[col]

    # ADM add the HEALPixel number to the output array.
    done["HPXNUM"] = hp.ang2pix(nside, done["RA"], done["DEC"],
                                nest=True, lonlat=True)

    # ADM sort on HEALPixel before returning.
    ii = np.argsort(done["HPXNUM"])
    done = done[ii]

    log.info(f"Done...t={time()-start:.1f}s")

    return done


# ADM always start by running sanity checks on maparray.
sanity_check_map_array()
