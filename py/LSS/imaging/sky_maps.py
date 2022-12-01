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

from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.wcs import WCS

from desitarget.io import read_targets_header
from desitarget.geomask import match

# ADM the DESI default logger.
from desiutil.log import get_logger

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
    ('HALPHA',     'Halpha', 'Halpha_fwhm06_0512.fits',          512, 'PIXMAP',  'TEMPERATURE',   '', False, True),
    ('HALPHA_ERROR',  'Halpha', 'Halpha_error_fwhm06_0512.fits', 512, 'PIXMAP',  'ERROR',         '', False, True),
    ('HALPHA_MASK', 'Halpha', 'Halpha_mask_fwhm06_0512.fits',    512, 'PIXMASK', 'MASK',       '> 1', False, True),
    ('CALIB_G',     'calibration', 'decam-ps1-0128-g.fits',      128, 'PIXMAP',  'NONE-IMAGE',    '', False, False),
    ('CALIB_R',     'calibration', 'decam-ps1-0128-r.fits',      128, 'PIXMAP',  'NONE-IMAGE',    '', False, False),
    ('CALIB_Z',     'calibration', 'decam-ps1-0128-z.fits',      128, 'PIXMAP',  'NONE-IMAGE',    '', False, False),
    ('CALIB_G_MASK', 'calibration', 'decam-ps1-0128-g.fits',     128, 'PIXMASK', 'NONE-IMAGE', '==0', False, False),
    ('CALIB_R_MASK', 'calibration', 'decam-ps1-0128-r.fits',     128, 'PIXMASK', 'NONE-IMAGE', '==0', False, False),
    ('CALIB_Z_MASK', 'calibration', 'decam-ps1-0128-z.fits',     128, 'PIXMASK', 'NONE-IMAGE', '==0', False, False),
    ('EBV_GAIA_FW15',      'EBV', 'recon_fw15.fits',            2048, 'PIXMAP',  'Recon_Mean',    '', False, True),
    ('EBV_GAIA_FW6P1',     'EBV', 'recon_fw6-1.fits',           2048, 'PIXMAP',  'Recon_Mean',    '', False, True),
    ('EBV_SGF14',          'EBV', 'ps1-ebv-4.5kpc.fits',         512, 'PIXMAP',  'ebv',           '', False, True),
    ('EBV_SGF14_MASK',     'EBV', 'ps1-ebv-4.5kpc.fits',         512, 'PIXMASK', 'status',     '< 0', False, True),
    ('KAPPA_PLANCK',     'kappa', 'dat_klm.fits',               2048, 'ALMMAP',  'NONE-3col',     '', False, True),
    ('KAPPA_PLANCK_MASK', 'kappa', 'mask.fits.gz',              2048, 'PIXMASK', 'I',          '==0', False, True),
    ], dtype=mapdt)


def sanity_check_map_array():
    """Convenience function to check the format of the map_array global.
    """
    log.info("Running sanity checks on maparray...")

    for skymap in maparray:

        mapname = skymap['MAPNAME']

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


def rancat_name_to_pixweight_name(rancatname, lssmapdir=None):
    """Convert random catalog name to corresponding pixweight filename.

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
        The full path to the corresponding pixweight filename in the
        pixweight_maps_all directory, which is expected to exist one
        directory below the lssmapdir directory.
    """
    outfn = os.path.basename(rancatname).replace(".fits", "-pixweight.fits")

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
    else:
        mapdata = fitsio.read(fn, columns=pixmap["COLNAME"])
        # ADM if we're dealing with a 2-D map, use hp.read_map.
        if len(mapdata.shape) > 1:
            colnames = fitsio.read(fn, rows=0).dtype.names
            w = np.where([pixmap["COLNAME"] in i for i in colnames])
            # ADM guard against a common incorrect-column-name error.
            if len(w) == 0:
                msg = "is the specified column name ({}) wrong for (2-D) map {}?"
                log.critical(msg.format(pixmap["COLNAME"]), mapname)
                raise ValueError(msg.format(pixmap["COLNAME"]), mapname)
            else:
                mapdata = hp.read_map(fn, field=w[0][0])

    return mapdata


def generate_mask(rancatname, lssmapdir=None, write=True):
    """Generate a file of mask values and TARGETID for a random catalog.

    Parameters
    ----------
    rancatname : :class:`str`
        Full path to a random catalog.
    lssmapdir : :class:`str`, optional, defaults to $LSS_MAP_DIR
        Location of the directory that hosts all of the sky maps. If
       `lssmapdir` is ``None`` (or not passed), $LSS_MAP_DIR is used.
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
    outfn = rancat_name_to_mask_name(rancatname, lssmapdir=lssmapdir)

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
            log.info("Using Galactic coordinates for {} map".format(mx["MAPNAME"]))
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
                          lssmapdir=None, outfn=None, write=True):
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
        :func:`rancat_name_to_pixweight_name()` is used.
    write : :class:`bool`, optional, defaults to ``True``
        If ``True`` then also write the output to file.

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
        outfn = rancat_name_to_pixweight_name(rancatname, lssmapdir=lssmapdir)
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

    ##########
    # MMM get columns/dtype from first file + associated skymap/skymask.
    # MMM Reading just one line to get names of columns and header.
    randomcat = randomcatlist[0]
    skymapvaluescat = rancat_name_to_map_name(randomcat, lssmapdir=lssmapdir)
    skymapmaskcat = rancat_name_to_mask_name(randomcat, lssmapdir=lssmapdir)

    # ADM need chxhdr to check random catalogs generated at same density.
    stdfield, chxhdr = fitsio.read(randomcat, rows=[0], header=True)
    skyfield = fitsio.read(skymapvaluescat, rows=[0])

    # MMM check if there are no foreign or misspelled items in fieldlist
    # MMM check all items in fieldlist belong to the stdfcol or skyfcol intersections

    foreign = [fieldslist[i] for i, x in enumerate(fieldslist)
               if x not in stdfcol+skyfcol]
    if foreign:
        msg = "You have some wrong or misspelled items in the field list\n \
        They are {} \n".format(foreign)
        raise raise_myerror(msg)

    # MMM select unique columns by matching to field list
    # MMM (standard field, sky_image, and mask).
    stdfcol = list(set(fieldslist).intersection(stdfield.dtype.names))
    skyfcol = list(set(fieldslist).intersection(skyfield.dtype.names))
    maskcol = ['SKYMAP_MASK']

    # MMM sanity check on ra dec.
    if not {"RA", "DEC"}.issubset(set(stdfield.dtype.names)):
        raise_myerror("RA or DEC field not found in randoms")

    ##########
    # MMM create header for later.
    hdr = {field: bitmask for field, bitmask in zip(fieldslist, bitmasklist)}
    # ADM should document the input random catalogs...
    hdr["INFILES"] = randomcatlist
    # ADM ...and the directory from which we read the LSS maps.
    hdr["LSSMAPDIR"] = lssmapdir

    # MMM create healpix rec arrays for output pixweight table.
    npix = hp.nside2npix(nside_out)
    counts = np.zeros(npix, dtype=[(field, '>i4') for field in fieldslist])
    wcounts = np.zeros(npix, dtype=[(field, '>f4') for field in fieldslist])

    # MMM loop over sets of files.
    for randomcat in randomcatlist:
        # MMM log file we are reading
        log.info("Reading in random catalog {} and associated files...t = {:.1f}s"
                 .format(randomcat, time()-start))

        # MMM names of the sky-map field and mask values.
        skymapvaluescat = rancat_name_to_map_name(randomcat, lssmapdir=lssmapdir)
        skymapmaskcat = rancat_name_to_mask_name(randomcat, lssmapdir=lssmapdir)

        # MMM read RA DEC and SKYMAP_MASK for each random.
        # ADM read ALL needed columns from randomcat here as a speed-up.
        ranvalues, ranhdr = fitsio.read(randomcat, columns=stdfcol+['RA', 'DEC'],
                                        header=True)

        # MMM read field values; only if need be.
        if skyfcol:
            skymapvalues = fitsio.read(skymapvaluescat, columns=skyfcol)
        else:
            skymapvalues = []

        skymapmask = fitsio.read(skymapmaskcat, columns=maskcol)

        # ADM check all random catalogs were generated at same density.
        if ranhdr["DENSITY"] != chxhdr["DENSITY"]:
            raise_myerror("Random catalogs {} and {} made at different densities"
                          .format(randomcat, randomcatlist[0]))

        # MMM Don't check targetids match (they should by construction).
        # MMM find nested HEALPixel in the passed nside for each random.
        theta, phi = np.radians(90-ranvalues['DEC']), np.radians(ranvalues['RA'])
        randpixnums = hp.ang2pix(nside_out, theta, phi, nest=True)

        # MMM if all bitmasks are same, no need to set mask every time.
        # MMM mask-in (i.e., list selected) randoms.
        need2setmask = True
        if bitmasklist.count(bitmasklist[0]) == len(bitmasklist):
            need2setmask = False
            maskin = (skymapmask['SKYMAP_MASK'] & bitmask) == 0
            uniq, ii, cnt = np.unique(randpixnums[maskin], return_inverse=True,
                                      return_counts=True)

        ############################
        # MMM ----- read all fields at once ----
        for col, values in zip([stdfcol, skyfcol], [ranvalues, skymapvalues]):
            if len(col) > 0:
                for field, bitmask in zip(fieldslist, bitmasklist):
                    if need2setmask:
                        maskin = (skymapmask['SKYMAP_MASK'] & bitmask) == 0
                        uniq, ii, cnt = np.unique(
                            randpixnums[maskin], return_inverse=True,
                            return_counts=True)
                    wcnt = np.bincount(ii, values[field][maskin])
                    counts[field][uniq] += cnt
                    wcounts[field][uniq] += wcnt

    ##########################
    # MMM compute weighted means.
    # MMM healpix unseen pixel value is -1.6375e+30.
    for field in fieldslist:
        ii = counts[field] > 0
        wcounts[ii][field] = wcounts[ii][field] / counts[ii][field]
        wcounts[counts[field] == 0][field] = hp.UNSEEN

    # MMM Write atomically (sanity check done before).
    if write:
        write_atomically(outfn, wcounts, extname='PIXWEIGHT', header=hdr)

    return wcounts


def generate_map_values(rancatname, lssmapdir=None, write=True):
    """Generate a file of map values and TARGETID for a random catalog.

    Parameters
    ----------
    rancatname : :class:`str`
        Full path to a random catalog.
    lssmapdir : :class:`str`, optional, defaults to $LSS_MAP_DIR
        Location of the directory that hosts all of the sky maps. If
       `lssmapdir` is ``None`` (or not passed), $LSS_MAP_DIR is used.
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

    # ADM read the random catalog.
    randoms, hdr, ident = read_randoms(rancatname)

    # ADM store the Galactic coordinates for the randoms.
    c = SkyCoord(randoms["RA"]*u.degree, randoms["DEC"]*u.degree)
    lgal, bgal = c.galactic.l.value, c.galactic.b.value

    # ADM grab the output filename.
    outfn = rancat_name_to_map_name(rancatname, lssmapdir=lssmapdir)

    # ADM limit to just the maps that correspond to pixel-maps.
    maps = maparray[maparray["MAPTYPE"] == "PIXMAP"]

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
    # ADM limit to just the map are we working with.
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


# ADM always start by running sanity checks on maparray.
sanity_check_map_array()
