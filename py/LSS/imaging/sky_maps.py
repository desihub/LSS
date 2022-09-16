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
    ('COLNAME', 'O'), ('MAPTYPE', 'O'), ('MASKVALUE', 'i'), ('NESTED', '?'),
    ('GALACTIC', '?')
]
# ADM update with new maps in this global array on a new line.
# ADM 'NESTED' is True for nested HEALPix and False for the ring scheme.
# ADM `GALACTIC' is True for a map in Galactic coords, False for RA/Dec.
# -------
# MMM 'MaskValue' is the value of the masked bits;
# MMM if Maskvalue is nonzero degraded pixel values may have contributions from masked areas.
# MMM if Maskvalue does not exist (no formally masked pixels), it defaults to zero for PIXMAP
# MMM WARNING: not all columns have the same dtype and range
# MMM WARNING: Need to decide on colum with MaskValue
# MMM Provisional notes:
#     Halpha -   Temperature and error: any value f4, mask: u1 1 (ok) or above (masked)
#     Calibration - no column names as it is not a table but an image fits file; mask data == 0 same file.
#     EBVdustGaia - colum names Recon_Mean, Recon_Variance, Recon_VarianceCorr f8 ; no masks
#     EBV_SFG     - colum names: ebv f4, mask: status i4 (masked area if status > 0)
#     kappa       - values: index, real, imag;  mask :
maparray = np.array([
    ('HALPHA',     'Halpha', 'Halpha_fwhm06_0512.fits',       512, 'PIXMAP', 'TEMPERATURE', 0, False, True),
    ('HALPHAErr',  'Halpha', 'Halpha_error_fwhm06_0512.fits', 512, 'PIXMAP', 'ERROR', 0, False, True),
    ('HALPHA_MASK', 'Halpha', 'Halpha_mask_fwhm06_0512.fits', 512, 'PIXMASK', 'MASK', 1,  False, True),
    ('CALIBG',     'calibration', 'decam-ps1-0128-g.fits',    512, 'PIXMAP', 'NONE-IMAGE', 0, False, False),
    ('CALIBR',     'calibration', 'decam-ps1-0128-r.fits',    512, 'PIXMAP', 'NONE-IMAGE', 0, False, False),
    ('CALIBZ',     'calibration', 'decam-ps1-0128-z.fits',    512, 'PIXMAP', 'NONE-IMAGE', 0, False, False),
    ('CALIBMASKG', 'calibration', 'decam-ps1-0128-g.fits',    512, 'PIXMASK', 'NONE-IMAGE', 0, False, False),
    ('CALIBMASKR', 'calibration', 'decam-ps1-0128-r.fits',    512, 'PIXMASK', 'NONE-IMAGE', 0, False, False),
    ('CALIBMASKZ', 'calibration', 'decam-ps1-0128-z.fits',    512, 'PIXMASK', 'NONE-IMAGE', 0, False, False),
    ('EBV_GAIA_FW15',      'EBV', 'recon_fw15.fits',         2048, 'PIXMASK', 'Recon_mean', 0, False, True),
    ('EBV_GAIA_FW6p1',     'EBV', 'recon_fw6-1.fits',        2048, 'PIXMASK', 'Recon_mean', 0, False, True),
    ('EBV_SGF14',          'EBV', 'ps1-ebv-4.5kpc.fits',      512, 'PIXMAP', 'ebv',        0, False, True),
    ('EBV_SGF14_mask',     'EBV', 'ps1-ebv-4.5kpc.fits',      512, 'PIXMASK', 'EBVcolumn', 0, False, True),
    ('KAPPA_PLANCK',     'kappa', 'dat_klm.fits',            2048, 'ALMMAP', 'NONE-3col',  0, False, True),
    ('KAPPA_PLANCK_MASK', 'kappa', 'mask.fits.gz',           2048, 'PIXMASK', 'I',        0, False, True),
    ], dtype=mapdt)


def get_lss_map_dir(lssmapdir=None):
    """Convenience function to get the $LSS_MAP_DIR environment variable.

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
    """Grab the LSS directory from the $LSS_MAP_DIR environment variable.

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
    lssmapdir = get_lss_map_dir(lssmapdir)

    # ADM convert between upper/lower case as needed.
    if "s" in survey:
        surv = survey.upper()
    else:
        surv = survey.lower()

    if surv not in ["SV1", "SV2", "SV3", "main"]:
        msg = "Survey: {} not recognized".format(survey)
        log.critical(msg)
        raise ValueError(msg)

    return os.path.join(os.path.dirname(get_lss_map_dir(lssmapdir)), surv, "LSS")


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
    lssmapdir = get_lss_map_dir(lssmapdir)

    # ADM make sure the tracer names are upper-case
    tracer = tracer.upper()

    if tracer not in ["LRG", "ELG"]:
        msg = "Tracer: {} not recognized".format(tracer)
        log.critical(msg)
        raise ValueError(msg)

    version = "v1"
    if tracer == "LRG":
        version = "v1.1"

    basedir = (os.path.dirname(get_lss_map_dir(lssmapdir)))

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

    # ADM set up an output array.
    done = np.zeros(len(randoms), dtype="u8")

    # ADM need to process once per different ident/random catalog.
    sidents = sorted(list(set(ident["IDENT"])))
    for sident in sidents:
        log.info("Working on random catalog {}".format(sident))
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
                trans = np.zeros(len(lsmx), dtype="u8")
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
            log.info("Transposing bits from {} to skymap_mask".format(mxname))
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
                # ADM be careful, here, to cast as an unsigned integer.
                trans |= np.array(is_bit_set * skymap_mask[snom], dtype='u8')
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

    pixmap, survey = wrap_pixmap(randoms, targets, nside=nside, gaialoc=gaialoc,
                                 test=test)

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
       `lssmapdir` is ``None`` (or not passed), $LSS_MAxP_DIR is used.
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
    lssmapdir = get_lss_map_dir(lssmapdir)
    fn = os.path.join(lssmapdir, pixmap["SUBDIR"], pixmap["FILENAME"])
    mapdata = fitsio.read(fn, columns=pixmap["COLNAME"])

    # ADM derive the nside of the map from its length.
    nsidemap = hp.npix2nside(len(mapdata))

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
