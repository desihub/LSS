# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
LSS.sky_maps
============

Routines for building weight maps from randoms, etc., for systematics
"""
import os
import fitsio
import numpy as np
from time import time

from desitarget.randoms import pixmap

# ADM the DESI default logger.
from desiutil.log import get_logger

# ADM initialize the DESI default logger.
log = get_logger()

# ADM start the clock.
start = time()


# ADM note to remove: It's convenient to have an environment variable
# ADM from which to read external maps and write new maps. We can
# ADM remove this note once we decide where that should be at NERSC.
def get_lss_dir(lssdir=None):
    """Convenience function to grab the $LSS_DIR environment variable.

    Parameters
    ----------
    lssdir : :class:`str`, optional, defaults to $LSS_DIR
        If `lssdir` is passed, it is returned from this function. If it's
        not passed, the $LSS_DIR environment variable is returned.

    Returns
    -------
    :class:`str`
        If `lssdir` is passed, it is returned from this function. If it's
        not passed, the directory stored in the $LSS_DIR environment
        variable is returned.
    """
    if lssdir is None:
        lssdir = os.environ.get('LSS_DIR')
        # ADM check that the $LSS_DIR environment variable is set.
        if lssdir is None:
            msg = "Pass mtldir or set $LSS_DIR environment variable!"
            log.critical(msg)
            raise ValueError(msg)

    return lssdir


def wrap_pixmap(randoms, targets, nside=512, gaialoc=None, test=False):
    """HEALPix map from randoms (wrapper on desitarget.randoms.pixmap)

    Parameters
    ----------
    randoms : :class:`str` or :class:`list`
        Filename of random catalog or list of filenames. Files must have
        columns 'RA', 'DEC', 'EBV', 'PSFDEPTH_W1/W2/G/R/Z', 'NOBS_G/R/Z'
        'GALDEPTH_G/R/Z', 'PSFSIZE_G/R/Z', 'MASKBITS' and have been
        generated at the same density. If `randoms` is a list files will
        be concatenated in list-order.
    targets : :class:`~numpy.ndarray` or `str`
        Corresponding (i.e. same Data Release) file of targets, or the
        name of a directory containing HEALPixel-split targets that
        can be read by :func:`desitarget.io.read_targets_in_box()`.
    nside : :class:`int`, optional, defaults to nside=256
        Resolution (HEALPix nside) at which to build the (NESTED) map.
    gaialoc : :class:`str`, optional, defaults to ``None``
        Name of a FITS file that already contains a column "STARDENS",
        which is simply read in. If ``None``, the stellar density is
        constructed from files in $GAIA_DIR.
    test : :class:`bool`, optional, defaults to ``False``
        If ``True`` then only read the first 100,000 entries in each
        random catalog. Useful for testing the code.

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

    Notes
    -----
    - If `gaialoc` is ``None`` then $GAIA_DIR must be set.
    - Docstring mostly stolen from :func:`desitarget.randoms.pixmap()`.
    """
    # ADM if we're testing, only read in a subset of randoms.
    rows = None
    if test:
        rows = np.arange(100000)

    # ADM if a file name was passed for the random catalog, read it in...
    if isinstance(randoms, str):
        log.info("Reading in random catalog...t = {:.1f}s".format(time()-start))
        # ADM also need to know the density of randoms in the catalog.
        dens = fitsio.read_header(randoms, "RANDOMS")["DENSITY"]
        randoms = fitsio.read(randoms, rows=rows)
    # ADM ...otherwise if a list was passed, concatenate the randoms in
    # ADM the list and check they were generated at the same density.
    elif isinstance(randoms, list):
        randomsall = []
        densall = []
        for fn in randoms:
            log.info("Reading random catalog {}...t = {:.1f}s".format(
                fn, time()-start))
            randomsall.append(fitsio.read(fn, rows=rows))
            # ADM also need to know the density of randoms in the catalog.
            densall.append(fitsio.read_header(fn, "RANDOMS")["DENSITY"])
            # ADM check all of the densities are the same.
            if not len(set(densall)) == 1:
                msg = "Densities in random catalogs do not match."
                log.critical(msg)
                for r, d in zip(randoms, densall):
                    log.info("{}: {}".format(r, d))
                raise ValueError(msg)
        # ADM concatenate randoms and store density.
        randoms = np.concatenate(randomsall)
        dens = np.sum(densall)
    else:
        msg = "randoms must be passed as either a list or a string!"
        log.critical(msg)
        raise ValueError

    log.info("Read {} total randoms at density {}".format(len(randoms), dens))

    skymap, survey = pixmap(randoms, targets, dens, nside=nside, gaialoc=gaialoc)

    return skymap
