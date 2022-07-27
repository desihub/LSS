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

from desitarget.io import read_targets_header

# ADM the DESI default logger.
from desiutil.log import get_logger

# ADM initialize the DESI default logger.
log = get_logger()

# ADM start the clock.
start = time()


def get_lss_map_dir(lssmapdir=None):
    """Convenience function to get the $LSS_MAP_DIR environment variable.

    Parameters
    ----------
    lssdir : :class:`str`, optional, defaults to $LSS_MAP_DIR
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

    pixmap, survey = wrap_pixmap(randoms, targets, nside=nside, gaialoc=gaialoc, test=test)

    hdr["SURVEY"] = survey

    #ADM write out the map.
    write_atomically(outfile, pixmap, extname='PIXWEIGHTS', header=hdr)
    log.info('wrote map of HEALPixel weights to {}...t={:.1f}s'.format(
        outfile, time()-start))


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

    Notes
    -----
    - The header of each filename in `infiles` must include the keyword
      "DENSITY" to establish the density used to make the random catalog.
    - If a list of filenames is passed, then the associated catalogs must
      all have been generated at the same density.
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
    # ADM ...otherwise if a list was passed, concatenate the randoms in
    # ADM the list and check they were generated at the same density.
    elif isinstance(infiles, list):
        randomsall = []
        densall = []
        for fn in infiles:
            log.info("Reading random catalog {}...t = {:.1f}s".format(
                fn, time()-start))
            randoms, hdr = fitsio.read(fn, rows=rows, header=True)
            randomsall.append(randoms)
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
        hdr["DENSITY"] = np.sum(densall)
    else:
        msg = "randoms must be passed as either a list or a string!"
        log.critical(msg)
        raise ValueError

    log.info("Read {} total randoms at density {}".format(
        len(randoms), hdr["DENSITY"]))

    return randoms, hdr
