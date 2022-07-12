# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
LSS.pixel_nobs
==============

Routines for adding pixel-level nobs values to target files.
"""
import os
import fitsio
import numpy as np
from time import time
import astropy.io.fits as fits
from astropy.wcs import WCS
from glob import glob

from desitarget.randoms import dr_extension
from desitarget.internal import sharedmem
from desitarget.geomask import match_to
from desitarget.io import write_with_units
from desitarget.geomask import is_in_hp

# ADM the DESI default logger.
from desiutil.log import get_logger

# ADM initialize the DESI default logger.
log = get_logger()

# ADM start the clock.
start = time()


def nexp_at_positions_in_a_brick(ras, decs, brickname, nors, drdir):
    """NEXP (per-band) at positions in a Legacy Surveys brick.

    Parameters
    ----------
    ras : :class:`~numpy.array`
        Right Ascensions of interest (degrees).
    decs : :class:`~numpy.array`
        Declinations of interest (degrees).
    brickname : :class:`str`
        Name of brick which contains RA/Dec positions, e.g., '1351p320'.
    nors : :class:`str`
        Pass "north" for northern bricks, "south" for southern bricks.
    drdir : :class:`str`
        The root directory pointing to a Legacy Surveys Data Release,
        e.g. /global/project/projectdirs/cosmo/data/legacysurvey/dr8.

    Returns
    -------
    :class:`dictionary`
        The number of observations (`nobs_x`) at each position/brickname.

    Notes
    -----
    - Adapted from :func:`desitarget.quantities_at_positions_in_a_brick()`
    - h/t to Rongpu Zhou for noticing that the targets and randoms used
      slightly different definitions of NOBS.
    """
    # ADM check nors is one of the allowed strings.
    if nors not in ["north", "south"]:
        msg = 'nors not "north"/"south" for brick '.format(brickname)
        log.critical(msg)
        raise ValueError(msg)

    # ADM expand the drdir to include "north"/"south".
    drdir = os.path.join(drdir, nors)

    # ADM determine whether the coadd files have extension .gz or .fz
    # based on the DR directory.
    extn, extn_nb = dr_extension(drdir)

    # ADM the output array.
    dt = [('PIXEL_NOBS_G', 'i2'), ('PIXEL_NOBS_R', 'i2'), ('PIXEL_NOBS_Z', 'i2')]
    nexp = np.zeros(len(ras), dtype=dt)

    # as a speed up, assume images in different filters for the brick have
    # the same WCS -> if read once (iswcs=True), use this info.
    # ADM this approach isn't strictly faster, I just included it for
    # ADM consistency with how the randoms are generated.
    iswcs = False

    # ADM {} are brick name and filter name (g/r/z).
    rootdir = os.path.join(drdir, 'coadd', brickname[:3], brickname)
    fileform = os.path.join(rootdir, 'legacysurvey-{}-nexp-{}.fits.') + str(extn)

    # ADM loop through the filters and store the number of observations
    # ADM at the RA and Dec positions of the passed points.
    for filt in ['g', 'r', 'z']:
        col = 'PIXEL_NOBS_' + filt.upper()
        fn = fileform.format(brickname, filt)
        if os.path.exists(fn):
            img = fits.open(fn)[extn_nb]
            # ADM if we've yet to succeed read the wcs information.
            if not iswcs:
                w = WCS(img.header)
                x, y = w.all_world2pix(ras, decs, 0)
                iswcs = True
            nexp[col] = img.data[y.round().astype("int"),
                                 x.round().astype("int")]
        # ADM if the file doesn't exist, set quantities to zero.
        else:
            nexp[col] = np.zeros(len(ras), dtype='i2')

    return nexp


def make_nexp_for_target_file(targfile, drdir, numproc=1):
    """Look up pixel-level NOBS (from coadds/stacks) for one target file.

    Parameters
    ----------
    targfile : :class:`str`
        Full path to a target file.
    drdir : :class:`str`
        Root directory for a Data Release from the Legacy Surveys
        e.g. /global/project/projectdirs/cosmo/data/legacysurvey/dr9.
    numproc : :class:`int`, optional, defaults to 1
        The number of processes to parallelize across. The default is
        to run the code in serial.

    Returns
    -------
    :class:`~numpy.array`
        The targets in the input `targfile` in the original order with
        the standard quantities RA, DEC, BRICKNAME, TARGETID and PHOTSYS
        and added quantities PIXEL_NOBS_G/R/Z, which are the pixel-level
        number of observations in each band.

    Notes
    -----
    - Useful as the NOBS listed in the target files differs from the
      pixel-level NOBS assigned to the DESI random catalogs.
    """
    # ADM read in the needed target columns...
    targs = fitsio.read(targfile, columns=
                        ["RA", "DEC", "BRICKNAME", "TARGETID", "PHOTSYS"])
    # ADM ...and compile the list of brick names in the file.
    bricknames = list(set(targs["BRICKNAME"]))
    nbricks = len(bricknames)

    # ADM wrapper to facilitate parallelization.
    def _get_nexp(brickname):
        """wrapper on nexp_at_positions_in_a_brick() for a brick name"""
        # ADM extract the information for one brick.
        ii = targs["BRICKNAME"] == brickname
        brick = targs[ii]

        # ADM determine if this is a northern or southern brick.
        if set(brick["PHOTSYS"]) == {'S'}:
            nors = "south"
        elif set(brick["PHOTSYS"]) == {'N'}:
            nors = "north"
        else:
            msg = 'PHOTSYS not all "N" or "S" for brick: {}'.format(brickname)
            log.critical(msg)
            raise ValueError(msg)            

        # ADM call the actual code.
        nexp =  nexp_at_positions_in_a_brick(brick["RA"], brick["DEC"], brickname,
                                             nors, drdir)

        # ADM make a table of all of the required information...
        dt = brick.dtype.descr + nexp.dtype.descr
        done = np.zeros(len(brick), dtype=dt)
        for col in brick.dtype.names:
            done[col] = brick[col]
        for col in nexp.dtype.names:
            done[col] = nexp[col]

        return done

    # ADM this is just to count bricks files in _update_status.
    nbrick = np.zeros((), dtype='i8')
    t0 = time()

    def _update_status(result):
        """wrapper function for the main parallel process"""
        if nbrick % 250 == 0 and nbrick > 0:
            elapsed = (time()-t0)/60.
            rate = nbrick/elapsed/60.
            log.info('{}/{} bricks; {:.1f} bricks/sec...t = {:.1f} mins'
                     .format(nbrick, nbricks, rate, elapsed))
        nbrick[...] += 1

        return result

    # ADM to parallelize or not to parallelize.
    if numproc > 1 and nbricks > 0:
        pool = sharedmem.MapReduce(np=numproc)
        with pool:
            nexp = pool.map(_get_nexp, bricknames, reduce=_update_status)
    else:
        nexp = []
        for brickname in bricknames:
            nexp.append(_update_status(_get_nexp(brickname)))

    if len(nexp) > 0:
        nexp = np.concatenate(nexp)

    # ADM match back on TARGETID to maintain original order of targets.
    ii = match_to(nexp["TARGETID"], targs["TARGETID"])

    return nexp[ii]


def write_nexp_for_target_file(targfile, drdir, outdir, numproc=1,
                               overwrite=True):
    """Write a file of pixel-level NOBS for one target file.

    Parameters
    ----------
    targfile : :class:`str`
        Full path to a target file.
    drdir : :class:`str`
        Root directory for a Data Release from the Legacy Surveys
        e.g. /global/project/projectdirs/cosmo/data/legacysurvey/dr9.
    outdir : :class:`str`
        The directory to which to write output files. This will be
        created if it doesn't yet exist. The actual file is written
        to <outdir> + pixel-nobs-<targfile>.
    numproc : :class:`int`, optional, defaults to 1
        The number of processes to parallelize across. The default is
        to run the code in serial.
    overwrite : :class:`bool`, optional, defaults to ``True``
        If ``False`` then don't overwrite the output file if it exists.
        The code will just proceed as if nothing happened. This is useful
        for "mopping up" missing files by running all `targfile`s with
        `overwrite`=``False`` and letting the code skip existing files.

    Returns
    -------
    Nothing, but a file of RA/DEC/BRICKNAME/TARGETID/PIXEL_NOBS_G/R/Z is
    written to the `outdir`. The filename is the same as the input
    `targfile` filename, but prepended with pixel-nobs.

    Notes
    -----
    - Useful as the NOBS listed in the target files differs from the
      pixel-level NOBS assigned to the DESI random catalogs.
    """
    # ADM create the name of the output file.
    outfn = os.path.join(
        outdir, "pixel-nobs-{}".format(os.path.basename(targfile)))

    # ADM only return if overwriting is turned off and the file exists.
    if os.path.isfile(outfn) and not overwrite:
        log.info("Refusing to overwrite {}".format(outfn))
        return

    # ADM look up the nexp information.
    nexp = make_nexp_for_target_file(targfile, drdir, numproc=numproc)

    # ADM copy the header of the input file.
    hdr = fitsio.read_header(targfile, "TARGETS")

    # ADM make the output directory if it doesn't exist.
    os.makedirs(os.path.dirname(outfn), exist_ok=True)

    # ADM write the results.
    log.info("Writing to {}".format(outfn))
    write_with_units(outfn, nexp, extname='PIXEL_NOBS', header=hdr)

    return


def write_nexp_in_healpix(targdir, drdir, outdir, nside=None, pixlist=None,
                          numproc=1, overwrite=True):
    """Write files of pixel-level NOBS for each target file in a directory.

    Parameters
    ----------
    targdir : :class:`str`
        Full path to a directory that contains target files.
    drdir : :class:`str`
        Root directory for a Data Release from the Legacy Surveys
        e.g. /global/project/projectdirs/cosmo/data/legacysurvey/dr9.
    outdir : :class:`str`
        The directory to which to write output files. This will be
        created if it doesn't yet exist. Each file in `targdir` is
        written to <outdir> + pixel-nobs-<targfile>.
    nside : :class:`int`, optional, defaults to `None`
        (NESTED) HEALPix `nside` to use with `pixlist`.
    pixlist : :class:`list` or `int`, optional, defaults to `None`
        Only process files for which the ZEROTH source in the file is
        in a list of (NESTED) HEALpixels at the supplied `nside`.
    numproc : :class:`int`, optional, defaults to 1
        The number of processes to parallelize across. The default is
        to run the code in serial.
    overwrite : :class:`bool`, optional, defaults to ``True``
        If ``False`` then don't overwrite any output files encountered
        if they already exist. This is useful for quickly "mopping up"
        missing files by letting the code skip existing work.

    Returns
    -------
    :class:`int`
        The number of files written to the outdir. The files contain
        columns RA/DEC/BRICKNAME/TARGETID/PIXEL_NOBS_G/R/Z. The filenames
        are the same as the input target filenames in `targdir`, but
        prepended with pixel-nobs.

    Notes
    -----
    - Useful as the NOBS listed in the target files differs from the
      pixel-level NOBS assigned to the DESI random catalogs.
    - Pass `pixlist`=``None`` to process ALL files in `targdir`.
    """
    # ADM make an array of all input files.
    targfiles = np.array(sorted(glob(os.path.join(targdir, "targets*fits"))))
    log.info("Processing {} files from {}...t={:.1f}s".format(
        len(targfiles), targdir, time()-start))

    # ADM and limit to passed HEALPixel list, if requested.
    if pixlist is not None:
        inhp = []
        for targfile in targfiles:
            zeroth = fitsio.read(targfile, rows=0)
            inhp.append(is_in_hp(zeroth, nside, pixlist))
        inhp = np.concatenate(inhp)
        targfiles = targfiles[inhp]
        log.info("Limiting to {} files in pixlist={}...t={:.1f}s".format(
            len(targfiles), pixlist, time()-start))

    # ADM find the nexp information and write it for each file.
    for targfile in targfiles:
        write_nexp_for_target_file(targfile, drdir, outdir, numproc=numproc,
                                   overwrite=overwrite)
    log.info("Done...t={:.1f}s".format(time()-start))

    return
