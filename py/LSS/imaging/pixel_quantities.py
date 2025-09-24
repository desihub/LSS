# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
LSS.pixel_quantities
====================

Routines for adding pixel-level quantities to target files.

.. _`bitmasks page`: https://www.legacysurvey.org/dr9/bitmasks
"""
import os
import fitsio
import numpy as np
from time import time
import astropy.io.fits as fits
from astropy.wcs import WCS
from glob import glob
import healpy as hp

from desitarget.randoms import dr_extension, quantities_at_positions_in_a_brick
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


def make_slurm_script(targdir, drdir, outdir, nside=2, numproc=60, mopup=False):
    """Write an example slurm script for parallelization

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
        written to <outdir> + pixel-<targfile>.
    nside : :class:`int`, optional, defaults to 2
        (NESTED) HEALPix `nside` to use with `pixlist`.
    numproc : :class:`int`, optional, defaults to 60
        The number of processes to parallelize across.
    mopup : :class:`bool`, optional, defaults to ``False``
        If ``True`` then do NOT overwrite existing output files. This is
        useful for "mopping up" failed or missing files.

    Returns
    -------
    Nothing, but a bash script that can be used to parallelize
    pixel-level lookups for targets is written to screen.
    """
    npix = hp.nside2npix(nside)

    print('#!/bin/bash -l')
    print('#SBATCH -q regular')
    print('#SBATCH -N 24')
    print('#SBATCH -t 04:00:00')
    print('#SBATCH -L SCRATCH,project')
    print('#SBATCH -C haswell')
    print('')
    for pixnum in range(npix):
        former = 'srun -N 1 get_pixel_quantities_for_targets '
        former += '--nside {} --healpixels {} --numproc {} {} {} {}'
        msg = former.format(nside, pixnum, numproc, targdir, drdir, outdir)
        if mopup:
            msg += ' --mopup'
        msg += " &"
        print(msg)

    return


def at_locations_in_a_brick(ras, decs, brickname, nors, drdir):
    """Pixel-level quantities at locations in a Legacy Surveys brick.

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
    :class:`~numpy.array`
        With the following quantities at each location in the brick:

        PIXEL_NOBS_G, R, Z:
            Number of observations in g, r, z-band.
        PIXEL_PSFDEPTH_G, R, Z:
            PSF depth at this location in g, r, z.
        PIXEL_GALDEPTH_G, R, Z:
            Galaxy depth in g, r, z.
        PIXEL_PSFDEPTH_W1, W2:
            (PSF) depth in W1, W2 (AB mag system).
        PIXEL_PSFSIZE_G, R, Z:
            Weighted average PSF FWHM (arcsec).
        PIXEL_APFLUX_G, R, Z:
            Sky background in a 0.75-arcsec-radius DESI aperture.
        PIXEL_APFLUX_IVAR_G, R, Z:
            Inverse variance of sky background.
        PIXEL_MASKBITS:
            Mask information. See the Legacy Surveys `bitmasks page`_.
        PIXEL_WISEMASK_W1:
            Mask information. See the Legacy Surveys `bitmasks page`_.
        PIXEL_WISEMASK_W2:
            Mask information. See the Legacy Surveys `bitmasks page`_.
x
    Notes
    -----
    - Wraps :func:`desitarget.quantities_at_positions_in_a_brick()`
    """
    # ADM check nors is one of the allowed strings.
    if nors not in ["north", "south"]:
        msg = 'nors not "north"/"south" for brick '.format(brickname)
        log.critical(msg)
        raise ValueError(msg)

    # ADM expand the drdir to include "north"/"south".
    drdir = os.path.join(drdir, nors)

    # ADM output from desitarget.quantities_at_positions_in_a_brick()
    q = quantities_at_positions_in_a_brick(ras, decs, brickname, drdir)

    # ADM restructure the output dictionary.
    dt = [('NOBS_G', 'i2'), ('NOBS_R', 'i2'), ('NOBS_Z', 'i2'),
          ('PSFDEPTH_G', 'f4'), ('PSFDEPTH_R', 'f4'), ('PSFDEPTH_Z', 'f4'),
          ('GALDEPTH_G', 'f4'), ('GALDEPTH_R', 'f4'), ('GALDEPTH_Z', 'f4'),
          ('PSFDEPTH_W1', 'f4'), ('PSFDEPTH_W2', 'f4'),
          ('PSFSIZE_G', 'f4'), ('PSFSIZE_R', 'f4'), ('PSFSIZE_Z', 'f4'),
          ('APFLUX_G', 'f4'), ('APFLUX_R', 'f4'), ('APFLUX_Z', 'f4'),
          ('APFLUX_IVAR_G', 'f4'), ('APFLUX_IVAR_R', 'f4'), ('APFLUX_IVAR_Z', 'f4'),
          ('MASKBITS', 'i2'), ('WISEMASK_W1', '|u1'), ('WISEMASK_W2', '|u1')]
    dt = [("PIXEL_{}".format(col), form) for col, form in dt]

    qout = np.zeros(len(ras), dtype=dt)
    for col in qout.dtype.names:
        key = col.split("PIXEL_")[-1].lower()
        # ADM there's a slight syntactical discrepancy between the output
        # ADM dictionary and the desired structure.
        key = key.replace("psfdepth_w", "psfdepth_W")
        qout[col] = q[key]

    return qout


def look_up_for_target_file(targfile, drdir, numproc=1):
    """Look up pixel-level quantities (from coadds) for one target file.

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
        the standard quantities RA, DEC, BRICKNAME, TARGETID, DESI_TARGET
        and PHOTSYS and the added quantities described in the docstring
        of :func:`at_locations_in_a_brick()`.
    """
    # ADM read in the needed target columns...
    targs = fitsio.read(targfile, columns=["RA", "DEC", "BRICKNAME",
                                           "DESI_TARGET", "TARGETID", "PHOTSYS"])
    # ADM ...and compile the list of brick names in the file.
    bricknames = list(set(targs["BRICKNAME"]))
    nbricks = len(bricknames)

    # ADM wrapper to facilitate parallelization.
    def _get_quantities(brickname):
        """wrapper on at_locations_in_a_brick() for a brick name"""
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
        q = at_locations_in_a_brick(brick["RA"], brick["DEC"], brickname,
                                    nors, drdir)

        # ADM make a table of all of the required information...
        dt = brick.dtype.descr + q.dtype.descr
        done = np.zeros(len(brick), dtype=dt)
        for col in brick.dtype.names:
            done[col] = brick[col]
        for col in q.dtype.names:
            done[col] = q[col]

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
            q = pool.map(_get_quantities, bricknames, reduce=_update_status)
    else:
        q = []
        for brickname in bricknames:
            q.append(_update_status(_get_quantities(brickname)))

    if len(q) > 0:
        q = np.concatenate(q)

    # ADM match back on TARGETID to maintain original order of targets.
    ii = match_to(q["TARGETID"], targs["TARGETID"])

    return q[ii]


def write_for_target_file(targfile, drdir, outdir, numproc=1,
                          overwrite=True):
    """Write a FITS file of pixel-level quantities for one target file.

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
    Nothing, but a file with the quantities documented in the docstring
    of :func:`look_up_for_target_file()` is written to the `outdir`. The
    filename is the same as the input `targfile` filename, but prepended
    with "pixel-".
    """
    # ADM create the name of the output file.
    outfn = os.path.join(
        outdir, "pixel-{}".format(os.path.basename(targfile)))

    # ADM only return if overwriting is turned off and the file exists.
    if os.path.isfile(outfn) and not overwrite:
        log.info("Refusing to overwrite {}".format(outfn))
        return

    # ADM look up the pixel-level information.
    q = look_up_for_target_file(targfile, drdir, numproc=numproc)

    # ADM copy the header of the input file.
    hdr = fitsio.read_header(targfile, "TARGETS")

    # ADM make the output directory if it doesn't exist.
    os.makedirs(os.path.dirname(outfn), exist_ok=True)

    # ADM write the results.
    log.info("Writing to {}".format(outfn))
    write_with_units(outfn, q, extname='PIXEL_TARGETS', header=hdr)

    return


def write_in_healpix(targdir, drdir, outdir, nside=None, pixlist=None,
                     numproc=1, overwrite=True):
    """Write pixel-level quantities for each target file in a directory.

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
        written to <outdir> + pixel-<targfile>.
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
        The number of files written to the `outdir`. Files contain the
        quantities returned by :func:`look_up_for_target_file()`. Each
        file is written to the `outdir`. The filename is the same as the
        input `targfile` filename, but prepended with "pixel-".

    Notes
    -----
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

    # ADM find the pixel-level quantities and write them for each file.
    nfiles = len(targfiles)
    for nfile, targfile in enumerate(targfiles):
        write_for_target_file(targfile, drdir, outdir, numproc=numproc,
                              overwrite=overwrite)
        if nfile % 10 == 0 and nfile > 0:
            elapsed = (time()-start)/60.
            rate = 60.*elapsed/nfile
            log.info('PROCESSED {}/{} files; {:.1f} secs/file...t = {:.1f} mins'
                     .format(nfile, nfiles, rate, elapsed))

    log.info("Done...t={:.1f}s".format(time()-start))

    return nfiles
