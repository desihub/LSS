#!/usr/bin/env python

import os
import argparse
from multiprocessing import Pool
import numpy as np

import fitsio
from astropy.table import Table
import healpy

import desispec.io
import desiutil.healpix
from desiutil.log import get_logger


def _filter_spectra_single_pixout(pix_out, outdir, skip_resolution,
                        catalog, specprod_dir, pixel_scheme_in, 
                        include_single_exp, nside_in, available_upix):
    """
    Produce a (single) skimmed spectra file.

    Parameters
    ----------
    pix_out : int
        pixel number for the output file
    outdir : str
        output base directory
    skip_resolution : bool
        if true, do not incluse resolution matrix in new coadd files
    catalog : Table
        input catalog, must include 'TARGETID' and 'OUTPUT_PIX' columns
    specprod_dir : str
        input base directory for spectra
    pixel_scheme_in : str
        HEALPIX64 or UNIQPIX
    include_single_exp : bool
        if true, spectra files containing individual exposures are also written
    nside_in, available_upix
        information associated to the input pixellization
    """       
    m = (catalog['OUTPUT_PIX'] == pix_out)
    pix_cat = catalog[m]

    #- Get list of input pixels which contribute to pix_out
    if pixel_scheme_in == 'HEALPIX64':
        if 'HPXPIXEL' in pix_cat.keys():  # faster
            pixels_in = pix_cat['HPXPIXEL']
        else:  # slower, but should give the same
            pixels_in = healpy.ang2pix(nside_in,
                                pix_cat['TARGET_RA'], pix_cat['TARGET_DEC'],
                                lonlat=True, nest=True)
    elif pixel_scheme_in == 'UNIQPIX':
        if 'UNIQPIX' in pix_cat.keys():  # faster
            pixels_in = pix_cat['UNIQPIX']
        else:  # slower, but should give the same
            pixels_in = desiutil.healpix.find_upix(pix_cat['TARGET_RA'], 
                                pix_cat['TARGET_DEC'], available_upix)
    pixels_in = np.unique(pixels_in)

    #- Read coadd spectra files and produce a single skimmed coadd file
    coadd_list = []
    if include_single_exp:
        spectra_list = []
    # list of which hdus to skip
    skip_hdus = None
    if skip_resolution:
        skip_hdus = ['RESOLUTION']
    for pix_in in pixels_in:
        coadd_file = os.path.join(specprod_dir, str(pix_in//100), str(pix_in),
                        f'coadd-main-dark-{pix_in}.fits')
        coadd = desispec.io.read_spectra(coadd_file,
                                     targetids=pix_cat['TARGETID'],
                                     skip_hdus=skip_hdus)
        coadd_list.append(coadd)

        #- Optionally read spectra files and produce single skimmed spectra file
        if include_single_exp:
            spectra_file = os.path.join(specprod_dir, str(pix_in//100), str(pix_in),
                            f'spectra-main-dark-{pix_in}.fits')
            #- Assume same skipped hdus as coadd files
            spectra = desispec.io.read_spectra(spectra_file,
                                         targetids=pix_cat['TARGETID'],
                                         skip_hdus=skip_hdus)
            spectra_list.append(spectra)


    if len(coadd_list)>0:
        coadd_out = desispec.spectra.stack(coadd_list)
        subdir_out = os.path.join(outdir, str(pix_out//100), str(pix_out))
        os.makedirs(subdir_out, exist_ok=True)
        desispec.io.write_spectra(f'{subdir_out}/coadd-main-dark-{pix_out}.fits',
                                    coadd_out)

        if include_single_exp:
            spectra_out = desispec.spectra.stack(spectra_list)
            desispec.io.write_spectra(f'{subdir_out}/spectra-main-dark-{pix_out}.fits',
                                    spectra_out)

    return 0


def _pixels_from_catalog(catalog, nside_out, nside_in=0):
    """
    Return the array of HEALPIX pixels (NESTED ordering with nside_out) 
    matching entries of an input catalog.

    Parameters
    ----------
    catalog : Table
        input catalog, must have a column 'HPXPIXEL' or columns 'TARGET_RA/DEC'.
    nside_out : int
        output NSIDE
    nside_in : int
        input NSIDE (optional, to be used with the catalog's 'HPXPIXEL' column)
    
    Returns
    -------
    array
        array of HEALPIX pixels, same size as input catalog
    """

    # Fast method:
    if 'HPXPIXEL' in catalog.keys() and nside_in>=nside_out:
        if nside_in == nside_out:
            return catalog['HPXPIXEL']
        else:
            ratio = nside_in / nside_out
            if not ratio.is_integer():
                raise ValueError("nside_in/nside_out isn't an integer.")
            factor = int(ratio**2)
            return catalog['HPXPIXEL'] // factor
    # Slow method:
    else:
        out_pixels = healpy.ang2pix(nside_out, 
                              catalog['TARGET_RA'],
                              catalog['TARGET_DEC'],
                              lonlat=True, nest=True)
        return out_pixels


def _parse(options=None):

    parser = argparse.ArgumentParser(
        description='Filter DESI spectra from a catalog. Primarily intended for quasars (Lya WG)')

    parser.add_argument('--catalog', type=str,
        help='input catalog file (must have TARGETID columns)')
    parser.add_argument('--specprod-dir', type=str,
        help='input data reduction directory (eg. $DESI_SPECTRO_REDUX/loa/healpix/main/dark)')
    parser.add_argument('--pixel-scheme-in', type=str,
        help='pixellization scheme used in the data reduction: must be HEALPIX64 or UNIQPIX')
    parser.add_argument('--uniqpix-filename', type=str, default='uniqpix-main-dark.fits',
        help='FITS file with list of uniqpixels in the input data reduction directory')
    parser.add_argument('-o', '--outdir', type=str,
        help='output spectra directory tree')
    parser.add_argument('--nside-out', type=int, default=16,
        help='NSIDE for the output spectra directory tree (default: 16)')
    parser.add_argument('--ncpu', type=int, default=1,
        help='number of CPUs for pool (default: 1)')
    parser.add_argument('--skip-resolution', action="store_true",
        help='do not include resolution matrices in skimmed spectra')
    parser.add_argument('--include-single-exposures', action='store_true',
        help='produce skimmed files from both coadds and (single-exposure) spectra files')

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)

    return args


if __name__ == "__main__":

    args = _parse()
    log = get_logger()

    if args.pixel_scheme_in == 'HEALPIX64':
        nside_in = 64
        available_upix = None
    elif args.pixel_scheme_in == 'UNIQPIX':
        nside_in = 0
        uniqpixfile = os.path.join(args.specprod_dir, args.uniqpix_filename)
        available_upix = Table.read(uniqpixfile)['UNIQPIX'].data
    else:
        raise ValueError('Wrong input value for --pixel-scheme-in')

    skimming_catalog = Table(fitsio.read(args.catalog, ext=1))
    log.info("Catalog loaded.")

    if args.include_single_exposures:
        log.info("Both skimmed spectra and coadd files will be produced")
    else:
        log.info("Only skimmed coadd files will be produced")

    skimming_catalog['OUTPUT_PIX'] = _pixels_from_catalog(skimming_catalog,
                                        args.nside_out, nside_in=nside_in)
    output_pixels = np.unique(skimming_catalog['OUTPUT_PIX'])
    #output_pixels = output_pixels[12:18]  # debug

    if args.ncpu==1:
        for pix_out in output_pixels:
            res = _filter_spectra_single_pixout(pix_out, args.outdir,
                        args.skip_resolution, skimming_catalog, args.specprod_dir,
                        args.pixel_scheme_in, args.include_single_exposures, 
                        nside_in, available_upix)
    else:
        list_args = [
            [pix_out, args.outdir,
            args.skip_resolution, skimming_catalog, args.specprod_dir,
            args.pixel_scheme_in, args.include_single_exposures, nside_in, 
            available_upix]
            for pix_out in output_pixels
        ]
        with Pool(args.ncpu) as pool:
            res = pool.starmap(_filter_spectra_single_pixout, list_args)

    log.info('Done.')

