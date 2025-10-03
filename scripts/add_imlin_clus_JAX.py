"""
Script to apply linear regression and get corrective weights for the imaging systematics.

The script is designed to mimic existing LSS scripts (in particular ``add_imlin_clus.py``), and to run on clustering catalogs.
"""

import argparse
import logging
import os
from pathlib import Path
import sys

import fitsio
import LSS.common_tools as common
import numpy as np
from astropy.table import Table, join, vstack
from LSS.globals import main

## some globals
nside = 256
nest = True

## Logging setup

logname = "mkCat"
logger = logging.getLogger(logname)
logger.setLevel(logging.INFO)
# create console handler and set level to debug
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
# create formatter
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
# add formatter to ch
ch.setFormatter(formatter)
# add ch to logger
logger.addHandler(ch)

# Define the argument parser
parser = argparse.ArgumentParser()

parser.add_argument("--type", help="Tracer type to be selected")
parser.add_argument(
    "--basedir",
    help="Base directory for output, default is SCRATCH",
    default=os.environ["SCRATCH"],
)
parser.add_argument(
    "--version",
    help="Catalog version; use 'test' unless you know what you are doing!",
    default="test",
)
parser.add_argument(
    "--survey", help="e.g., main (for all), DA02, any future DA", default="DA2"
)
parser.add_argument("--verspec", help="Version for redshifts", default="loa-v1")
parser.add_argument(
    "--usemaps",
    help="List of maps to use; defaults to what is set by globals",
    type=str,
    nargs="*",
    default=None,
)

# Updated flags
parser.add_argument(
    "--exclude_debv", help="Exclude EBV_DIFF_* maps from fit_maps", action="store_true"
)
parser.add_argument(
    "--relax_zbounds", help="Use less restricted redshift bounds", action="store_true"
)
parser.add_argument(
    "--imsys_finezbin",
    help="Perform imaging regressions in dz=0.1 bins",
    action="store_true",
)
parser.add_argument(
    "--imsys_1zbin",
    help="Perform imaging regressions in just 1 z bin",
    action="store_true",
)
parser.add_argument(
    "--imsys_clus",
    help="Add weights for imaging systematics using eboss method, applied to clustering catalogs",
    action="store_true",
)
parser.add_argument(
    "--imsys_clus_ran",
    help="Add imaging weights to random catalogs",
    action="store_true",
)
parser.add_argument(
    "--replace_syscol",
    help="Replace any existing WEIGHT_SYS with new weights",
    action="store_true",
)
parser.add_argument(
    "--add_syscol2blind",
    help="Add the new weight column to the blinded catalogs",
    action="store_true",
)
parser.add_argument(
    "--Y1_mode",
    help="use this with LRGs to reproduce the Y1/DR2 BAO behavior",
    action="store_true",
)

# Other arguments
parser.add_argument(
    "--imsys_zbin",
    help="if yes, do imaging systematic regressions in z bins",
    default="y",
)
parser.add_argument(
    "--par", help="run different random number in parallel?", default="y"
)
parser.add_argument(
    "--extra_clus_dir",
    help="Optional extra directory structure for clustering catalog",
    default="fNL/",
)
parser.add_argument("--minr", help="Minimum random file index", default=0, type=int)
parser.add_argument(
    "--maxr", help="Maximum random file index (18 available)", default=18, type=int
)
parser.add_argument(
    "--syscol",
    help="Name for new systematic column; determined automatically if None",
    default=None,
)
parser.add_argument(
    "--nran4imsys",
    help="Number of random files to use for linear regression",
    default=10,
    type=int,
)


args = parser.parse_args()
logger.info("Parsed arguments: %s", args)

# Use arguments to define various paths and get the redshift range of the tracer

tracer_type = args.type
basedir = args.basedir
version = args.version
specrel = args.verspec
rm = int(args.minr)
rx = int(args.maxr)

mainp = main(
    args.type,
    args.verspec,
    survey=args.survey,
    relax_zbounds="y" if args.relax_zbounds else "n",
)
zmin = mainp.zmin
zmax = mainp.zmax

maindir = os.path.join(basedir, args.survey, "LSS")

ldirspec = os.path.join(maindir, specrel)

dirout = os.path.join(ldirspec, "LSScats", version)

dirin = dirout
lssmapdirout = os.path.join(dirout, "hpmaps")

# define column name
if args.syscol is None:
    if args.imsys_zbin == "y":
        syscol = "WEIGHT_IMLIN"
    if args.imsys_1zbin:
        syscol = "WEIGHT_IMLIN_1ZBIN"
    if args.imsys_finezbin:
        syscol = "WEIGHT_IMLIN_FINEZBIN"
else:
    syscol = args.syscol

# define maps to use for the regression
if args.usemaps is None:
    fit_maps = mainp.fit_maps
    # if args.imsys_finezbin: # disabled for now
    #     fit_maps = mainp.fit_maps_all
elif args.usemaps[0] == "all":
    fit_maps = mainp.fit_maps_all
    syscol += "_ALL"
elif args.usemaps[0] == "allebv":
    fit_maps = mainp.fit_maps_allebv
    syscol += "_ALLEBV"
elif args.usemaps[0] == "allebvcmb":
    fit_maps = mainp.fit_maps_allebvcmb
    syscol += "_ALLEBVCMB"
else:
    fit_maps = [mapn for mapn in args.usemaps]

# Remove EBV_DIFF_* if exclude_debv is set
if args.exclude_debv:
    fit_maps = [m for m in fit_maps if not m.startswith("EBV_DIFF_")]
    syscol = (syscol or "_WEIGHT_IMLIN") + "_NODEBV"

if args.Y1_mode:
    syscol = args.syscol
    fit_maps = mainp.fit_maps
    logger.info("Fit maps are set to Y1 choices; weight column will be %s", syscol)

logger.info("Using %s as column name for the computed weights.", syscol)
logger.info("Using the following systematics maps for the regression: %s", fit_maps)

# Define the redshift ranges for the regression
zl = (zmin, zmax)  # General redshift range

tpstr = args.type
if "BGS" in tracer_type:
    tpstr = "BGS_BRIGHT"
if "LRG" in tracer_type:
    tpstr = "LRG"
inds = np.arange(rm, rx)

if tracer_type[:3] == "ELG":
    if args.imsys_zbin == "y":
        zrl = [(0.8, 1.1), (1.1, 1.6)]
    elif args.imsys_1zbin:
        zrl = [(0.8, 1.6)]
    zsysmin = 0.8
    zsysmax = 1.6


if tracer_type[:3] == "QSO":
    if args.imsys_zbin == "y":
        zrl = [(0.8, 1.3), (1.3, 2.1), (2.1, 3.5)]
    elif args.imsys_1zbin:
        zrl = [(0.8, 3.5)]
    zsysmin = 0.8
    zsysmax = 3.5
if tracer_type[:3] == "LRG":
    if args.imsys_zbin == "y":
        zrl = [(0.4, 0.6), (0.6, 0.8), (0.8, 1.1)]
    elif args.imsys_1zbin:
        zrl = [(0.4, 1.1)]
    zsysmin = 0.4
    zsysmax = 1.1
    if args.relax_zbounds:
        zsysmax = 1.2
        zsysmin = 0.3

if "BGS_BRIGHT-" in tracer_type:
    zrl = [(0.1, 0.4)]
elif tracer_type[:3] == "BGS":
    zrl = [(0.01, 0.5)]
    zmin = 0.01
    zmax = 0.5

logger.info("The added weight column will be named %s", syscol)

debv = common.get_debv()
zcmb = common.mk_zcmbmap()
sky_g, sky_r, sky_z = common.get_skyres()

# Do the regression
if args.imsys_clus:
    from LSS.imaging.systematics_linear_regression import make_fit_maps_dictionary, produce_imweights

    # define the paths for the input files (loading is deferred to ``produce_imweights``)
    fname_ngc_out = os.path.join(
        dirout, args.extra_clus_dir, f"{tracer_type}_NGC_clustering.dat.fits"
    )

    fname_sgc_out = os.path.join(
        dirout, args.extra_clus_dir, f"{tracer_type}_SGC_clustering.dat.fits"
    )

    # get paths for random catalogs (loading is deferred to ``produce_imweights``)
    randoms_fnames_out = [
        os.path.join(
            dirout,
            args.extra_clus_dir,
            f"{tracer_type}_NGC_{i}_clustering.ran.fits",
        )
        for i in range(args.nran4imsys)
    ] + [
        os.path.join(
            dirout,
            args.extra_clus_dir,
            f"{tracer_type}_SGC_{i}_clustering.ran.fits",
        )
        for i in range(args.nran4imsys)
    ]

    # If on NERSC and applicable, for the INPUT files, switch to dvs_ro
    def global_to_dvs_ro(path: str) -> str:
        """
        If the root directory of the path is ``global``, returns the same path using ``dvs_ro`` instead. Otherwise, returns the original path unchanged.

        Parameters
        ----------
        path : str
            Any path.

        Returns
        -------
        str
            Same as input path with leading ``global`` switched to ``dvs_ro`` if applicable.
        """
        posixpath = Path(path)
        path_parts = posixpath.parts
        if (
            len(path_parts) > 2
            and posixpath.is_absolute()
            and path_parts[1] == "global"
        ):
            posixpath = Path("/dvs_ro").joinpath(*path_parts[2:])
        return str(posixpath)

    fname_sgc_in = global_to_dvs_ro(fname_sgc_out)
    fname_ngc_in = global_to_dvs_ro(fname_ngc_out)
    randoms_fnames_in = [
        global_to_dvs_ro(randoms_fname) for randoms_fname in randoms_fnames_out
    ]

    # Get redshift ranges
    if args.imsys_finezbin:
        dz = 0.1
        zm = zsysmin
        zx = zm + dz
        redshift_ranges = [(zm, zx)]
        while zm < zsysmax:
            zx = zm + dz
            zx = round(zx, 1)
            redshift_ranges += [(zm, zx)]
            zm = zx
        use_maps = fit_maps

    elif args.imsys_1zbin:
        redshift_ranges = [(zmin, zmax)]
        use_maps = fit_maps

    elif args.imsys_zbin == "y":
        redshift_ranges = zrl
        if tracer_type == "LRG" and args.Y1_mode:
            use_maps = make_fit_maps_dictionary(
                default=fit_maps,
                except_when=[
                    ("S", (0.4, 0.6), mainp.fit_maps46s),
                    ("S", (0.6, 0.8), mainp.fit_maps68s),
                    ("S", (0.8, 1.1), mainp.fit_maps81s),
                ],
                regions=["N", "S"],
                redshift_range=zrl,
            )
        else:
            use_maps = fit_maps
    else:
        sys.exit("no valid z binning choice in arguments, exiting...")

    # perform regression
    weights = produce_imweights(
        data_catalog_paths=[fname_sgc_in, fname_ngc_in],
        random_catalogs_paths=randoms_fnames_in,
        is_clustering_catalog=True,
        weight_scheme=None,
        tracer_type=tracer_type,
        redshift_range=redshift_ranges,
        templates_maps_path_S=os.path.join(
            lssmapdirout, f"{tpstr}_mapprops_healpix_nested_nside{nside}_S.fits"
        ),
        templates_maps_path_N=os.path.join(
            lssmapdirout, f"{tpstr}_mapprops_healpix_nested_nside{nside}_N.fits"
        ),
        fit_maps=use_maps,
        output_directory=os.path.join(dirout, args.extra_clus_dir),
        output_catalog_path=None,  # writing to disk will be done later to handle SGC/NGC separately
        output_column_name=syscol,
        save_summary_plots=True,
        nbins=10,  # is the default
        tail=0.5,  # is the default
        logger=logger,
        loglevel="INFO",
        templates_maps_nside=nside,
        templates_maps_nested=nest,
    )

    ## attach data to NGC/SGC catalogs, write those out

    # Need to load the data individual data catalogs again
    data_sgc = Table.read(fname_sgc_in)
    data_ngc = Table.read(fname_ngc_in)
    # Catalogs are just concatenated in the order SGC, NGC
    # so this is enough to assign weights to the correct one
    transition_index = len(data_sgc)
    assert transition_index + len(data_ngc) == len(weights), "Shape mismatch!"
    # add custom column to catalog
    data_sgc[syscol] = weights[:transition_index]
    data_ngc[syscol] = weights[transition_index:]
    # overwrite the WEIGHT columns
    if args.replace_syscol:
        data_sgc["WEIGHT"] /= data_sgc["WEIGHT_SYS"]
        data_sgc["WEIGHT_SYS"] = data_sgc[syscol]
        data_sgc["WEIGHT"] *= data_sgc["WEIGHT_SYS"]

        data_ngc["WEIGHT"] /= data_ngc["WEIGHT_SYS"]
        data_ngc["WEIGHT_SYS"] = data_ngc[syscol]
        data_ngc["WEIGHT"] *= data_ngc["WEIGHT_SYS"]
    # write out everything
    common.write_LSS_scratchcp(
        data_sgc,
        fname_sgc_out,
        logger=logger,
    )
    common.write_LSS_scratchcp(
        data_ngc,
        fname_ngc_out,
        logger=logger,
    )

# Optionally also write the weights in the randoms
if args.imsys_clus_ran:
    fname = os.path.join(
        dirout, args.extra_clus_dir, f"{tracer_type}_NGC_clustering.dat.fits"
    )
    dat_ngc = Table(fitsio.read(fname, columns=["TARGETID", syscol]))
    fname = os.path.join(
        dirout, args.extra_clus_dir, f"{tracer_type}_SGC_clustering.dat.fits"
    )
    dat_sgc = Table(fitsio.read(fname, columns=["TARGETID", syscol]))
    dat = vstack([dat_sgc, dat_ngc])
    dat.rename_column("TARGETID", "TARGETID_DATA")
    regl = ["NGC", "SGC"]
    syscolr = syscol

    # if args.replace_syscol == 'y':
    #    syscolr = 'WEIGHT_SYS'
    def _add2ran(rann):
        for reg in regl:
            ran_fn = os.path.join(
                dirout,
                args.extra_clus_dir,
                f"{tracer_type}_{reg}_{rann}_clustering.ran.fits",
            )
            ran = Table(fitsio.read(ran_fn))
            if syscolr in ran.colnames:
                ran.remove_column(syscolr)
            ran = join(ran, dat, keys=["TARGETID_DATA"])
            if args.replace_syscol:
                ran["WEIGHT"] /= ran["WEIGHT_SYS"]
                ran["WEIGHT_SYS"] = ran[syscolr]
                ran["WEIGHT"] *= ran["WEIGHT_SYS"]
            common.write_LSS_scratchcp(ran, ran_fn, logger=logger)

    if args.par == "y":
        from multiprocessing import Pool

        with Pool() as pool:
            res = pool.map(_add2ran, inds)
    else:
        for rn in inds:  # range(rm,rx):
            _add2ran(rn)

# Optionally also write the weights in the blinded catalogs
if args.add_syscol2blind:
    syscolr = syscol
    # if args.replace_syscol == 'y':
    #    syscolr = 'WEIGHT_SYS'

    dats = []
    for reg in ["NGC", "SGC"]:
        fname = os.path.join(
            dirout,
            args.extra_clus_dir,
            f"{tracer_type}_{reg}_clustering.dat.fits",
        )
        dati = Table(fitsio.read(fname, columns=["TARGETID", syscol]))
        dats.append(dati)
        fname_blind = os.path.join(
            dirout,
            args.extra_clus_dir,
            "blinded",
            f"{tracer_type}_{reg}_clustering.dat.fits",
        )
        dat_blind = Table(fitsio.read(fname_blind))
        if syscol in list(dat_blind.colnames):
            dat_blind.remove_column(syscol)

        dat_blind = join(dat_blind, dati, keys=["TARGETID"])
        if args.replace_syscol:
            dat_blind["WEIGHT"] /= dat_blind["WEIGHT_SYS"]
            dat_blind["WEIGHT_SYS"] = dat_blind[syscol]
            dat_blind["WEIGHT"] *= dat_blind["WEIGHT_SYS"]

        common.write_LSS_scratchcp(dat_blind, fname_blind, logger=logger)
    dat = vstack(dats)
    dat.rename_column("TARGETID", "TARGETID_DATA")
    regl = ["NGC", "SGC"]

    def _add2ranblind(rann):
        for reg in regl:
            ran_fn = os.path.join(
                dirout,
                args.extra_clus_dir,
                "blinded",
                f"{tracer_type}_{reg}_{rann}_clustering.ran.fits",
            )
            ran = Table(fitsio.read(ran_fn))
            if syscolr in ran.colnames:
                ran.remove_column(syscolr)
            ran = join(ran, dat, keys=["TARGETID_DATA"])
            if args.replace_syscol:
                ran["WEIGHT"] /= ran["WEIGHT_SYS"]
                ran["WEIGHT_SYS"] = ran[syscolr]
                ran["WEIGHT"] *= ran["WEIGHT_SYS"]

            common.write_LSS_scratchcp(ran, ran_fn, logger=logger)

    if args.par == "y":
        from multiprocessing import Pool

        with Pool() as pool:
            res = pool.map(_add2ranblind, inds)
    else:
        for rn in inds:  # range(rm,rx):
            _add2ranblind(rn)
