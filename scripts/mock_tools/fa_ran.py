#!/usr/bin/env python

import os
import sys
from glob import glob
import numpy as np
import astropy.io.fits as fits
from astropy.table import Table, vstack
from astropy import units
from astropy.coordinates import SkyCoord
import fitsio
from datetime import datetime, timezone
from time import time

# AR desitarget
import desitarget
from desitarget.targets import encode_targetid
from desitarget.targets import zcut, midzcut # zcut=2.1 (lya), midzcut=1.6 (loz)
from desitarget.targetmask import obsconditions
from desitarget.io import read_targets_in_tiles
from desitarget.geomask import pixarea2nside, match  # AR to match arrays with no repeats
from desitarget.mtl import make_mtl
from desitarget.internal import sharedmem

# AR fiberassign
import fiberassign
from fiberassign.scripts.assign import parse_assign, run_assign_full
from fiberassign.utils import Logger
from fiberassign.fba_launch_io import get_desitarget_paths, assert_isoformat_utc

# AR desimodel
import desimodel
from desimodel.footprint import is_point_in_desi, tiles2pix
from desiutil.redirect import stdouterr_redirected
import healpy as hp
from argparse import ArgumentParser


# AR created file names
def get_fn(outdir, flavor, steps=None):
    # AR flavor = settings, log, tiles, sky, targ, fastats-targ, fastats-sky
    # AR pass = passid
    if flavor == "tiles":
        fn = os.path.join(
            outdir,
            "faruns",
            "{}.fits".format(flavor),
            )
    elif flavor == "settings":
        fn = os.path.join(outdir, "inputs", "{}-settings.asc".format(steps.replace(",", "-")))
    elif flavor == "log":
        fn = os.path.join(outdir, "outputs", "{}.log".format(steps.replace(",", "-")))
    elif flavor in ["fastats-targ", "fastats-sky", "fastats-pairs"]:
        fn = os.path.join(outdir, "outputs", "{}.asc".format(flavor))
    else:
        fn = os.path.join(outdir, "inputs", "{}.fits".format(flavor))
    return fn


# AR creates one tiles file for with all passes
# AR and one tiles file per pass
# AR formatting assumed to be that of tiles-main.ecsv
def create_tiles(obsfn, program, outdir,infn='/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-main.ecsv'):
    d_orig = Table.read(infn)
    # AR cut on program
    sel = np.array([program.lower() for program in d_orig["PROGRAM"]]) == program
    d_obs = Table.read(obsfn)
    sel &= np.isin(d_orig['TILEID'],d_obs['TILEID'])
    d_orig = d_orig[sel]
    # AR
    d = Table()
    for key in d_orig.dtype.names:
        d[key] = d_orig[key]
    # AR update PROGRAM with args.program, and add OBSCONDITIONS
    d["PROGRAM"] = program
    d["OBSCONDITIONS"] = obsconditions.mask(program.upper())
    # AR cut on the number of passes (first pass is 0)
    fitsio.write(get_fn(outdir, "tiles"), d.as_array(), clobber=True)
    return True


# AR wrapper to read targets
# AR ! ABACUS !
def wrapper_read_targets():
    tiles = fits.open(get_fn(args.outdir, "tiles"))[1].data
    log.info("{:.1f}s\twrapper_read_targets\thpdir={}".format(time() - start, args.infn))
    d = Table.read(args.infn)
    log.info("{:.1f}s\twrapper_read_targets\tread {} targets".format(time() - start, len(d)))
    return d


# AR sky file


# AR internal function to run fiber assignment
# AR on a tile, when args.numproc > 1
# AR uses global variables:
# AR - targ_ras, targ_decs, targ_pixs : input target coordinates
def _do_run_assign_full(intargfn_fadir_footfn_skyfn_targfn):
    global targ_ras, targ_decs, targ_pixs
    # AR decode folder/filenames
    intargfn, fadir, footfn, targfn = intargfn_fadir_footfn_skyfn_targfn.split(
        ","
    )
    # AR targ
    tiles = fits.open(footfn)[1].data
    nside, nest = pixarea2nside(7.), True # AR be careful those are the same as below
    pixlist = tiles2pix(nside, tiles=tiles)
    ii = np.where(np.in1d(targ_pixs, pixlist))[0]
    jj = is_point_in_desi(tiles, targ_ras[ii], targ_decs[ii])
    rows = ii[jj]
    d = fitsio.read(intargfn, rows=rows)
    fitsio.write(targfn, d, clobber=True)
    # AR run fiber assignment
    opts = [
        "--rundate",
        args.rundate,
        "--overwrite",
        "--write_all_targets",
        "--footprint",
        footfn,
        "--dir",
        fadir,
        "--targets",
        targfn,
        "--sky_per_petal",
        args.sky_per_petal,
        "--standards_per_petal",
        args.standards_per_petal,
        "--sky_per_slitblock",
        str(args.sky_per_slitblock),
        "--margin-pos",
        str(args.margin_pos),
        "--margin-gfa",
        str(args.margin_gfa),
        "--margin-petal",
        str(args.margin_petal)
    ]
    log.info(
        "{:.1f}s\t_do_run_assign_full\trun_assign_full for {}".format(
            time() - start, targfn
        )
    )
    ag = parse_assign(opts)
    run_assign_full(ag)
    # AR clean input files
    for fn in [footfn, skyfn, targfn]:
        os.remove(fn)
    return True







def main():

    # AR print start time
    log.info(
        "{:.1f}s\tsettings\tstarting at utc_time={}".format(
            time() - start, utc_time_now_str
        )
    )

    # AR printing settings
    tmpstr = " , ".join(
        [kwargs[0] + "=" + str(kwargs[1]) for kwargs in args._get_kwargs()]
    )
    log.info("{:.1f}s\tsettings\targs: {}".format(time() - start, tmpstr))

    # AR machine
    log.info(
        "{:.1f}s\tsettings\tHOSTNAME={}".format(time() - start, os.getenv("HOSTNAME"))
    )
    # AR fiberassign, desitarget, desimodel code version, path
    for module, name in zip(
        [fiberassign, desitarget, desimodel], ["fiberassign", "desitarget", "desimodel"]
    ):
        log.info(
            "{:.1f}s\tsettings\trunning with {} code version: {}".format(
                time() - start, name, module.__version__
            )
        )
        log.info(
            "{:.1f}s\tsettings\trunning with {} code path: {}".format(
                time() - start, name, module.__path__
            )
        )

    # AR create the tiles files
    if "tiles" in args.steps.split(","):
        # AR then create tiles file
        _ = create_tiles(args.tilesfn, args.program, args.outdir)

    # AR safe "box" area
    tiles = fits.open(get_fn(args.outdir, "tiles"))[1].data


    # AR create the targets file
    # AR ! ABACUS !
    if "targ" in args.steps.split(","):
        log.info("{:.1f}s\ttarg\tstart".format(time() - start))
        tiles = fits.open(get_fn(args.outdir, "tiles"))[1].data
        d = wrapper_read_targets()
        log.info("{:.1f}s\ttarg\tdone reading targets".format(time() - start))
        # AR mtl
        myd = Table(d)
        # AR to store fa stats
        # AR no make_mtl()
        myd["NUMOBS"] = np.zeros(len(myd), dtype=int)
        myd.meta["EXTNAME"] = "MTL"
        # AR no lya, midz_not_loz
        log.info(
            "{:.1f}s\ttarg\tstart writing {}".format(
                time() - start, get_fn(args.outdir, "targ")
            )
        )
        myd.write(get_fn(args.outdir, "targ"), overwrite=True)
        # AR propagating some settings into the PRIMARY header
        fd = fitsio.FITS(get_fn(args.outdir, "targ"), "rw")
        fd["MTL"].write_key("DTVER", args.dtver)
        fd["MTL"].write_key("SURVEY", args.survey)
        fd["MTL"].write_key("PROGRAM", args.program)
        fd["MTL"].write_key("TILESFN", args.tilesfn)
        fd["MTL"].write_key("RUNDATE", args.rundate)
        fd["MTL"].write_key("NSKYPET", args.sky_per_petal)
        fd["MTL"].write_key("NSTDPET", args.standards_per_petal)
        fd["MTL"].write_key("NSKYPSL", args.sky_per_slitblock)
        fd["MTL"].write_key("MARGPOS", args.margin_pos)
        fd["MTL"].write_key("MARGGFA", args.margin_gfa)
        fd["MTL"].write_key("MARGPET", args.margin_petal)
        fd.close()
        log.info(
            "{:.1f}s\ttarg\tdone writing {}".format(
                time() - start, get_fn(args.outdir, "targ")
            )
        )

    # AR run the fiber assignment + update the target file
    if "fa" in args.steps.split(","):
        # AR listing passids
        tiles = fits.open(get_fn(args.outdir, "tiles"))[1].data
        input_targ = get_fn(args.outdir, "targ")
        # AR targ_ras, targ_decs for as global variable if args.numproc > 1
        # AR same for all passids
        if args.numproc > 1:
            global targ_ras, targ_decs, targ_pixs
            d = fitsio.read(input_targ, columns=["RA", "DEC"])
            targ_ras, targ_decs = d["RA"], d["DEC"]
            nside, nest = pixarea2nside(7.), True
            targ_pixs = hp.ang2pix(nside, np.radians((90. - targ_decs)), np.radians(targ_ras), nest=nest)
        fadir = os.path.join(args.outdir, "faruns")
            
        # AR clean fadir folder
        fns = glob(os.path.join(fadir, "fba-*.fits*"))
        if len(fns) > 0:
            for fn in fns:
                os.remove(fn)
        # AR run fiber assignment
        if args.numproc == 1:
            opts = [
                "--rundate",
                args.rundate,
                "--overwrite",
                "--write_all_targets",
                "--footprint",
                get_fn(args.outdir, "tiles", passid=passid),
                "--dir",
                fadir,
                "--sky",
                get_fn(args.outdir, "sky"),
                "--targets",
                input_targ,
                "--sky_per_petal",
                args.sky_per_petal,
                "--standards_per_petal",
                args.standards_per_petal,
            ]
            log.info(
                "{:.1f}s\tfa\trun_assign_full with 1 processor".format(
                    time() - start
                )
            )
            ag = parse_assign(opts)
            run_assign_full(ag)
        else:
            # AR tileids for that passid
            tileids = tiles["TILEID"]
            # AR tiles file for that passid
            foot_pass = get_fn(args.outdir, "tiles")
            # AR preparing files per tileid
            intargfn_fadir_footfn_skyfn_targfn = []
            for tileid in tileids:
                # AR tileid filenames
                footfn = foot_pass.replace(
                    ".fits", "-{:06d}.fits".format(int(tileid))
                )
                targfn = footfn.replace("tiles-", "targ-")
                d = fitsio.read(foot_pass)
                d = d[d["TILEID"] == tileid]
                fitsio.write(footfn, d, clobber=True)
                # AR
                intargfn_fadir_footfn_skyfn_targfn += [
                    ",".join([input_targ, fadir, footfn, targfn])
                ]
            pool = sharedmem.MapReduce(np=args.numproc)
            with pool:
                _ = pool.map(
                    _do_run_assign_full, intargfn_fadir_footfn_skyfn_targfn
                )

    # AR print end time
    log.info(
        "{:.1f}s\tsettings\tfinished at utc_time={}".format(
            time() - start, datetime.now(tz=timezone.utc).isoformat(timespec="seconds")
        )
    )


if __name__ == "__main__":

    # AR allowed steps in desi_fa_smallrun
    steps_all = ["tiles", "targ", "fa"]

    # AR reading arguments
    parser = ArgumentParser()
    parser.add_argument(
        "--infn",
        help="full path to the input target file",
        type=str,
        default=None,
    )
    parser.add_argument(
        "--outdir",
        help="root directory; all other paths are relative to that directory",
        type=str,
        default=None,
    )
    parser.add_argument(
        "--program",
        help="dark or bright or backup (default=dark)",
        type=str,
        default="dark",
        choices=["dark", "bright"],
    )
    parser.add_argument(
        "--survey",
        help="main",
        type=str,
        default="main",
        choices=["main"],
    )
    parser.add_argument(
        "--rundate",
        help="yyyy-mm-ddThh:mm:ss+00:00 rundate for focalplane with UTC timezone formatting (default=current UTC time)",
        type=str,
        default=None,
        required=False,
    )
    parser.add_argument(
        "--dtver",
        help="desitarget version (default=1.1.1)",
        type=str,
        default="1.1.1",
    )
    parser.add_argument(
        "--tilesfn",
        help="file full path containing TILEID of data to use for footprint",
        type=str,
        default=None,
    )
    parser.add_argument(
        "--standards_per_petal",
        help="required number of standards per petal (default=10)",
        type=str,
        default="0",
    )
    parser.add_argument(
        "--sky_per_petal",
        help="required number of sky targets per petal (default=40)",
        type=str,
        default="0",
    )
    parser.add_argument(
        "--sky_per_slitblock",
        help="Required number of sky targets per fiber slitblock (default=1)",
        type=int,
        default=0,
        required=False,
    )
    parser.add_argument("--margin-pos", "--margin_pos", type=float, required=False, default=0.05,
                        help="Add margin (in mm) around positioner keep-out polygons (default: 0.05)")
    parser.add_argument("--margin-petal", "--margin_petal", type=float, required=False, default=0.4,
                        help="Add margin (in mm) around petal-boundary keep-out polygons (default: 0.4)")
    parser.add_argument("--margin-gfa", "--margin_gfa", type=float, required=False, default=0.4,
                        help="Add margin (in mm) around GFA keep-out polygons (default: 0.4)")
    parser.add_argument(
        "--nolog", help="do not print log file (default=n)", type=str, default="n",
    )
    parser.add_argument(
        "--numproc",
        help="number of concurrent processes to use (default=1)",
        type=int,
        default=1,
    )
    parser.add_argument(
        "--steps",
        help="comma-separated list of steps to accomplish, amongst {} (defaults={})".format(steps_all, ",".join(steps_all)),
        type=str,
        default=",".join(steps_all),
    )
    #
    args = parser.parse_args()
    log = Logger.get()
    start = time()


    # AR steps
    wrong_steps = [step for step in args.steps.split(",") if step not in steps_all]
    if len(wrong_steps) > 0:
        log.error("args.steps have the following not authorized steps : {}; exiting".format(
            wrong_steps
            )
        )
        sys.exit(1)

    # AR utc_time_now, rundate
    utc_time_now = datetime.now(tz=timezone.utc)
    utc_time_now_str = utc_time_now.isoformat(timespec="seconds")
    if args.rundate is None:
        args.rundate = utc_time_now_str
    # AR rundate correctly formatted?
    if not assert_isoformat_utc(args.rundate):
        log.error(
            "args.rundate={} is not yyyy-mm-ddThh:mm:ss+00:00; exiting".format(
                args.rundate
            )
        )
        sys.exit()

    for kwargs in args._get_kwargs():
        print(kwargs)

    # AR outdir
    if not os.path.isdir(args.outdir):
        os.mkdir(args.outdir)

    # AR subdirs
    for subdir in ["inputs", "faruns", "outputs"]:
        if not os.path.isdir(os.path.join(args.outdir, subdir)):
            os.mkdir(os.path.join(args.outdir, subdir))


    # AR saving settings
    fn = open(get_fn(args.outdir, "settings", steps=args.steps), "w")
    for kwargs in args._get_kwargs():
        fn.write("{} = {}\n".format(kwargs[0], kwargs[1]))
    fn.write("\n")
    fn.write(
        "python {} {}\n".format(
            sys.argv[0],
            " ".join(
                [
                    "--" + kwargs[0] + " " + str(kwargs[1])
                    for kwargs in args._get_kwargs()
                    if kwargs[1] is not None
                ]
            ),
        )
    )
    fn.write("\n")
    fn.close()

    # AR log file
    if args.nolog == "n":
        logfn = get_fn(args.outdir, "log", steps=args.steps)
    else:
        logfn = None

    # AR reproducible random seed
    np.random.seed(12345)

    from desitarget.targetmask import desi_mask, bgs_mask

    dtkey, btkey = "DESI_TARGET", "BGS_TARGET"

    # AR: log filename
    if args.nolog == "n":
        if os.path.isfile(logfn):
            os.remove(logfn)
        with stdouterr_redirected(to=logfn):
            main()
    else:
        main()
