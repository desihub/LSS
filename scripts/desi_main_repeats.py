#!/usr/bin/env python

import os
import multiprocessing
import itertools
import fitsio
import numpy as np
from astropy.io import fits
from astropy import constants
from astropy.table import Table, vstack, hstack
from desitarget.targetmask import desi_mask, bgs_mask
from desitarget.sv3.sv3_targetmask import desi_mask as sv3_desi_mask
from desitarget.sv3.sv3_targetmask import bgs_mask as sv3_bgs_mask
from desitarget.targetmask import zwarn_mask as zmtl_zwarn_mask
from desitarget.geomask import match_to
from desispec.validredshifts import validate
from matplotlib import pyplot as plt
from matplotlib import gridspec
from argparse import ArgumentParser

allowed_prods = ["iron", "jura", "kibo"]
allowed_steps = ["parent", "pairs", "plot"]


def parse():
    parser = ArgumentParser()
    parser.add_argument(
        "--outroot",
        help="output root",
        required=True,
        type=str,
        default=None,
    )
    parser.add_argument(
        "--prod",
        help="spectro prod (default=kibo)",
        choices=allowed_prods,
        type=str,
        default="kibo",
    )
    parser.add_argument(
        "--prog",
        help="desi program (default=dark)",
        choices=["bright", "dark"],
        type=str,
        default="dark",
    )
    parser.add_argument(
        "--steps",
        help="comma-separated list of steps (default={})".format(
            ",".join(allowed_steps)
        ),
        type=str,
        default=",".join(allowed_steps),
    )
    parser.add_argument(
        "--numproc",
        help="number of concurrent processes to use; set to 0 to not process (default=1)",
        type=int,
        default=1,
    )
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    for kwargs in args._get_kwargs():
        print(kwargs)
    return args


def get_prod_dirs(prod):

    prod_dir = os.path.join(os.getenv("DESI_ROOT"), "spectro", "redux", prod)
    prod_dir = prod_dir.replace("/global/", "/dvs_ro/")

    if prod in ["iron", "jura", "kibo"]:

        zcat_dir = os.path.join(prod_dir, "zcatalog", "v1")

    return prod_dir, zcat_dir


def get_snr2time(tracer):

    ensemble_dir = os.path.join(os.getenv("DESIMODEL"), "data", "tsnr")
    fn = os.path.join(ensemble_dir, "tsnr-ensemble-{}.fits".format(tracer.lower()))

    return fits.getheader(fn, 0)["SNR2TIME"]


def get_prog_props(prog):

    survey_dict = {
        "bright": {
            "tracers": ["BGS_BRIGHT", "BGS_FAINT"],
            "sv3mask": sv3_bgs_mask,
            "mainmask": bgs_mask,
            "sv3key": "SV3_BGS_TARGET",
            "mainkey": "BGS_TARGET",
        },
        "dark": {
            "tracers": ["LRG", "ELG", "QSO"],
            "sv3mask": sv3_desi_mask,
            "mainmask": desi_mask,
            "sv3key": "SV3_DESI_TARGET",
            "mainkey": "DESI_TARGET",
        },
    }

    return survey_dict[prog]


# AR main: restrict to:
# AR - tracers
def read_main(zcat_dir, prog, keys):

    prog_props = get_prog_props(prog)
    tracers = prog_props["tracers"]
    mask = prog_props["mainmask"]
    mask_key = prog_props["mainkey"]

    fn = os.path.join(zcat_dir, "ztile-main-{}-cumulative.fits".format(prog))
    main_keys = ["DESI_TARGET", "BGS_TARGET"]
    d = Table(fitsio.read(fn, "ZCATALOG", columns=keys + main_keys))
    print(len(d))
    d = d[keys + main_keys]  # AR re-ordering columns

    d["SURVEY"] = "main"

    sel = np.zeros(len(d), dtype=bool)
    for tracer in tracers:
        sel |= (d[mask_key] & mask[tracer]) > 0
    d = d[sel]

    return d


# AR sv3: restrict to
# AR - tracers
# AR - first-time observations
# AR - keep a unique TARGETID (<0.02% repeats remaining...), with the highest TSNR2_LRG
def read_sv3(zcat_dir, prog, keys):

    prog_props = get_prog_props(prog)
    tracers = prog_props["tracers"]
    mask = prog_props["sv3mask"]
    mask_key = prog_props["sv3key"]

    fn = os.path.join(zcat_dir, "ztile-sv3-{}-cumulative.fits".format(prog))
    sv3_keys = ["SV3_DESI_TARGET", "SV3_BGS_TARGET", "PRIORITY_INIT", "PRIORITY"]
    d = Table(fitsio.read(fn, "ZCATALOG", columns=keys + sv3_keys))
    d = d[keys + sv3_keys]  # AR re-ordering columns
    d["SURVEY"] = "sv3"

    sel = np.zeros(len(d), dtype=bool)
    for tracer in tracers:
        sel |= (d[mask_key] & mask[tracer]) > 0
    sel &= d["PRIORITY"] == d["PRIORITY_INIT"]
    d = d[sel]

    # AR unique TARGETID
    if prog == "bright":
        ii = np.lexsort([-d["TSNR2_BGS"], d["TARGETID"]])
    if prog == "dark":
        ii = np.lexsort([-d["TSNR2_LRG"], d["TARGETID"]])
    d = d[ii]
    _, ii = np.unique(d["TARGETID"], return_index=True)
    d = d[ii]

    return d


def get_validredshifts_zmtl(rrfn):

    extra_columns = [
        "TARGETID",
        "OII_FLUX",
        "OII_FLUX_IVAR",
        "Z_NEW",
        "ZERR_NEW",
        "IS_QSO_QN_NEW_RR",
    ]

    d = validate(
        rrfn,
        fiberstatus_cut=False,
        return_target_columns=False,
        extra_columns=extra_columns,
    )

    d["RRFN"] = rrfn

    # AR add zmtl columns
    zmtlfn = rrfn.replace("redrock-", "zmtl-")
    d2 = fitsio.read(zmtlfn, "ZMTL")
    assert np.all(d["TARGETID"] == d2["TARGETID"])
    for key in ["ZWARN", "Z_QN", "Z_QN_CONF", "IS_QSO_QN"]:
        d["ZMTL_{}".format(key)] = d2[key]

    return d


def get_valid_zmtl(rrfns, numproc):

    pool = multiprocessing.Pool(numproc)
    with pool:
        ds = pool.map(get_validredshifts_zmtl, rrfns)

    d = vstack(ds)

    return d


def create_parent(prod, prog, numproc):

    prod_dir, zcat_dir = get_prod_dirs(prod)
    keys = [
        "TARGETID",
        "TARGET_RA",
        "TARGET_DEC",
        "PHOTSYS",
        "COADD_FIBERSTATUS",
        "TILEID",
        "LASTNIGHT",
        "PETAL_LOC",
        "FIBER",
        "TSNR2_BGS",
        "TSNR2_LRG",
        "TSNR2_ELG",
        "TSNR2_QSO",
        "TSNR2_LYA",
        "Z",
        "ZERR",
        "ZWARN",
        "DELTACHI2",
    ]

    # AR read
    main_d = read_main(zcat_dir, prog, keys)
    sv3_d = read_sv3(zcat_dir, prog, keys)

    # AR now, add sv3 repeats of main, then merge
    sel = np.in1d(sv3_d["TARGETID"], main_d["TARGETID"])
    sv3_d = sv3_d[sel]
    sv3_d["DESI_TARGET"], sv3_d["BGS_TARGET"] = 0, 0
    sv3_d = sv3_d[main_d.colnames]
    assert np.all(main_d.colnames == sv3_d.colnames)

    # AR fill in DESI_TARGET, BGS_TARGET
    _, ii = np.unique(main_d["TARGETID"], return_index=True)
    tmp_d = main_d[ii]
    ii = match_to(tmp_d["TARGETID"], sv3_d["TARGETID"])
    assert np.all(tmp_d["TARGETID"][ii] == sv3_d["TARGETID"])
    for key in ["DESI_TARGET", "BGS_TARGET"]:
        sv3_d[key] = tmp_d[key][ii]

    # AR vstack
    d = vstack([main_d, sv3_d])

    # AR restrict to repeats
    tids, counts = np.unique(d["TARGETID"], return_counts=True)
    sel = np.in1d(d["TARGETID"], tids[counts > 1])
    d = d[sel]

    # AR redrock path
    d["RRFN"] = [
        os.path.join(
            prod_dir,
            "tiles",
            "cumulative",
            str(tileid),
            str(lastnight),
            "redrock-{}-{}-thru{}.fits".format(petal, tileid, lastnight),
        )
        for tileid, lastnight, petal in zip(d["TILEID"], d["LASTNIGHT"], d["PETAL_LOC"])
    ]

    # AR valid redshifts + zmtl_zwarn
    rrfns = np.unique(d["RRFN"])
    print(len(rrfns))
    valid_d = get_valid_zmtl(rrfns, numproc)
    print(len(valid_d))

    # AR speed up
    sel = np.in1d(valid_d["TARGETID"], d["TARGETID"])
    valid_d = valid_d[sel]
    print(len(valid_d))

    # AR match
    unqids = np.array(
        ["{}-{}".format(tid, rrfn) for tid, rrfn in zip(d["TARGETID"], d["RRFN"])]
    )
    valid_unqids = np.array(
        [
            "{}-{}".format(tid, rrfn)
            for tid, rrfn in zip(valid_d["TARGETID"], valid_d["RRFN"])
        ]
    )
    ii = match_to(valid_unqids, unqids)
    assert np.all(valid_unqids[ii] == unqids)
    for key in valid_d.colnames:
        if key not in d.colnames:
            d[key] = valid_d[key][ii]


    # AR header infos
    d.meta["PRODDIR"], d.meta["ZCATDIR"] = prod_dir, zcat_dir
    d.meta["PROG"] = prog
    for tracer in ["BGS", "LRG", "ELG", "QSO", "LYA"]:
        d.meta["{}SNR2T".format(tracer)] = get_snr2time(tracer)

    return d


def get_combinations(ii):

    jjs = itertools.combinations(ii, 2)
    j0s, j1s = [], []
    for jj in jjs:
        j0s.append(jj[0])
        j1s.append(jj[1])

    return j0s, j1s


def create_pairs(parent_d, numproc):

    # AR keys in common
    comm_keys = [
        "TARGETID",
        "TARGET_RA",
        "TARGET_DEC",
        "PHOTSYS",
        "DESI_TARGET",
        "BGS_TARGET",
    ]

    # AR keys specific to each spectrum
    rp_keys = [
        "SURVEY",
        "TILEID",
        "LASTNIGHT",
        "PETAL_LOC",
        "FIBER",
        "COADD_FIBERSTATUS",
        "TSNR2_BGS",
        "TSNR2_LRG",
        "TSNR2_ELG",
        "TSNR2_QSO",
        "TSNR2_LYA",
        "Z",
        "ZERR",
        "ZWARN",
        "DELTACHI2",
        "OII_FLUX",
        "OII_FLUX_IVAR",
        "Z_NEW",
        "ZERR_NEW",
        "IS_QSO_QN_NEW_RR",
        "ZMTL_ZWARN",
        "ZMTL_Z_QN",
        "ZMTL_Z_QN_CONF",
        "ZMTL_IS_QSO_QN",
        "GOOD_BGS",
        "GOOD_LRG",
        "GOOD_ELG",
        "GOOD_QSO",
    ]

    _, counts = np.unique(parent_d["TARGETID"], return_counts=True)
    assert counts.min() > 1
    print(
        "using  {} unique TARGETIDs to build pairs from {} observations".format(
            counts.size, len(parent_d)
        )
    )

    # AR sort by TARGETID
    ii = parent_d["TARGETID"].argsort()
    parent_d = parent_d[ii]

    # AR ii: index of first occurence of a TARGETID
    tids, ii = np.unique(parent_d["TARGETID"], return_index=True)

    # AR n_repeats: nb of occurences for a TARGETID
    n_repeats = np.diff(ii)
    n_repeats = np.append(n_repeats, len(parent_d) - ii[-1])

    # AR list of lists of indexes for each tid
    iis = [[i + _ for _ in range(n_repeat)] for i, n_repeat in zip(ii, n_repeats)]

    # AR get all the pair combinations
    pool = multiprocessing.Pool(numproc)
    with pool:
        outs = pool.map(get_combinations, iis)
    ii0 = np.hstack([_[0] for _ in outs])
    ii1 = np.hstack([_[1] for _ in outs])

    # AR build two arrays
    d0, d1 = parent_d[ii0], parent_d[ii1]
    assert np.all(d0["TARGETID"] == d1["TARGETID"])
    assert np.all(d0["TILEID"] != d1["TILEID"])

    # AR hstack
    d0 = d0[comm_keys + rp_keys]
    d1 = d1[rp_keys]
    for key in rp_keys:
        d0[key].name = "{}_0".format(key)
        d1[key].name = "{}_1".format(key)
    d = hstack([d0, d1])

    # AR dv
    d["DV"] = constants.c.to("km/s").value * (d["Z_1"] - d["Z_0"]) / (1 + d["Z_0"])
    d["DV_NEW"] = (
        constants.c.to("km/s").value
        * (d["Z_NEW_1"] - d["Z_NEW_0"])
        / (1 + d["Z_NEW_0"])
    )

    return d


# TODO: add qsos...
def make_plot(outroot, d):

    prog, prod_dir = d.meta["PROG"], d.meta["PRODDIR"]
    prod = os.path.basename(prod_dir)

    if prog == "bright":
        tracers = ["BGS_BRIGHT", "BGS_FAINT"]
        mask, mask_key = bgs_mask, "BGS_TARGET"
        effkey, effmin, effmax, effxlim = "TSNR2_BGS", 0.85 * 180, 1.5 * 180, (0, 500)
    if prog == "dark":
        tracers = ["LRG", "ELG_LOP", "ELG_VLO"]
        mask, mask_key = desi_mask, "DESI_TARGET"
        effkey, effmin, effmax, effxlim = (
            "TSNR2_LRG",
            0.85 * 1000,
            1.5 * 1000,
            (500, 1500),
        )

    for tracer in tracers:

        # AR per-tracer setting
        if tracer[:3] == "BGS":
            zmin, zmax, goodkey = 0.1, 0.4, "GOOD_BGS"
        else:
            if tracer == "LRG":
                zmin, zmax, goodkey = 0.4, 1.1, "GOOD_LRG"
            else:
                zmin, zmax, goodkey = 0.6, 1.6, "GOOD_ELG"
        # AR
        if tracer == "QSO":
            dv_thresh = 3000
        else:
            dv_thresh = 1000
        # AR
        if tracer[:3] == "ELG":
            xs = np.log10(
                d["OII_FLUX_0"] * np.sqrt(d["OII_FLUX_IVAR_0"])
            ) + 0.2 * np.log10(d["DELTACHI2_0"])
            xlim = (0, 3)
            xlab = "log10(FOII_SNR) + 0.2 * log10(DELTACHI2)"
        else:
            xs = np.log10(d["DELTACHI2_0"])
            xlab = "log10(DELTACHI2)"
            xlim = (0, 6)

        # AR efftime_spec
        snr2time = d.meta["{}SNR2T".format(effkey.split("_")[1])]
        efftime0s = snr2time * d["{}_0".format(effkey)]
        efftime1s = snr2time * d["{}_1".format(effkey)]

        # AR for both elements of a pair, cut on:
        # AR - tracer
        # AR - coadd_fiberstatus
        # AR - zmtl_zwarn_mask nodata + bad
        # AR - efftime_spec
        # AR cut on redshift range with a OR
        sel = (d[mask_key] & mask[tracer]) > 0
        sel &= (d["COADD_FIBERSTATUS_0"] == 0) & (d["COADD_FIBERSTATUS_1"] == 0)
        nodata0 = (d["ZMTL_ZWARN_0"] & zmtl_zwarn_mask["NODATA"]) > 0
        nodata1 = (d["ZMTL_ZWARN_1"] & zmtl_zwarn_mask["NODATA"]) > 0
        badqa0 = (d["ZMTL_ZWARN_0"] & zmtl_zwarn_mask.mask("BAD_SPECQA|BAD_PETALQA")) > 0
        badqa1 = (d["ZMTL_ZWARN_1"] & zmtl_zwarn_mask.mask("BAD_SPECQA|BAD_PETALQA")) > 0
        sel &= (~nodata0) & (~nodata1)
        sel &= (~badqa0) & (~badqa1)
        sel &= (efftime0s > effmin) & (efftime1s > effmin)
        sel &= (efftime0s < effmax) & (efftime1s < effmax)
        sel &= (d["{}_0".format(goodkey)]) & (d["{}_1".format(goodkey)])
        selz0 = (d["Z_0"] > zmin) & (d["Z_0"] < zmax)
        selz1 = (d["Z_1"] > zmin) & (d["Z_1"] < zmax)
        sel &= (selz0) | (selz1)
        fail = (sel) & (np.abs(d["DV"]) > dv_thresh)

        # AR plot
        fig = plt.figure(figsize=(30, 5))
        gs = gridspec.GridSpec(1, 5, wspace=0.3)

        title = "Main / {} {} : {}".format(prog, prod, tracer)
        sellab = "All selected pairs ({})".format(sel.sum())
        faillab = "Failures ({} / {} = {:.2f}%".format(
            fail.sum(), sel.sum(), 100.0 * fail.sum() / sel.sum()
        )

        # AR sky
        d["TARGET_RA"][d["TARGET_RA"] > 300] -= 360
        ax = plt.subplot(gs[0])
        ax.scatter(
            d["TARGET_RA"][sel],
            d["TARGET_DEC"][sel],
            c="0.5",
            s=1,
            alpha=0.1,
            label=sellab,
            rasterized=True,
        )
        ax.scatter(
            d["TARGET_RA"][fail],
            d["TARGET_DEC"][fail],
            c="r",
            s=1,
            alpha=0.5,
            label=faillab,
            rasterized=True,
        )
        ax.set_title(title)
        ax.set_xlabel("R.A. [deg]")
        ax.set_ylabel("Dec. [deg]")
        ax.set_xlim(300, -60)
        ax.set_ylim(-30, 90)
        ax.grid()
        ax.legend(loc=1, markerscale=5)

        # AR efftime_spec
        ax = plt.subplot(gs[1])
        bins = np.linspace(0, 2000, 500)
        _ = ax.hist(
            efftime0s[sel],
            bins=bins,
            density=True,
            histtype="stepfilled",
            color="0.5",
            alpha=0.5,
            label=sellab,
        )
        _ = ax.hist(
            efftime0s[fail],
            bins=bins,
            density=True,
            histtype="step",
            color="r",
            alpha=1.0,
            label=faillab,
        )
        for efftime in [effmin, effmax]:
            ax.axvline(efftime, color="k", lw=2, ls="--", label="{}s".format(efftime))
        ax.set_title(
            "Main/{} {} : {} (failures: {} / {} = {:.2f}%)".format(
                prog,
                prod,
                tracer,
                fail.sum(),
                sel.sum(),
                100.0 * fail.sum() / sel.sum(),
            )
        )
        ax.set_title(title)
        ax.set_xlabel("EFFTIME_SPEC [s]")
        ax.set_ylabel("Norm. counts")
        ax.set_xlim(effxlim)
        ax.grid()
        ax.legend(loc=2)

        # AR z
        ax = plt.subplot(gs[2])
        bins = np.linspace(0, 2, 200)
        _ = ax.hist(
            d["Z_0"][sel],
            bins=bins,
            density=True,
            histtype="stepfilled",
            color="0.5",
            alpha=0.5,
            label=sellab,
        )
        _ = ax.hist(
            d["Z_0"][fail],
            bins=bins,
            density=True,
            histtype="step",
            color="r",
            alpha=1.0,
            label=faillab,
        )
        for z in [zmin, zmax]:
            ax.axvline(z, color="k", lw=2, ls="--", label=str(z))
        ax.set_title(
            "Main/{} {} : {} (failures: {} / {} = {:.2f}%)".format(
                prog,
                prod,
                tracer,
                fail.sum(),
                sel.sum(),
                100.0 * fail.sum() / sel.sum(),
            )
        )
        ax.set_title(title)
        ax.set_xlabel("Z")
        ax.set_ylabel("Norm. counts")
        ax.set_xlim(bins[0], bins[-1])
        ax.grid()
        ax.legend(loc=1)

        # AR dchi2
        ax = plt.subplot(gs[3])
        ax.scatter(
            xs[sel],
            d["DV"][sel],
            c="0.5",
            s=1,
            alpha=0.1,
            rasterized=True,
            label=sellab,
        )
        ax.scatter(
            xs[fail],
            d["DV"][fail],
            c="r",
            s=1,
            alpha=0.5,
            rasterized=True,
            label=faillab,
        )
        ax.set_title(title)
        ax.set_xlabel(xlab)
        ax.set_ylabel("|Dv| [km/s]")
        ax.set_yscale("log")
        ax.set_xlim(xlim)
        ax.set_ylim(1e-3, 1e6)
        ax.grid()
        ax.axhline(
            dv_thresh, color="k", lw=2, ls="--", label="{} km/s".format(dv_thresh)
        )
        ax.legend(loc=3, markerscale=5)

        # AR dv
        ax = plt.subplot(gs[4])
        bins = np.linspace(-200, 200, 1000)
        _ = ax.hist(
            d["DV"][sel],
            bins=bins,
            density=True,
            histtype="stepfilled",
            color="0.5",
            alpha=0.5,
            label=sellab,
        )
        ax.set_title(title)
        ax.set_xlabel("Dv [km/s]")
        ax.set_ylabel("Norm. counts")
        ax.set_xlim(bins[0], bins[-1])
        ax.set_ylim(0, 0.05)
        ax.grid()
        ax.legend(loc=2)

        plt.savefig("{}-{}.png".format(outroot, tracer), bbox_inches="tight")
        plt.close()


def main():

    args = parse()

    if "parent" in args.steps:

        parent_d = create_parent(args.prod, args.prog, args.numproc)
        parent_d.write("{}-parent.fits".format(args.outroot), overwrite=args.overwrite)

    if "pairs" in args.steps:

        parent_d = Table.read("{}-parent.fits".format(args.outroot))
        pairs_d = create_pairs(parent_d, args.numproc)
        for key in parent_d.meta:
            pairs_d.meta[key] = parent_d.meta[key]
        pairs_d.write("{}-pairs.fits".format(args.outroot), overwrite=args.overwrite)

    if "plot" in args.steps:

        d = Table.read("{}-pairs.fits".format(args.outroot))
        make_plot("{}-pairs".format(args.outroot), d)


if __name__ == "__main__":
    main()
