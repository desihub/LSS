#copied from https://github.com/desihub/desispec/blob/master/py/desispec/tile_qa_plot.py#L957
#small updates to print basic prediction for lya numbers

from desitarget.io import read_targets_in_tiles
from desispec.maskbits import fibermask
from desispec.io import read_fibermap, findfile
from desispec.tsnr import tsnr2_to_efftime
from desimodel.focalplane.geometry import get_tile_radius_deg
from desimodel.footprint import is_point_in_desi
from desiutil.log import get_logger
from desiutil.dust import ebv as dust_ebv
from astropy.table import Table, vstack
from astropy.io import fits
from astropy import units
from astropy.coordinates import SkyCoord
tile_radius_deg = get_tile_radius_deg()

def get_quantz_cmap(name, n, cmin=0, cmax=1):
    """
    Creates a quantized colormap.
    Args:
        name: matplotlib colormap name (e.g. "tab20") (string)
        n: number of colors
        cmin (optional, defaults to 0): first color of the original colormap to use (between 0 and 1) (float)
        cmax (optional, defaults to 1): last color of the original colormap to use (between 0 and 1) (float)
    Returns:
        A matplotlib cmap object.
    Notes:
        https://matplotlib.org/examples/api/colorbar_only.html
    """
    cmaporig = matplotlib.cm.get_cmap(name)
    mycol = cmaporig(np.linspace(cmin, cmax, n))
    cmap = matplotlib.colors.ListedColormap(mycol)
    cmap.set_under(mycol[0])
    cmap.set_over (mycol[-1])
    return cmap


def get_tilecov(
    tileid,
    surveys="main",
    programs=None,
    lastnight=None,
    indesi=True,
    outpng=None,
    plot_tiles=False,
    verbose=False,
):
    """
    Computes the average number of observed tiles covering a given tile.
    Args:
        tileid: tileid (int)
        surveys (optional, defaults to "main"): comma-separated list of surveys to consider (reads the tiles-SURVEY.ecsv file) (str)
        programs (optional, defaults to None): comma-separated list of programs (case-sensitive) to consider in the tiles-SURVEY.ecsv file (str)
        lastnight (optional, defaults to today): only consider tiles observed up to lastnight (int)
        surveys (optional, defaults to "main"): comma-separated list of surveys to consider (reads the tiles-SURVEY.ecsv file) (str)
        indesi (optional, defaults to True): restrict to IN_DESI=True tiles? (bool)
        outpng (optional, defaults to None): if provided, output file with a plot (str)
        plot_tiles (optional, defaults to False): plot overlapping tiles? (bool)
        verbose (optional, defaults to False): print log.info() (bool)
    Returns:
        ntilecov: average number of observed tiles covering the considered tile (float)
        outdict: a dictionary, with an entry for each observed, overlapping tile, containing the list of observed overlapping tiles (dict)
    Notes:
        If the tile is not covered by randoms, ntilecov=np.nan, tileids=[] (and no plot is made).
        The "regular" use is to provide a single PROGRAM in programs (e.g., programs="DARK").
        This function relies on the following files:
            $DESI_SURVEYOPS/ops/tiles-{SURVEY}.ecsv for SURVEY in surveys (to get the tiles to consider)
            $DESI_ROOT/spectro/redux/daily/exposures-daily.fits (to get the existing observations up to lastnight)
            $DESI_TARGET/catalogs/dr9/2.4.0/randoms/resolve/randoms-1-0/
        If one wants to consider the latest observations, one should wait the 10am pacific update of exposures-daily.fits.
    """
    # AR lastnight
    if lastnight is None:
        lastnight = int(datetime.now().strftime("%Y%m%d"))
    # AR files
    allowed_surveys = ["sv1", "sv2", "sv3", "main", "catchall"]
    sel = ~np.in1d(surveys.split(","), allowed_surveys)
    if sel.sum() > 0:
        msg = "surveys={} not in allowed_surveys={}".format(
            ",".join([survey for survey in np.array(surveys.split(","))[sel]]),
            ",".join(allowed_surveys),
        )
        log.error(msg)
        raise ValueError(msg)
    tilesfns = [
        os.path.join(os.getenv("DESI_SURVEYOPS"), "ops", "tiles-{}.ecsv".format(survey))
        for survey in surveys.split(",")
    ]
    expsfn = os.path.join(os.getenv("DESI_ROOT"), "spectro", "redux", "daily", "exposures-daily.fits")
    # AR we need that specific version which is healpix-split, hence readable by read_targets_in_tile(quick=True))
    randdir = os.path.join(os.getenv("DESI_TARGET"), "catalogs", "dr9", "2.4.0", "randoms", "resolve", "randoms-1-0")

    # AR exposures with EFFTIME_SPEC>0 and NIGHT<=LASTNIGHT
    exps = Table.read(expsfn, "EXPOSURES")
    sel = (exps["EFFTIME_SPEC"] > 0) & (exps["NIGHT"] <= lastnight)
    exps = exps[sel]

    # AR read the tiles
    ds = []
    for tilesfn in tilesfns:
        if verbose:
            log.info("reading {}".format(tilesfn))
        d = Table.read(tilesfn)
        if d["RA"].unit == "deg":
            d["RA"].unit, d["DEC"].unit = None, None
        if "sv2" in tilesfn:
            d["IN_DESI"] = d["IN_DESI"].astype(bool)
        ds.append(d)
    tiles = vstack(ds, metadata_conflicts="silent")

    # AR first, before any cut:
    # AR - get the considered tile
    # AR - read the randoms inside that tile
    sel = tiles["TILEID"] == tileid
    if sel.sum() == 0:
        msg = "no TILEID={} found in {}".format(tileid, tilesfn)
        log.error(msg)
        raise ValueError(msg)
    if programs is None:
        log.warning("programs=None, will consider *all* kind of tiles")
    else:
        if tiles["PROGRAM"][sel][0] not in programs.split(","):
            log.warning(
                "TILEID={} has PROGRAM={}, not included in the programs={} used for computation".format(
                    tileid, tiles["PROGRAM"][sel][0], programs,
                )
            )
    c = SkyCoord(
        ra=tiles["RA"][sel][0] * units.degree,
        dec=tiles["DEC"][sel][0] * units.degree,
        frame="icrs"
    )
    d = read_targets_in_tiles(randdir, tiles=tiles[sel], quick=True)
    if len(d) == 0:
        log.warning("found 0 randoms in TILEID={}; cannot proceed; returning np.nan, empty_dictionary".format(tileid))
        return np.nan, {}
    if verbose:
        log.info("found {} randoms in TILEID={}".format(len(d), tileid))

    # AR then cut on:
    # AR - PROGRAM, IN_DESI: to get the tiles to consider
    # AR - exposures: to get the observations with NIGHT <= LASTNIGHT
    sel = np.ones(len(tiles), dtype=bool)
    if verbose:
        log.info("starting from {} tiles".format(len(tiles)))
    if programs is not None:
        sel = np.in1d(tiles["PROGRAM"], programs.split(","))
        if verbose:
            log.info("considering {} tiles after cutting on PROGRAM={}".format(sel.sum(), programs))
    if indesi:
        sel &= tiles["IN_DESI"]
        if verbose:
            log.info("considering {} tiles after cutting on IN_DESI".format(sel.sum()))
    sel &= np.in1d(tiles["TILEID"], exps["TILEID"])
    if verbose:
        log.info("considering {} tiles after cutting on NIGHT <= {}".format(sel.sum(), lastnight))
    tiles = tiles[sel]
    # AR overlap
    cs = SkyCoord(ra=tiles["RA"] * units.degree, dec=tiles["DEC"] * units.degree, frame="icrs")
    sel = cs.separation(c).value <= 2 *  tile_radius_deg
    tiles = tiles[sel]
    if verbose:
        log.info("selecting {} overlapping tiles: {}".format(len(tiles), tiles["TILEID"].tolist()))

    # AR get exposures
    outdict = {
        tileid : exps["EXPID"][exps["TILEID"] == tileid].tolist()
        for tileid in tiles["TILEID"]
    }

    # AR count the number of tile coverage
    ntile = np.zeros(len(d), dtype=int)
    for i in range(len(tiles)):
        sel = is_point_in_desi(tiles[[i]], d["RA"], d["DEC"])
        if verbose:
            log.info("fraction of TILEID={} covered by TILEID={}: {:.2f}".format(tileid, tiles[i]["TILEID"], sel.mean()))
        ntile[sel] += 1
    ntilecov = ntile.mean()
    if verbose:
        log.info("mean coverage of TILEID={}: {:.2f}".format(tileid, ntilecov))
    nlya = 0
    ntot = len(ntile)
    mod_nlya = [0,300,150,75]
    for nt in range(1,4):
        sel = ntile == nt
        nlya += len(ntile[sel])/ntot*mod_nlya[nt]
    print('predicted number of lya '+str(nlya))
    # AR plot?
    if plot_tiles:
        # AR cbar settings
        cmin = 0
        # AR for "regular" programs, setting cmax to the 
        # AR    designed max. npass (though considering future possibility
        # AR    to have more pass, e.g. for mainBRIGHT, hence the np.max())
        refcmaxs = {
            "sv3BACKUP" : 5, "sv3BRIGHT" : 11, "sv3DARK" : 14,
            "mainBACKUP" : 1, "mainBRIGHT" : 4, "mainDARK" : 7,
        }
        if "{}{}".format(surveys, programs) in refcmaxs:
            cmax = np.max([refcmaxs["{}{}".format(surveys, programs)], ntile.max()])
        else:
            cmax = ntile.max()
        cmap = get_quantz_cmap(matplotlib.cm.jet, cmax - cmin + 1, 0, 1)
        # AR case overlap Dec.=0
        if d["RA"].max() - d["RA"].min() > 100:
            dowrap = True
        else:
            dowrap = False
        if dowrap:
            d["RA"][d["RA"] > 300] -= 360
        #
        fig, ax = plt.subplots()
        sc = ax.scatter(d["RA"], d["DEC"], c=ntile, s=1, cmap=cmap, vmin=cmin, vmax=cmax)
        # AR plot overlapping tiles?
        if plot_tiles:
            angs = np.linspace(0, 2 * np.pi, 100)
            dras = tile_radius_deg * np.cos(angs)
            ddecs = tile_radius_deg * np.sin(angs)
            for i in range(len(tiles)):
                if tiles["TILEID"][i] != tileid:
                    ras = tiles["RA"][i] + dras / np.cos(np.radians(tiles["DEC"][i] + ddecs))
                    if dowrap:
                        ras[ras > 300] -= 360
                    decs = tiles["DEC"][i] + ddecs
                    ax.plot(ras, decs, label="TILEID={}".format(tiles["TILEID"][i]))
            ax.legend(loc=2, ncol=2, fontsize=8)
        #
        ax.set_title("Mean coverage of TILEID={} on {}: {:.2f}".format(tileid, lastnight, ntilecov))
        ax.set_xlabel("R.A. [deg]")
        ax.set_ylabel("Dec. [deg]")
        dra = 1.1 * tile_radius_deg / np.cos(np.radians(d["DEC"].mean()))
        ddec = 1.1 * tile_radius_deg
        ax.set_xlim(d["RA"].mean() + dra, d["RA"].mean() - dra)
        ax.set_ylim(d["DEC"].mean() - ddec, d["DEC"].mean() + ddec)
        ax.grid()
        cbar = plt.colorbar(sc, ticks=np.arange(cmin, cmax + 1, dtype=int))
        cbar.set_label("Number of observed tiles on {}".format(lastnight))
        cbar.mappable.set_clim(cmin, cmax)
        plt.show()
        #plt.savefig(outpng, bbox_inches="tight")
        #plt.close()
    return ntilecov, ntile