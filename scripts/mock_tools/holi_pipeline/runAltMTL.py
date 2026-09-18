#!/usr/bin/env python

import argparse

from LSS.SV3 import altmtltools as amt

from desiutil.log import get_logger


log = get_logger()


parser = argparse.ArgumentParser(
    prog="runAltMTL",
    description="Progresses alternate MTLs through the MTL.",
)

parser.add_argument(
    "-a",
    "--altMTLBaseDir",
    dest="altMTLBaseDir",
    required=True,
    type=str,
    help="the path to the location where alt MTLs are stored, up to, but not including survey and obscon information.",
)

parser.add_argument(
    "-obscon",
    "--obscon",
    dest="obscon",
    default="DARK",
    help="observation conditions, either BRIGHT or DARK.",
    required=False,
    type=str,
)
parser.add_argument(
    "-mockid",
    "--mockid",
    dest="mockid",
    help="ID of the mock/seed",
    required=True,
    type=int,
)
parser.add_argument(
    "-s",
    "--survey",
    dest="survey",
    default="sv3",
    help="DESI survey to create Alt MTLs for. Either sv3 or main.",
    required=False,
    type=str,
)
parser.add_argument(
    "-sec",
    "--secondary",
    dest="secondary",
    default=False,
    action="store_true",
    help="set flag to incorporate secondary targets.",
)
parser.add_argument(
    "-mock",
    "--mock",
    dest="mock",
    default=True,
    action="store_true",
    help="set flag if running pipeline on mocks.",
)
parser.add_argument(
    "-tf",
    "--targfile",
    dest="targfile",
    required=False,
    default=None,
    type=str,
    help="Location for target file for mocks or data. Only required if mocks are being processed.",
)
parser.add_argument(
    "-v",
    "--verbose",
    dest="verbose",
    default=False,
    action="store_true",
    help="set flag to enter verbose mode",
)
parser.add_argument(
    "-qr",
    "--quickRestart",
    dest="quickRestart",
    default=False,
    action="store_true",
    help="set flag to remove any AMTL updates that have already been performed. Useful for rapidfire debugging of steps in this part of the pipeline.",
)
parser.add_argument(
    "-rep",
    "--reproducing",
    action="store_true",
    dest="reproducing",
    default=False,
    help="WARNING: THIS FLAG SHOULD ONLY BE USED FOR DEBUGGING. Pass this flag to confirm to the alt mtl code that you are trying to reproduce real MTLs. This option should (must?) be used in conjunction with --shuffleSubpriorities.",
    required=False,
)
parser.add_argument(
    "-prof",
    "--profile",
    dest="profile",
    default=False,
    action="store_true",
    help="set flag to profile code time usage. This flag may not profile all components of any particular stage of the AMTL pipeline. ",
)
parser.add_argument(
    "-d",
    "--debug",
    dest="debug",
    default=False,
    action="store_true",
    help="set flag to enter debug mode.",
)
parser.add_argument(
    "-nfl",
    "--NumObsNotFromLedger",
    dest="numobs_from_ledger",
    default=True,
    action="store_false",
    help="If True (flag is NOT set) then inherit the number of observations so far from the ledger rather than expecting it to have a reasonable value in the zcat.",
)

parser.add_argument(
    "-redoFA",
    "--redoFA",
    dest="redoFA",
    default=False,
    action="store_true",
    help="pass this flag to regenerate already existing fiber assignment files.",
)

parser.add_argument(
    "-getosubp",
    "--getosubp",
    action="store_true",
    dest="getosubp",
    default=False,
    help="WARNING: THIS FLAG SHOULD ONLY BE USED FOR DEBUGGING AND NEVER FOR MOCKS. Pass this flag to grab subpriorities directly from the real survey MTLs for fiberassignment.",
    required=False,
)
parser.add_argument(
    "-md",
    "--multiDate",
    action="store_true",
    dest="multiDate",
    default=False,
    help="Currently this flag is being debugged. In the future, it will switch between interactive submission of each date as a separate job (True) and of all nights to be looped through until a single job`s time runs out. ",
    required=False,
)
parser.add_argument(
    "-ppn",
    "--ProcPerNode",
    dest="ProcPerNode",
    default=None,
    help="Number of processes to spawn per requested node. If not specified, determined automatically from NERSC_HOST.",
    required=False,
    type=int,
)
parser.add_argument(
    "-rmbd",
    "--realMTLBaseDir",
    dest="mtldir",
    default="/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/",
    help="Location of the real (or mock) MTLs that serve as the basis for the alternate MTLs. Defaults to location of data MTLs. Do NOT include survey or obscon information here. ",
    required=False,
    type=str,
)
parser.add_argument(
    "-zcd",
    "--zCatDir",
    dest="zcatdir",
    default="/global/cfs/cdirs/desi/spectro/redux/daily/",
    help="Location of the real redshift catalogs for use in alt MTL loop.  Defaults to location of survey zcatalogs.",
    required=False,
    type=str,
)
parser.add_argument(
    "-zfix",
    "--zfix",
    dest="zfix",
    required=False,
    default=None,
    type=str,
    help="Filename with redshifts to fix the update_ledger altZCat",
)

args = parser.parse_args()

if args.mock:
    log.info("args.getosubp: {args.getosubp}")
    assert not (args.getosubp)

# Leave confirmation file in output directory if using original subpriorities
if args.getosubp:
    from pathlib import Path

    Path(args.altMTLBaseDir + "/GETOSUBPTRUE").touch()


singleDate = not (args.multiDate)

if args.zfix is None:
    zfix = None
else:
    zfix = args.zfix.format(mock_number=args.mockid)

retval = amt.loop_alt_ledger(
    args.obscon,
    survey=args.survey,
    mtldir=args.mtldir,
    zcatdir=args.zcatdir,
    altmtlbasedir=args.altMTLBaseDir.format(mock_number=args.mockid),
    ndirs=None,
    numobs_from_ledger=args.numobs_from_ledger,
    secondary=args.secondary,
    getosubp=args.getosubp,
    quickRestart=args.quickRestart,
    multiproc=False,
    nproc=args.mockid,
    singleDate=singleDate,
    redoFA=args.redoFA,
    mock=args.mock,
    targets=None,
    debug=args.debug,
    verbose=args.verbose,
    reproducing=args.reproducing,
    debugOrig=True,
    zfix=zfix,
)

if args.verbose:
    log.debug("retval: {0}".format(retval))
    log.debug("finished with one iteration of procFunc")
