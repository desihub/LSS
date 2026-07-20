#!/usr/bin/env python3

from desitarget import mtl

import glob
import os
import errno
from LSS.SV3 import altmtltools as amtl
import argparse


def create_dirs(value):
    if not os.path.exists(value):
        try:
            os.makedirs(value, 0o755)
            print("made " + value)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise


def init_amtl(concate_tracers, out_dir, obscon, paral=True, recreate_tileTrack=False):
    """
    File access out parameters
    
    * aux_datapath = {
        "DARK": "/global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummit_v4_1/aux/mainsurvey-DARKobscon-TileTracker.ecsv",
        "BRIGHT": "/global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummitBGS_v2/aux_data/mainsurvey-BRIGHTobscon-TileTracker.ecsv",
        }
    * "/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-specstatus.ecsv"
    
    """
    initledger_path = os.path.join(out_dir, "initled")
    altmtl_path = os.path.join(out_dir, "Univ000")
    if not os.path.isdir(initledger_path):
        print("Running initial ledgers")
        if paral:
            mtl.make_ledger(concate_tracers, initledger_path, obscon=obscon.upper(), numproc=4)
        else:
            mtl.make_ledger(concate_tracers, initledger_path, obscon=obscon.upper())
    print("Creating list of tiles to be processed by AltMTL mock production")
    path = os.path.join(initledger_path, "main", obscon.lower())
    ff = glob.glob(
        os.path.join(path, "mtl-{obscon}-hp-*.ecsv".format(obscon=obscon.lower()))
    )
    dd = []
    for f in ff:
        dd.append(int(f.split("hp-")[-1].split(".ecsv")[0]))
    tosave = ",".join(map(str, sorted(dd)))
    savepath = os.path.join(
        initledger_path, "hpxlist_{obscon}.txt".format(obscon=obscon.lower())
    )
    ff = open(savepath, "w")
    ff.write(tosave)
    ff.close()
    print("saving list of HP ledgers in " + savepath)
    path_to_altmtl = os.path.join(altmtl_path, "main", obscon.lower())
    print("Copying initial ledgers to altmtl directory ", path_to_altmtl)
    create_dirs(path_to_altmtl)
    os.system("cp %s %s" % (os.path.join(path, "*"), path_to_altmtl))
    print("Creating tileTracker file and tilestatus file")
    startDateShort = 19990101
    endDate = "20240418"  # 2024-04-18T00:00:00+00:00
    if ("T" in endDate) & ("-" in endDate):
        endDateShort = int(endDate.split("T")[0].replace("-", ""))
    else:
        endDateShort = int(endDate)

    aux_datapath = {
        "DARK": "/global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummit_v4_1/aux/mainsurvey-DARKobscon-TileTracker.ecsv",
        "BRIGHT": "/global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummitBGS_v2/aux_data/mainsurvey-BRIGHTobscon-TileTracker.ecsv",
    }
    if recreate_tileTrack:
        amtl.makeTileTracker(
            altmtl_path,
            survey="main",
            obscon=obscon.upper(),
            startDate=startDateShort,
            endDate=endDateShort,
            overwrite=True,
        )
    else:
        os.system("cp %s %s" % (aux_datapath[obscon.upper()], altmtl_path))
    ztilefile = (
        "/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-specstatus.ecsv"
    )
    ztilefn = ztilefile.split("/")[-1]
    if not os.path.isfile(os.path.join(altmtl_path, ztilefn)):
        amtl.processTileFile(
            ztilefile, os.path.join(altmtl_path, ztilefn), None, endDate
        )


#
# Main
#
# initialize_amtl_mocks_da2.py
#    /pscratch/sd/d/desica/DA2/mocks/holi_v3/forFA$i.fits
#    /pscratch/sd/d/desica/DA2/mocks/holi_v3/altmtl$i/
#    DARK

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    # file concatenated forFA0.fits
    parser.add_argument("--inputs", nargs="+", required=True)
    # directory to save altmtl and initledgers
    parser.add_argument("--outputs", nargs="+", required=True)
    parser.add_argument('--noskip', action='store_true', help='Do not skip existing files')
    parser.add_argument("--obscon", choices=["DARK", "BRIGHT"], required=True)
    args = parser.parse_args()

    print("\n\nInit amtl for mock catalogs")
    print(f"inputs : {args.inputs}")
    print(f"outputs: {args.outputs}")

    concate_tracers = args.inputs[0]  # Input mock
    out_dir = args.outputs[0]  # Output path
    obscon = args.obscon  # DARK or BRIGHT

    init_amtl(concate_tracers, out_dir, obscon)
