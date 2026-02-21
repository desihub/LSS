#!/bin/bash

set -e

source /global/common/software/desi/desi_environment.sh main
export LSSCODE=$HOME
PYTHONPATH=$PYTHONPATH:$LSSCODE/LSS

TRACER='ELG_LOP'

srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_ran.py --type $TRACER --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec iron --survey Y1 --add_tl y --apply_veto y --maxr 18 --version $1

python $LSSCODE/LSS/scripts/main/mkCat_main.py --type $TRACER --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec iron --survey Y1  --mkHPmaps y --version $1

python $LSSCODE/LSS/scripts/main/mkCat_main.py --type $TRACER --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --ranonly y --apply_map_veto y --imsys y --survey Y1 --verspec iron --imsys_zbin y --use_map_veto _HPmapcut --version $1

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi scripts/Y1_sysnetELG_zbins_new.sh $1 #note, this loads cosmodesi to run

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/main/mkCat_main.py --type $TRACER --notqso y  --fulld n --survey Y1 --verspec iron --version $1 --clusd y --clusran y --splitGC y --nz y --par y --imsys_colname WEIGHT_SN --basedir /global/cfs/cdirs/desi/survey/catalogs/

