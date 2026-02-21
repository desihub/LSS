#!/bin/bash

set -e

source /global/common/software/desi/desi_environment.sh main
export LSSCODE=$HOME
PYTHONPATH=$PYTHONPATH:$LSSCODE/LSS

TRACER='QSO'

srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_ran.py --type $TRACER --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec iron --survey Y1 --add_tl y --apply_veto y --maxr 18 --version $1

python $LSSCODE/LSS/scripts/main/mkCat_main.py --type $TRACER --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec iron --survey Y1  --mkHPmaps y --version $1

python $LSSCODE/LSS/scripts/main/mkCat_main.py --type $TRACER --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --ranonly y --apply_map_veto y --imsys y --survey Y1 --verspec iron --imsys_zbin y --use_map_veto _HPmapcut --version $1

source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
export LSSCODE=$HOME
PYTHONPATH=$PYTHONPATH:$LSSCODE/LSS

python scripts/main/mkCat_main.py --type $TRACER --version $1  --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --regressis y --add_regressis y --survey Y1 --verspec iron --imsys_zbin y

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/main/mkCat_main.py --type $TRACER  --fulld n --survey Y1 --verspec iron --version $1 --clusd y --clusran y --splitGC y --nz y --par y --imsys_colname WEIGHT_RF--basedir /global/cfs/cdirs/desi/survey/catalogs/

