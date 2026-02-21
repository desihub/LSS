#!/bin/bash

set -e

source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
export LSSCODE=$HOME
PYTHONPATH=$PYTHONPATH:$LSSCODE/LSS

python $LSSCODE/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT --verspec iron --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --add_fs y --survey Y1 --version $1

python $LSSCODE/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT-21.5 --verspec iron --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --survey Y1 --version $1 --redoBGS215 y

python $LSSCODE/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT-21.5  --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --imsys y --survey Y1 --verspec iron --imsys_zbin y --version $1 --use_map_veto _HPmapcut

python $LSSCODE/LSS/scripts/validation/validation_improp_full.py --tracers BGS_BRIGHT-21.5 --version $1 --weight_col WEIGHT_IMLIN

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSSCODE/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT-21.5  --fulld n --survey Y1 --verspec iron --version $1 --clusd y --clusran y --splitGC y --nz y --imsys_colname WEIGHT_IMLIN --par y --basedir /global/cfs/cdirs/desi/survey/catalogs/

