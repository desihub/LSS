#!/bin/bash

set -e

source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
export LSSCODE=$HOME
PYTHONPATH=$PYTHONPATH:$LSSCODE/LSS

#srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/main/mkCat_main.py --type ELG_LOP --notqso y --fulld n --survey Y1 --verspec iron --version $1 --clusd y --clusran y --splitGC y --nz y --par y --imsys_colname WEIGHT_SN --basedir /global/cfs/cdirs/desi/survey/catalogs/ 

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/main/mkCat_main.py --type ELG_LOP --notqso y --fulld n --survey Y1 --verspec iron --version ${1}pip --clusd y --clusran y --splitGC y --nz y --par y --imsys_colname WEIGHT_SN --basedir /global/cfs/cdirs/desi/survey/catalogs/ --compmd altmtl

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/main/mkCat_main.py --type LRG --fulld n --survey Y1 --verspec iron --version $1 --clusd y --clusran y --splitGC y --nz y --par y --imsys_colname WEIGHT_IMLIN --basedir /global/cfs/cdirs/desi/survey/catalogs/ 

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/main/mkCat_main.py --type LRG --fulld n --survey Y1 --verspec iron --version ${1}pip --clusd y --clusran y --splitGC y --nz y --par y --imsys_colname WEIGHT_IMLIN --basedir /global/cfs/cdirs/desi/survey/catalogs/ --compmd altmtl

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/main/mkCat_main.py --type BGS_BRIGHT-21.5 --fulld n --survey Y1 --verspec iron --version $1 --clusd y --clusran y --splitGC y --nz y --par y --imsys_colname WEIGHT_IMLIN --basedir /global/cfs/cdirs/desi/survey/catalogs/ 

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/main/mkCat_main.py --type BGS_BRIGHT-21.5 --fulld n --survey Y1 --verspec iron --version ${1}pip --clusd y --clusran y --splitGC y --nz y --par y --imsys_colname WEIGHT_IMLIN --basedir /global/cfs/cdirs/desi/survey/catalogs/ --compmd altmtl

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/main/mkCat_main.py --type QSO --fulld n --survey Y1 --verspec iron --version $1 --clusd y --clusran y --splitGC y --nz y --par y --imsys_colname WEIGHT_RF --basedir /global/cfs/cdirs/desi/survey/catalogs/ 

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/main/mkCat_main.py --type QSO --fulld n --survey Y1 --verspec iron --version ${1}pip --clusd y --clusran y --splitGC y --nz y --par y --imsys_colname WEIGHT_RF --basedir /global/cfs/cdirs/desi/survey/catalogs/ --compmd altmtl

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/main/mkCat_main.py --type BGS_BRIGHT --fulld n --survey Y1 --verspec iron --version $1 --clusd y --clusran y --splitGC y --nz y --par y --basedir /global/cfs/cdirs/desi/survey/catalogs/ 

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/main/mkCat_main.py --type BGS_BRIGHT --fulld n --survey Y1 --verspec iron --version ${1}pip --clusd y --clusran y --splitGC y --nz y --par y --basedir /global/cfs/cdirs/desi/survey/catalogs/ --compmd altmtl

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/main/mkCat_main.py --type BGS_ANY --fulld n --survey Y1 --verspec iron --version $1 --clusd y --clusran y --splitGC y --nz y --par y --basedir /global/cfs/cdirs/desi/survey/catalogs/ 

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/main/mkCat_main.py --type BGS_ANY --fulld n --survey Y1 --verspec iron --version ${1}pip --clusd y --clusran y --splitGC y --nz y --par y --basedir /global/cfs/cdirs/desi/survey/catalogs/ --compmd altmtl
