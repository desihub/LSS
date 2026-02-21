#!/bin/bash

set -e

source /global/common/software/desi/desi_environment.sh main
export LSSCODE=$HOME
PYTHONPATH=$PYTHONPATH:$LSSCODE/LSS

TRACER='QSO'
verspec=loa-v1
survey=DA2


srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main.py --type $TRACER  --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec $verspec --survey $survey --version $1

#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_ran.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec iron --type $TRACER  --combwspec n --fullr y --survey Y1 --maxr 18 --version $1

#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_ran.py --type $TRACER  --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec iron --survey Y1 --fillran y --apply_veto y --maxr 18 --version $1

#python $LSSCODE/LSS/scripts/main/mkCat_main.py --type $TRACER  --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec $verspec --survey $survey  --mkHPmaps y --version $1

python $LSSCODE/LSS/scripts/main/mkCat_main.py --type $TRACER  --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec $verspec --survey $survey --maxr 0 --version $1

#python $LSSCODE/LSS/scripts/main/mkCat_main_ran.py --type $TRACER  --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --add_tl y --maxr 18 --par n --version $1

srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main.py --type $TRACER  --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_map_veto y --verspec $verspec --survey $survey --version $1

python $LSSCODE/LSS/scripts/main/mkCat_main.py --type $TRACER  --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec $verspec --survey $survey --add_weight_zfail y  --use_map_veto _HPmapcut --version $1

#python $LSSCODE/LSS/scripts/main/mkCat_main.py --type $TRACER   --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --swap20211212 y --verspec $verspec --survey $survey --version $1

python $LSSCODE/LSS/scripts/main/mkCat_main.py --type $TRACER  --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --imsys y  --verspec $verspec --survey $survey --imsys_zbin y --use_map_veto _HPmapcut --version $1

source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
export LSSCODE=$HOME
PYTHONPATH=$PYTHONPATH:$LSSCODE/LSS

python scripts/main/mkCat_main.py --type $TRACER --version $1  --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --regressis y --add_regressis y  --verspec $verspec --survey $survey --imsys_zbin y

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/main/mkCat_main.py --type $TRACER  --fulld n  --verspec $verspec --survey $survey --version $1 --clusd y --clusran y --splitGC y --nz y --par y --imsys_colname WEIGHT_RF --basedir /global/cfs/cdirs/desi/survey/catalogs/ --extra_clus_dir 'nonKP/'

