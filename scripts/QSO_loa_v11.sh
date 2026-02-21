#!/bin/bash

set -e

source /global/common/software/desi/desi_environment.sh main
export LSSCODE=$HOME
PYTHONPATH=$PYTHONPATH:$LSSCODE/LSS

verspec=loa-v1
echo $verspec
survey=DA2
TRACER='QSO'

#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main.py --type $TRACER  --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec $verspec --survey $survey --version $1

#python $LSSCODE/LSS/scripts/main/mkCat_main.py --type $TRACER  --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec $verspec --survey $survey --maxr 0 --version $1

#python $LSSCODE/LSS/scripts/main/mkCat_main.py --type $TRACER  --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_map_veto y --maxr 0 --verspec $verspec --survey $survey --version $1

source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
export LSSCODE=$HOME
PYTHONPATH=$PYTHONPATH:$LSSCODE/LSS

#python $LSSCODE/LSS/scripts/main/mkCat_main.py --type $TRACER  --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec $verspec --survey $survey --add_weight_zfail y  --use_map_veto _HPmapcut --version $1

python $LSSCODE/LSS/scripts/main/mkCat_main.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --imsys y --survey $survey --verspec $verspec --imsys_zbin y --use_map_veto _HPmapcut --version $1

#~3 minutes
python $LSSCODE/LSS/scripts/main/patch_HPmapcut.py --tracers QSO --survey $survey --verspec $verspec --version $1

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSSCODE/LSS/scripts/main/mkCat_main.py --type QSO --fulld n --survey $survey --verspec $verspec --clusd y --clusran y --splitGC y --nz y --par y --imsys_colname WEIGHT_IMLIN --basedir /global/cfs/cdirs/desi/survey/catalogs/ --version $1 --extra_clus_dir 'nonKP/'
