#!/bin/bash

set -e

source /global/common/software/desi/desi_environment.sh main
export LSSCODE=$HOME
PYTHONPATH=$PYTHONPATH:$LSSCODE/LSS

verspec=loa-v1
survey=DA2

#make BGS -21.5 sample
python $LSSCODE/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT-21.5 --verspec $verspec --absmagmd 'nok' --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --survey $survey --version $1

python $LSSCODE/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT-21.5 --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --imsys y --survey $survey --verspec $verspec --imsys_zbin y --use_map_veto _HPmapcut --version $1

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSSCODE/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT-21.5 --fulld n --survey $survey --verspec $verspec --clusd y --clusran y --splitGC y --nz y --par y  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --version $1 --imsys_colname WEIGHT_IMLIN --extra_clus_dir 'nonKP/'

#make BGS -21.35 sample
python $LSSCODE/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT-21.35 --verspec $verspec --absmagmd 'nok' --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --survey $survey --version $1

python $LSSCODE/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT-21.35 --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --imsys y --survey $survey --verspec $verspec --imsys_zbin y --use_map_veto _HPmapcut --version $1

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSSCODE/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT-21.35 --fulld n --survey $survey --verspec  $verspec --clusd y --clusran y --splitGC y --nz y --par y  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --version $1 --imsys_colname WEIGHT_IMLIN --extra_clus_dir 'nonKP/'


#make BGS -20.2 sample
python $LSSCODE/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT-20.2 --verspec $verspec --absmagmd 'nok' --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --survey $survey --version $1

python $LSSCODE/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT-20.2 --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --imsys y --survey $survey --verspec $verspec --imsys_zbin y --use_map_veto _HPmapcut --version $1

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSSCODE/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT-20.2 --fulld n --survey $survey --verspec $verspec --clusd y --clusran y --splitGC y --nz y --par y  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --version $1 --imsys_colname WEIGHT_IMLIN --extra_clus_dir 'nonKP/' 


srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSSCODE/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT --fulld n --survey $survey --verspec $verspec --clusd y --clusran y --splitGC y --nz y --par y  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --version $1 --extra_clus_dir 'nonKP/'

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSSCODE/LSS/scripts/main/mkCat_main.py --type BGS_ANY --fulld n --survey $survey --verspec $verspec --clusd y --clusran y --splitGC y --nz y --par y  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --version $1 --extra_clus_dir 'nonKP/'

#mv /global/cfs/cdirs/desi/survey/catalogs/$survey/LSS/$verspec/LSScats/$1/*clustering* /global/cfs/cdirs/desi/survey/catalogs/$survey/LSS/$verspec/LSScats/$1/nonKP/