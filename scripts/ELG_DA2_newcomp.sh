#!/bin/bash

set -e

source /global/common/software/desi/desi_environment.sh main
module load LSS/main
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main


verspec=loa-v1
survey=DA2
LSS=$HOME/LSScode
PYTHONPATH=$LSS/LSS/py:$PYTHONPATH

srun -N 1 -C cpu -t 04:00:00 -q interactive  python $LSS/LSS/scripts/main/mkCat_main.py --type ELG_LOP --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec $verspec --survey $survey --version $1 #--mkgtl

srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type ELG_LOP --notqso y  --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --add_veto y --verspec $verspec --survey $survey --maxr 0 --version $1

srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type ELG_LOP --notqso y  --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec $verspec --survey $survey --maxr 0 --version $1

srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --type ELG_LOP --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --add_tl y --maxr 18 --version $1

srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type ELG_LOP --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_map_veto y --verspec $verspec --survey $survey --version $1

srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type ELG_LOP --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec $verspec --add_weight_zfail y --survey $survey  --use_map_veto _HPmapcut --version $1

#10 minutes
srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSS/LSS/scripts/main/mkCat_main.py --type ELG_LOP --notqso y  --fulld n --survey $survey --verspec $verspec --clusd y --clusran y --splitGC y --nz y --par y  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --version $1 --extra_clus_dir 'nonKP/' 


