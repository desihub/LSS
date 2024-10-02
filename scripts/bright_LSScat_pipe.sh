#!/bin/bash

set -e

source /global/common/software/desi/desi_environment.sh main
export LSSCODE=$HOME
PYTHONPATH=$PYTHONPATH:$LSSCODE/LSS/py

verspec=kibo-v1
survey=DA2

#srun -N 1 -C cpu -t 04:00:00 -q interactive  python $LSSCODE/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec $verspec --survey $survey --version $1

#srun -N 1 -C cpu -t 04:00:00 -q interactive  python $LSSCODE/LSS/scripts/main/mkCat_main.py --type BGS_ANY --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec $verspec --survey $survey --version $1

#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_ran.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec $verspec --type bright --combwspec n --fullr y --survey $survey --maxr 18 --version $1 --nproc 9

#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_ran.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --fillran y --maxr 18 --version $1

#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_ran.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --apply_veto y --maxr 18 --version $1

#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_ran.py --type BGS_ANY --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --apply_veto y --maxr 18 --version $1

#python $LSSCODE/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec $verspec --survey $survey  --mkHPmaps y --version $1

#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec $verspec --survey $survey --maxr 0 --version $1

#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main.py --type BGS_ANY --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec $verspec --survey $survey --maxr 0 --version $1

#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_ran.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --add_tl y --maxr 18 --version $1

#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_ran.py --type BGS_ANY --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --add_tl y --maxr 18 --version $1

#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_map_veto y --verspec $verspec --survey $survey --version $1

#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main.py --type BGS_ANY --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_map_veto y --verspec $verspec --survey $survey --version $1

#python $LSSCODE/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec $verspec --add_weight_zfail y --survey $survey  --use_map_veto _HPmapcut --version $1

srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main.py --type BGS_ANY --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec $verspec --add_weight_zfail y --survey $survey  --use_map_veto _HPmapcut --version $1

#make BGS -21.5 sample
python scripts/main/mkCat_main.py --type BGS_BRIGHT-21.5 --verspec $verspec --absmagmd 'nok' --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --survey $survey --version $1


python $LSSCODE/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT-21.5 --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --imsys y --survey $survey --verspec $verspec --imsys_zbin y --use_map_veto _HPmapcut --version $1

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSSCODE/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT --fulld n --survey $survey --verspec $verspec --clusd y --clusran y --splitGC y --nz y --par y  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --version $1 

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSSCODE/LSS/scripts/main/mkCat_main.py --type BGS_ANY --fulld n --survey $survey --verspec $verspec --clusd y --clusran y --splitGC y --nz y --par y  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --version $1 

mv /global/cfs/cdirs/desi/survey/catalogs/$survey/LSS/$verspec/LSScats/$1/*clustering* /global/cfs/cdirs/desi/survey/catalogs/$survey/LSS/$verspec/LSScats/$1/nonKP/