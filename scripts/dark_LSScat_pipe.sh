#!/bin/bash

set -e

source /global/common/software/desi/desi_environment.sh main
#export LSSCODE=$HOME ; do this before script, e.g., export LSSCODE=$HOME/LSScode for desica
PYTHONPATH=$PYTHONPATH:$LSSCODE/LSS/py

verspec=loa-v1
survey=DA2

#srun -N 1 -C cpu -t 04:00:00 -q interactive  python $LSSCODE/LSS/scripts/main/mkCat_main.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec $verspec --survey $survey --version $1

#srun -N 1 -C cpu -t 04:00:00 -q interactive  python $LSSCODE/LSS/scripts/main/mkCat_main.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec $verspec --survey $survey --version $1

#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main.py --type ELG_LOP --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec $verspec --survey $survey --version $1

#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main.py --type ELG --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec $verspec --survey $survey --version $1


#only should be done for each spec release 
#~50 minutes
#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_ran.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec $verspec --type dark --combwspec y --fullr y --survey $survey --maxr 18 --version $1 --nproc 9

#run this for new version of LSS catalogs
#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_ran.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec $verspec --type dark --combwspec n --fullr y --survey $survey --maxr 18 --version $1 --nproc 9


#python $LSSCODE/LSS/scripts/main/mkCat_main.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --add_veto y --verspec $verspec --survey $survey --maxr 0 --version $1

#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_ran.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --fillran y --apply_veto y --add_veto y --maxr 18 --version $1

#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_ran.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --apply_veto y --maxr 18 --version $1

#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_ran.py --type ELG_LOP --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --apply_veto y --maxr 18 --version $1

#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_ran.py --type ELG --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --apply_veto y --maxr 18 --version $1


#srun -N 1 -C cpu -t 04:00:00 -q interactive  python $LSSCODE/LSS/scripts/main/mkCat_main.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec $verspec --survey $survey  --mkHPmaps y --version $1

#srun -N 1 -C cpu -t 04:00:00 -q interactive  python $LSSCODE/LSS/scripts/main/mkCat_main.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec $verspec --survey $survey  --mkHPmaps y --version $1

#srun -N 1 -C cpu -t 04:00:00 -q interactive  python $LSSCODE/LSS/scripts/main/mkCat_main.py --type ELG_LOP --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec $verspec --survey $survey  --mkHPmaps y --version $1

#srun -N 1 -C cpu -t 04:00:00 -q interactive  python $LSSCODE/LSS/scripts/main/mkCat_main.py --type ELG --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec $verspec --survey $survey  --mkHPmaps y --version $1

#python $LSSCODE/LSS/scripts/main/mkCat_main.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec $verspec --survey $survey --maxr 0 --version $1

#python $LSSCODE/LSS/scripts/main/mkCat_main.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec $verspec --survey $survey --maxr 0 --version $1

#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main.py --type ELG_LOP --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec $verspec --survey $survey --maxr 0 --version $1

#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main.py --type ELG --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec $verspec --survey $survey --maxr 0 --version $1

#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_ran.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --add_tl y --maxr 18 --version $1

#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_ran.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --add_tl y --maxr 18 --version $1

#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_ran.py --type ELG_LOP --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --add_tl y --maxr 18 --version $1

#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_ran.py --type ELG --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --add_tl y --maxr 18 --version $1

#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_map_veto y --verspec $verspec --survey $survey --version $1

#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_map_veto y --verspec $verspec --survey $survey --version $1

#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main.py --type ELG_LOP --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_map_veto y --verspec $verspec --survey $survey --version $1

#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main.py --type ELG --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_map_veto y --verspec $verspec --survey $survey --version $1

#python $LSSCODE/LSS/scripts/main/mkCat_main.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec $verspec --add_weight_zfail y --survey $survey  --use_map_veto _HPmapcut --version $1

#python $LSSCODE/LSS/scripts/main/mkCat_main.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec $verspec --add_weight_zfail y --survey $survey  --use_map_veto _HPmapcut --version $1

#25 minutes, most of time because of multiple writes
#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main.py --type ELG_LOP --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec $verspec --add_weight_zfail y --survey $survey  --use_map_veto _HPmapcut --version $1

#25 minutes, most of time because of multiple writes
#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main.py --type ELG --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec $verspec --add_weight_zfail y --survey $survey  --use_map_veto _HPmapcut --version $1

source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
#export LSSCODE=$HOME
PYTHONPATH=$PYTHONPATH:$LSSCODE/LSS/py

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSSCODE/LSS/scripts/main/mkCat_main.py --type QSO   --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --imsys y --survey $survey --verspec $verspec --imsys_zbin y --use_map_veto _HPmapcut --version $1

python $LSSCODE/LSS/scripts/main/patch_HPmapcut.py --tracers QSO --survey $survey --verspec $verspec --version $1

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSSCODE/LSS/scripts/main/mkCat_main.py --type QSO  --fulld n --survey $survey --verspec $verspec --clusd y --clusran y --splitGC y --nz y --par y --imsys_colname WEIGHT_IMLIN --basedir /global/cfs/cdirs/desi/survey/catalogs/ --version $1 --extra_clus_dir 'nonKP/'


#narrow dz bins take a while
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --imsys y --survey $survey --verspec $verspec --imsys_zbin y --use_map_veto _HPmapcut --version $1

#~3 minutes
python $LSSCODE/LSS/scripts/main/patch_HPmapcut.py --tracers LRG --survey $survey --verspec $verspec --version $1

#~5 minutes
srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSSCODE/LSS/scripts/main/mkCat_main.py --type LRG --fulld n --survey $survey --verspec $verspec --clusd y --clusran y --splitGC y --nz y --par y --imsys_colname WEIGHT_IMLIN --basedir /global/cfs/cdirs/desi/survey/catalogs/ --version $1 --extra_clus_dir 'nonKP/'

#$LSSCODE/LSS/scripts/Y3_sysnetELG_zbins.sh $1

#srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSSCODE/LSS/scripts/main/mkCat_main.py --type ELG_LOP --notqso y  --fulld n --survey $survey --verspec $verspec --clusd y --clusran y --splitGC y --nz y --par y --imsys_colname WEIGHT_SN --basedir /global/cfs/cdirs/desi/survey/catalogs/ --version $1 --extra_clus_dir 'nonKP/'



#mkdir /global/cfs/cdirs/desi/survey/catalogs/$survey/LSS/$verspec/LSScats/$1/nonKP/
#mv /global/cfs/cdirs/desi/survey/catalogs/$survey/LSS/$verspec/LSScats/$1/*clustering* /global/cfs/cdirs/desi/survey/catalogs/$survey/LSS/$verspec/LSScats/$1/nonKP/
