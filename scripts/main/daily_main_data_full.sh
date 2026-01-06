#!/bin/bash
#script to run main/daily data through up to veto, do all tracer types

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSSCODE/scripts/main/mkCat_main.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec daily 

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSSCODE/scripts/main/mkCat_main.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec daily 
srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSSCODE/scripts/main/mkCat_main.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec daily 

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSSCODE/scripts/main/mkCat_main.py --type BGS_ANY --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec daily 
srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSSCODE/scripts/main/mkCat_main.py --type ELG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec daily 
srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSSCODE/scripts/main/mkCat_main.py --type ELG_LOP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec daily 
srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSSCODE/scripts/main/mkCat_main.py --type ELG_LOP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec daily --notqso y

