#!/bin/bash

set -e

PYTHONPATH=$PYTHONPATH:$LSSCODE/LSS
#module swap fiberassign/5.0.0

python $LSSCODE/LSS/scripts/main/combdata_main.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec daily --prog dark 
python $LSSCODE/LSS/scripts/main/combdata_main.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec daily --prog bright 

python $LSSCODE/LSS/scripts/main/combdata_main.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec daily --prog dark --doqso y --dospec n --combpix n > QSO.out

srun -N 1 -C cpu -t 02:00:00 -q interactive python $LSSCODE/LSS/scripts/mkemlin.py
python $LSSCODE/LSS/scripts/main/combdata_main.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec daily --prog dark --mkemlin y --dospec n --combpix n

$LSSCODE/LSS/scripts/main/daily_main_data_full.sh

srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_ran.py  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec daily --type dark --ranmtl y 
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/getpota_daily_ran.py --prog DARK
python $LSSCODE/LSS/scripts/main/mkCat_main_ran.py  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec daily --type dark --combwspec y --counttiles y --maxr 1

srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_ran.py  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec daily --type bright --ranmtl y 
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/getpota_daily_ran.py --prog BRIGHT
python $LSSCODE/LSS/scripts/main/mkCat_main_ran.py  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec daily --type bright --combwspec y --counttiles y  --maxr 1

$LSSCODE/LSS/scripts/main/daily_main_randoms_all_noveto_notpix0.sh

python $LSSCODE/LSS/scripts/main/mkCat_main.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --add_veto y --verspec daily --maxr 1

$LSSCODE/LSS/scripts/main/daily_main_data1ran_vetos.sh > compnum_afterveto.out

python $LSSCODE/LSS/scripts/main/write_daily_stats.py