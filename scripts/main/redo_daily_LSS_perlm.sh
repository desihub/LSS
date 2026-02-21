#!/bin/bash

set -e

PYTHONPATH=$PYTHONPATH:$LSSCODE/LSS

python $LSSCODE/LSS/scripts/main/combdata_main.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec daily --prog dark --redotarspec y
python $LSSCODE/LSS/scripts/main/combdata_main.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec daily --prog bright --redotarspec y

$LSSCODE/LSS/scripts/main/daily_main_data_full.sh

$LSSCODE/LSS/scripts/main/daily_main_randoms_all_noveto_redo.sh

python $LSSCODE/LSS/scripts/main/mkCat_main.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --add_veto y --verspec daily --maxr 1

$LSSCODE/LSS/scripts/main/daily_main_data1ran_vetos.sh > compnum_afterveto.out

python $LSSCODE/LSS/scripts/main/write_daily_stats.py