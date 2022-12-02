#!/bin/bash
#script to run main/daily data through up to veto, do all tracer types

set -e

python $LSSCODE/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --regressis y --add_regressis y --verspec $1 --survey Y1 --add_weight_zfail y
python $LSSCODE/LSS/scripts/main/mkCat_main.py --type BGS_ANY --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --regressis y --add_regressis y --verspec $1 --survey Y1 --add_weight_zfail y
python $LSSCODE/LSS/scripts/main/mkCat_main.py --type ELG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --regressis y --add_regressis y --verspec $1 --survey Y1 --add_weight_zfail y
python $LSSCODE/LSS/scripts/main/mkCat_main.py --type ELG_LOP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --regressis y --add_regressis y --verspec $1 --survey Y1 --add_weight_zfail y
python $LSSCODE/LSS/scripts/main/mkCat_main.py --type ELG_LOP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --regressis y --add_regressis y --verspec $1 --notqso y --survey Y1 --add_weight_zfail y
python $LSSCODE/LSS/scripts/main/mkCat_main.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --regressis y --add_regressis y --verspec $1 --survey Y1 --add_weight_zfail y
python $LSSCODE/LSS/scripts/main/mkCat_main.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --regressis y --add_regressis y --verspec $1 --survey Y1 --add_weight_zfail y
