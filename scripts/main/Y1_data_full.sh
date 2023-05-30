#!/bin/bash
#script to run main/daily data through up to veto, do all tracer types

python $LSSCODE/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec iron --survey Y1 --version $1
python $LSSCODE/LSS/scripts/main/mkCat_main.py --type ELG_LOP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec iron --notqso y --survey Y1 --version $1
python $LSSCODE/LSS/scripts/main/mkCat_main.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec iron --survey Y1 --version $1
python $LSSCODE/LSS/scripts/main/mkCat_main.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec iron --survey Y1 --version $1
