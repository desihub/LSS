#!/bin/bash
#script to run main/daily data through veto, do all tracer types
#to be run after getting interaction node (e.g., via salloc -N 1 -C haswell -t 04:00:00 --qos interactive --account desi)

python mkCat_main.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec daily 
python mkCat_main.py --type BGS_ANY --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec daily 
python mkCat_main.py --type ELG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec daily 
python mkCat_main.py --type ELG_LOP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec daily 
python mkCat_main.py --type ELG_LOP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec daily --notqso y
python mkCat_main.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec daily 
python mkCat_main.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec daily 
