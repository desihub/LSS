#!/bin/bash

python $LSSCODE/LSS/scripts/main/mkCat_main.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec daily --maxr 0 --survey Y1
python $LSSCODE/LSS/scripts/main/mkCat_main.py --type ELG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec daily --maxr 0 --survey Y1
python $LSSCODE/LSS/scripts/main/mkCat_main.py --type ELG_LOP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec daily --maxr 0 --survey Y1
python $LSSCODE/LSS/scripts/main/mkCat_main.py --type ELG_LOP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec daily --maxr 0 --notqso y --survey Y1
python $LSSCODE/LSS/scripts/main/mkCat_main.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --add_veto y --apply_veto y --verspec daily --maxr 0 --survey Y1
python $LSSCODE/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec daily --maxr 0 --survey Y1
python $LSSCODE/LSS/scripts/main/mkCat_main.py --type BGS_ANY --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec daily --maxr 0 --survey Y1
