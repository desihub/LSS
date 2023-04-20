#!/bin/bash

python scripts/main/mkCat_main.py --type ELG_LOP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec iron --survey Y1  --mkHPmaps y --version $1 --notqso y

python scripts/main/mkCat_main.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec iron --survey Y1  --mkHPmaps y --version $1 

python scripts/main/mkCat_main.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec iron --survey Y1  --mkHPmaps y --version $1

python scripts/main/mkCat_main.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec iron --survey Y1  --mkHPmaps y --version $1