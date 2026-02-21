#!/bin/bash
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main

PYTHONPATH=$PYTHONPATH:$HOME/LSS

python scripts/main/mkCat_main.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --type ELG_LOP --notqso y --regressis y --add_regressis y --fulld n --survey Y1 --verspec iron --version $1

python scripts/main/mkCat_main.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --type LRG --regressis y --add_regressis y --fulld n --survey Y1 --verspec iron --version $1

python scripts/main/mkCat_main.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --type QSO --regressis y --add_regressis y --fulld n --survey Y1 --verspec iron --version $1

python scripts/main/mkCat_main.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --type BGS_BRIGHT-21.5 --regressis y --add_regressis y --fulld n --survey Y1 --verspec iron --version $1
