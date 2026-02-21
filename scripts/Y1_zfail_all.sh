#!/bin/bash
source /global/common/software/desi/desi_environment.sh main

PYTHONPATH=$PYTHONPATH:$HOME/LSS

python scripts/main/mkCat_main.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --type ELG_LOP --notqso y --add_weight_zfail y --fulld n --survey Y1 --verspec iron --version $1

python scripts/main/mkCat_main.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --type LRG --add_weight_zfail y --fulld n --survey Y1 --verspec iron --version $1

python scripts/main/mkCat_main.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --type QSO --add_weight_zfail y --fulld n --survey Y1 --verspec iron --version $1

python scripts/main/mkCat_main.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --type BGS_BRIGHT --add_weight_zfail y --fulld n --survey Y1 --verspec iron --version $1
