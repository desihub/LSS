#!/bin/bash

source /global/common/software/desi/desi_environment.sh master
PYTHONPATH=$PYTHONPATH:$HOME/LSS
PYTHONPATH=$PYTHONPATH:/global/homes/a/ajross/.local/lib/python3.8/site-packages/


python main/mkCat_main.py --type LRG --verspec guadalupe --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --clusd y --clusran y --imsys n --nz y --regressis y --add_regressis y --survey DA02 --version $1 
python main/mkCat_main.py --type ELG_LOP --verspec guadalupe --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --clusd y --clusran y --imsys n --nz y --regressis y --add_regressis y --survey DA02 --notqso y --version $1 
python main/mkCat_main.py --type QSO --verspec guadalupe --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --clusd y --clusran y --imsys n --nz y --regressis y --add_regressis y --survey DA02 --version $1 
python main/mkCat_main.py --type BGS_BRIGHT --verspec guadalupe --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --clusd y --clusran y --imsys n --nz y --regressis y --add_regressis y --survey DA02 --version $1 


