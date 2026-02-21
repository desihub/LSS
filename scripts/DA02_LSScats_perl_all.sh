#!/bin/bash

source /global/common/software/desi/desi_environment.sh master
PYTHONPATH=$PYTHONPATH:$HOME/LSS
PYTHONPATH=$PYTHONPATH:/global/homes/a/ajross/.local/lib/python3.8/site-packages/

python main/mkCat_main.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec guadalupe --survey DA02 --version $1
python main/mkCat_main.py --type ELG_LOP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec guadalupe --survey DA02 --notqso y --version $1
python main/mkCat_main.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec guadalupe --survey DA02 --version $1
python main/mkCat_main.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec guadalupe --survey DA02 --version $1

srun -N 1 -C cpu -t 01:00:00 -q interactive python main/mkCat_main_ran.py  --type LRG  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec guadalupe --fullr y --survey DA02 --version $1
srun -N 1 -C cpu -t 01:00:00 -q interactive python main/mkCat_main_ran.py  --type ELG_LOP  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec guadalupe --fullr y --survey DA02 --notqso y --version $1
srun -N 1 -C cpu -t 01:00:00 -q interactive python main/mkCat_main_ran.py  --type QSO  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec guadalupe --fullr y --survey DA02 --version $1
srun -N 1 -C cpu -t 01:00:00 -q interactive python main/mkCat_main_ran.py  --type BGS_BRIGHT  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec guadalupe --survey DA02 --fullr y --version $1

python main/mkCat_main.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --add_veto y --apply_veto y --verspec guadalupe --survey DA02 --version $1 --maxr 0
python main/mkCat_main.py --type ELG_LOP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --add_veto y --apply_veto y --verspec guadalupe --survey DA02 --version $1 --maxr 0 --notqso y
python main/mkCat_main.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --apply_veto y --verspec guadalupe --survey DA02 --version $1 --maxr 0
python main/mkCat_main.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --apply_veto y --verspec guadalupe --survey DA02 --version $1 --maxr 0

srun -N 1 -C cpu -t 01:00:00 -q interactive python main/mkCat_main_ran.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --add_veto y --apply_veto y --verspec guadalupe --survey DA02 --version $1 
srun -N 1 -C cpu -t 01:00:00 -q interactive python main/mkCat_main_ran.py --type ELG_LOP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --add_veto y --apply_veto y --verspec guadalupe --survey DA02 --notqso y --version $1 
srun -N 1 -C cpu -t 01:00:00 -q interactive python main/mkCat_main_ran.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec guadalupe --survey DA02 --version $1 
srun -N 1 -C cpu -t 01:00:00 -q interactive python main/mkCat_main_ran.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec guadalupe --survey DA02 --version $1 


python main/mkCat_main.py --type LRG --verspec guadalupe --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --clusd y --clusran y --imsys n --nz y --regressis y --add_regressis y --survey DA02 --version $1 
python main/mkCat_main.py --type ELG_LOP --verspec guadalupe --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --clusd y --clusran y --imsys n --nz y --regressis y --add_regressis y --survey DA02 --notqso y --version $1 
python main/mkCat_main.py --type QSO --verspec guadalupe --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --clusd y --clusran y --imsys n --nz y --regressis y --add_regressis y --survey DA02 --version $1 
python main/mkCat_main.py --type BGS_BRIGHT --verspec guadalupe --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --clusd y --clusran y --imsys n --nz y --regressis y --add_regressis y --survey DA02 --version $1 


