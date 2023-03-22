#!/bin/bash

srun -N 1 -C cpu -t 04:00:00 -q interactive python scripts/main/mkCat_main_ran.py --type ELG_LOP --apply_veto y --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec iron --fillran y --survey Y1 --version $1 --notqso y
python scripts/main/mkCat_main.py --type ELG_LOP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec iron --survey Y1  --mkHPmaps y --version $1 --notqso y

srun -N 1 -C cpu -t 04:00:00 -q interactive python scripts/main/mkCat_main_ran.py --type BGS_BRIGHT --apply_veto y --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec iron --fillran y --survey Y1 --version $1 
python scripts/main/mkCat_main.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec iron --survey Y1  --mkHPmaps y --version $1 

srun -N 1 -C cpu -t 04:00:00 -q interactive python scripts/main/mkCat_main_ran.py --type LRG --apply_veto y --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec iron --fillran y --survey Y1 --version $1
python scripts/main/mkCat_main.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec iron --survey Y1  --mkHPmaps y --version $1

srun -N 1 -C cpu -t 04:00:00 -q interactive python scripts/main/mkCat_main_ran.py --type QSO --apply_veto y --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec iron --fillran y --survey Y1 --version $1
python scripts/main/mkCat_main.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec iron --survey Y1  --mkHPmaps y --version $1