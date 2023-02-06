#!/bin/bash
#script to run main/daily data through up to veto, do all tracer types

python $LSSCODE/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec $1 --survey Y1 --maxr 0
srun -N 1 -C cpu -t 04:00:00 -q interactive python mkCat_main_ran.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec $1 --survey Y1
python $LSSCODE/LSS/scripts/main/mkCat_main.py --type BGS_ANY --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec $1 --survey Y1 --maxr 0
srun -N 1 -C cpu -t 04:00:00 -q interactive python mkCat_main_ran.py --type BGS_ANY --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec $1 --survey Y1
python $LSSCODE/LSS/scripts/main/mkCat_main.py --type ELG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec $1 --survey Y1 --maxr 0
srun -N 1 -C cpu -t 04:00:00 -q interactive python mkCat_main_ran.py --type ELG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec $1 --survey Y1
python $LSSCODE/LSS/scripts/main/mkCat_main.py --type ELG_LOP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec $1 --survey Y1 --maxr 0
srun -N 1 -C cpu -t 04:00:00 -q interactive python mkCat_main_ran.py --type ELG_LOP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec $1 --survey Y1
python $LSSCODE/LSS/scripts/main/mkCat_main.py --type ELG_LOP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec $1 --survey Y1 --maxr 0 --notqso y
srun -N 1 -C cpu -t 04:00:00 -q interactive python mkCat_main_ran.py --type ELG_LOP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec $1 --survey Y1 --notqso y
python $LSSCODE/LSS/scripts/main/mkCat_main.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec $1 --survey Y1 --maxr 0
srun -N 1 -C cpu -t 04:00:00 -q interactive python mkCat_main_ran.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec $1 --survey Y1
python $LSSCODE/LSS/scripts/main/mkCat_main.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec $1 --survey Y1 --maxr 0
srun -N 1 -C cpu -t 04:00:00 -q interactive python mkCat_main_ran.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec $1 --survey Y1
