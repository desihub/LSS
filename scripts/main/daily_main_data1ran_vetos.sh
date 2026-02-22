#!/bin/bash
#to run as individual jobs `srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi <command>`

python $LSSCODE/scripts/main/mkCat_main.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec daily --maxr 1 
srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSSCODE/scripts/main/mkCat_main.py --type ELG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec daily --maxr 1 
srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSSCODE/scripts/main/mkCat_main.py --type ELG_LOP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec daily --maxr 1 
srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSSCODE/scripts/main/mkCat_main.py --type ELG_LOP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec daily --maxr 1 --notqso y
python $LSSCODE/scripts/main/mkCat_main.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec daily --maxr 1 
python $LSSCODE/scripts/main/mkCat_main.py --type LGE --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec daily --maxr 1 
python $LSSCODE/scripts/main/mkCat_main.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec daily --maxr 1 
srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSSCODE/scripts/main/mkCat_main.py --type BGS_ANY --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec daily --maxr 1 
