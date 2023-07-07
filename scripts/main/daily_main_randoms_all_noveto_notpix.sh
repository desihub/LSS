#!/bin/bash

srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_ran.py  --type QSO  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec daily  --fullr y
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_ran.py  --type LRG  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec daily  --fullr y
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_ran.py  --type ELG  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec daily  --fullr y
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_ran.py  --type ELG_LOP  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec daily  --fullr y
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_ran.py  --type ELG_LOP  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec daily  --fullr y --notqso y
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_ran.py  --type BGS_BRIGHT  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec daily  --fullr y
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_ran.py  --type BGS_ANY  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec daily  --fullr y


