#!/bin/bash

srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/scripts/main/mkCat_main_ran.py  --type dark1b  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec daily  --fullr y
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/scripts/main/mkCat_main_ran.py  --type dark  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec daily  --fullr y
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/scripts/main/mkCat_main_ran.py  --type bright  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec daily  --fullr y


