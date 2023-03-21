#!/bin/bash

srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_ran_px.py  --type BGS_ANY  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec daily --rfa n --combhp n --combfull y --fullr y --refullr y
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_ran_px.py  --type QSO  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec daily --rfa n --combhp n --combfull y --fullr y --refullr y
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_ran_px.py  --type LRG  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec daily --rfa n --combhp n --combfull y --fullr y --refullr y
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_ran_px.py  --type ELG  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec daily --rfa n --combhp n --combfull y --fullr y --refullr y
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_ran_px.py  --type ELG_LOP  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec daily --rfa n --combhp n --combfull y --fullr y --refullr y
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_ran_px.py  --type ELG_LOP  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec daily --rfa n --combhp n --combfull y --fullr y --notqso y --refullr y
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_ran_px.py  --type BGS_BRIGHT  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec daily --rfa n --combhp n --combfull y --fullr y --refullr y

