#!/bin/bash
#script to run main/daily randoms through for a given tracer type
#to be run after getting interactive job with at least two node (e.g., via salloc -N 2 -C haswell -t 04:00:00 --qos interactive --account desi)
#1st argument should be tracer type and 2nd should be whether or not to reject qso targets

srun -N 1 python mkCat_main_ran_px.py  --type $1  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec daily --rfa n --combhp n --combfull n --fullr y --notqso $2 --minr 9 &
srun -N 1 python mkCat_main_ran_px.py  --type $1  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec daily --rfa n --combhp n --combfull n --fullr y --notqso $2 --maxr 9 &
wait
srun -N 1 python  mkCat_main_ran_px.py  --type $1  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec daily --rfa n --combhp n --combfull y --fullr n --notqso $2 --minr 9 &
srun -N 1 python mkCat_main_ran_px.py  --type $1  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec daily --rfa n --combhp n --combfull y --fullr n --notqso $2 --maxr 9 &
