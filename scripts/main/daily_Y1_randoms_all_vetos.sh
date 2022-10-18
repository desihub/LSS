#!/bin/bash
#script to run main/daily randoms through for a given tracer type, in serial, one node for each tracer
#to be run after getting 4 interactive nodes (e.g., via salloc -N 7 -C haswell -t 04:00:00 --qos interactive --account desi)
#1st argument should be tracer type and 2nd should be whether or not to reject qso targets

srun -N 1 python mkCat_main_ran.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --notqso n --survey Y1
srun -N 1 python mkCat_main_ran.py --type ELG_LOP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --notqso y --survey Y1
srun -N 1 python mkCat_main_ran.py --type ELG_LOP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec --survey Y1
srun -N 1 python mkCat_main_ran.py --type ELG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --survey Y1
srun -N 1 python mkCat_main_ran.py --type BGS_ANY --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --survey Y1
srun -N 1 python mkCat_main_ran.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --survey Y1
srun -N 1 python mkCat_main_ran.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --add_veto y --apply_veto y --verspec daily --survey Y1
