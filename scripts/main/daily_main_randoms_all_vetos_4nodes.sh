#!/bin/bash
#script to run main/daily randoms through for a given tracer type
#to be run after getting 4 interactive nodes (e.g., via salloc -N 4 -C haswell -t 04:00:00 --qos interactive --account desi)
#1st argument should be tracer type and 2nd should be whether or not to reject qso targets

srun -N 1 python mkCat_main_ran.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --notqso n --maxr 5 &
srun -N 1 python mkCat_main_ran.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --notqso n --minr 5 --maxr 10 &
srun -N 1 python mkCat_main_ran.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --notqso n --minr 10 --maxr 14 &
srun -N 1 python mkCat_main_ran.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --notqso n --minr 14 &
wait

srun -N 1 python mkCat_main_ran.py --type ELG_LOP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --notqso y --maxr 5 &
srun -N 1 python mkCat_main_ran.py --type ELG_LOP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --notqso y --minr 5 --maxr 10 &
srun -N 1 python mkCat_main_ran.py --type ELG_LOP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --notqso y --minr 10 --maxr 14 &
srun -N 1 python mkCat_main_ran.py --type ELG_LOP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --notqso y --minr 14 &
wait 

srun -N 1 python mkCat_main_ran.py --type ELG_LOP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --maxr 5 &
srun -N 1 python mkCat_main_ran.py --type ELG_LOP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --minr 5 --maxr 10 &
srun -N 1 python mkCat_main_ran.py --type ELG_LOP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --minr 10 --maxr 14 & 
srun -N 1 python mkCat_main_ran.py --type ELG_LOP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --minr 14 & 
wait

srun -N 1 python mkCat_main_ran.py --type ELG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --maxr 5 &
srun -N 1 python mkCat_main_ran.py --type ELG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --minr 5 --maxr 10 &
srun -N 1 python mkCat_main_ran.py --type ELG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --minr 10 --maxr 14 & 
srun -N 1 python mkCat_main_ran.py --type ELG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --minr 14 & 
wait

srun -N 1 python mkCat_main_ran.py --type BGS_ANY --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --maxr 5 &
srun -N 1 python mkCat_main_ran.py --type BGS_ANY --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --minr 5 --maxr 10 &
srun -N 1 python mkCat_main_ran.py --type BGS_ANY --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --minr 10 --maxr 14 &
srun -N 1 python mkCat_main_ran.py --type BGS_ANY --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --minr 14 &
wait

srun -N 1 python mkCat_main_ran.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --maxr 5 &
srun -N 1 python mkCat_main_ran.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --minr 5 --maxr 10 &
srun -N 1 python mkCat_main_ran.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --minr 10 --maxr 14 &
srun -N 1 python mkCat_main_ran.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --minr 14 &
wait

srun -N 1 python mkCat_main_ran.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --maxr 5 &
srun -N 1 python mkCat_main_ran.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --minr 5 --maxr 10 &
srun -N 1 python mkCat_main_ran.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --minr 10 --maxr 14 &
srun -N 1 python mkCat_main_ran.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --minr 14 &
