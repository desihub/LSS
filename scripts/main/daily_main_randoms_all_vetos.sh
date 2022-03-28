#!/bin/bash
#script to run main/daily randoms through for a given tracer type
#to be run after getting interaction node (e.g., via salloc -N 1 -C haswell -t 04:00:00 --qos interactive --account desi)
#1st argument should be tracer type and 2nd should be whether or not to reject qso targets

python mkCat_main_ran.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --notqso n --maxr 6
python mkCat_main_ran.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --notqso n --minr 6 --maxr 12
python mkCat_main_ran.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --notqso n --minr 12 

python mkCat_main_ran.py --type ELG_LOP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --notqso y --maxr 6
python mkCat_main_ran.py --type ELG_LOP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --notqso y --minr 6 --maxr 12
python mkCat_main_ran.py --type ELG_LOP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --notqso y --minr 12 

python mkCat_main_ran.py --type ELG_LOP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --maxr 6
python mkCat_main_ran.py --type ELG_LOP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --minr 6 --maxr 12
python mkCat_main_ran.py --type ELG_LOP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --minr 12 

python mkCat_main_ran.py --type ELG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --maxr 6
python mkCat_main_ran.py --type ELG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --minr 6 --maxr 12
python mkCat_main_ran.py --type ELG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --minr 12 

python mkCat_main_ran.py --type BGS_ANY --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --maxr 6
python mkCat_main_ran.py --type BGS_ANY --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --minr 6 --maxr 12
python mkCat_main_ran.py --type BGS_ANY --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --minr 12 

python mkCat_main_ran.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --maxr 6
python mkCat_main_ran.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --minr 6 --maxr 12
python mkCat_main_ran.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --minr 12 

python mkCat_main_ran.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --maxr 6
python mkCat_main_ran.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --minr 6 --maxr 12
python mkCat_main_ran.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --minr 12 
