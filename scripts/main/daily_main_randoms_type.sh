#!/bin/bash
#script to run main/daily randoms through for a given tracer type
#1st argument should be tracer type and 2nd should be whether or not to reject qso targets

python mkCat_main_ran_px.py  --type $1  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec daily --rfa n --combhp n --combfull n --fullr y --notqso $2
python mkCat_main_ran_px.py  --type $1  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec daily --rfa n --combhp n --combfull y --fullr n --notqso $2
python mkCat_main_ran.py --type $1 --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec daily --notqso $2
