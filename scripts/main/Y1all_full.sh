#!/bin/bash
#make sure to set $LSSCODE to wherever the LSS git repo is (e.g., $HOME)
#provide the version, e.g. Y1all_full.sh 0.1
source /global/common/software/desi/desi_environment.sh main
PYTHONPATH=$PYTHONPATH:$HOME/LSS


#python $LSSCODE/LSS/scripts/main/combdata_main.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec iron --prog dark --survey Y1 --counts_only y --dospec n
#python $LSSCODE/LSS/scripts/main/combdata_main.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec iron --prog bright --survey Y1 --counts_only y --dospec n

#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_Y1ran_px.py  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec iron --type dark 
#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_Y1ran_px.py  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec iron --type bright 

#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_Y1ran_px.py  --type ELG_LOP  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec iron --combhp n --combfull y --fullr y --notqso y --version $1
#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_Y1ran_px.py  --type LRG  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec iron --combhp n --combfull y --fullr y  --version $1
#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_Y1ran_px.py  --type QSO  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec iron --combhp n --combfull y --fullr y  --version $1
#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_Y1ran_px.py  --type BGS_BRIGHT  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec iron --combhp n --combfull y --fullr y  --version $1

#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_ran.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr y --verspec iron --notqso n --survey Y1 --version $1
#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_ran.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --fillran y --apply_veto y --verspec iron --notqso n --survey Y1 --maxr 9 --version $1
#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_ran.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --fillran y --apply_veto y --verspec iron --notqso n --survey Y1 --minr 9 --version $1
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_ran.py --type ELG_LOP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr y --verspec iron --notqso y --survey Y1 --version $1
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_ran.py --type ELG_LOP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --fillran y --apply_veto y --verspec iron --notqso y --survey Y1 --maxr 9 --version $1
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_ran.py --type ELG_LOP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --fillran y --apply_veto y --verspec iron --notqso y --survey Y1 --minr 9 --version $1
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_ran.py --type BGS_BRIGHT  --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr y --verspec iron --notqso n --survey Y1 --version $1
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_ran.py --type BGS_BRIGHT  --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --fillran y --apply_veto y --verspec iron --notqso n --survey Y1 --maxr 9 --version $1
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_ran.py --type BGS_BRIGHT  --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --fillran y --apply_veto y --verspec iron --notqso n --survey Y1 --minr 9 --version $1
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_ran.py --type LRG  --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr y --verspec iron --notqso n --survey Y1 --version $1
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_ran.py --type LRG  --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --fillran y --add_veto y --apply_veto y --verspec iron --notqso n --survey Y1 --maxr 9 --version $1
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main_ran.py --type LRG  --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --fillran y --add_veto y --apply_veto y --verspec iron --notqso n --survey Y1 --minr 9 --version $1


source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main

python $LSSCODE/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y   --apply_veto y --mkHPmaps y --verspec iron --survey Y1 --maxr 0 --regressis y --add_regressis y --add_weight_zfail y --version $1
python $LSSCODE/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n   --verspec iron --survey Y1 --maxr 0 --regressis y --add_regressis y --version $1
python $LSSCODE/LSS/scripts/main/mkCat_main.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --add_veto y   --apply_veto y --mkHPmaps y --verspec iron --survey Y1 --maxr 0 --regressis y --add_regressis y --add_weight_zfail y --version $1
python $LSSCODE/LSS/scripts/main/mkCat_main.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --apply_veto y --mkHPmaps y --verspec iron --survey Y1 --maxr 0 --regressis y --add_regressis y --add_weight_zfail y --version $1
python $LSSCODE/LSS/scripts/main/mkCat_main.py --type ELG_LOP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --apply_veto y --mkHPmaps y --verspec iron --survey Y1 --notqso y --maxr 0 --regressis y --add_regressis y --add_weight_zfail y --version $1

srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSSCODE/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT --verspec iron --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --add_ke y --survey Y1 --version $1

python $LSSCODE/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n   --FKPfull y  --verspec iron --survey Y1 --maxr 0  --version $1
python $LSSCODE/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT-21.5 --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n   --FKPfull y  --verspec iron --survey Y1 --maxr 0  --version $1
python $LSSCODE/LSS/scripts/main/mkCat_main.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --FKPfull y  --verspec iron --survey Y1 --maxr 0 --version $1
python $LSSCODE/LSS/scripts/main/mkCat_main.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --FKPfull y  --verspec iron --survey Y1 --maxr 0  --version $1
python $LSSCODE/LSS/scripts/main/mkCat_main.py --type ELG_LOP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --FKPfull y --verspec iron --survey Y1 --notqso y --maxr 0 --version $1

python $LSSCODE/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n   --nzfull y  --verspec iron --survey Y1 --maxr 0  --version $1
python $LSSCODE/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT-21.5 --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n   --nzfull y  --verspec iron --survey Y1 --maxr 0  --version $1
python $LSSCODE/LSS/scripts/main/mkCat_main.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --nzfull y  --verspec iron --survey Y1 --maxr 0 --version $1
python $LSSCODE/LSS/scripts/main/mkCat_main.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --nzfull y  --verspec iron --survey Y1 --maxr 0  --version $1
python $LSSCODE/LSS/scripts/main/mkCat_main.py --type ELG_LOP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --nzfull y --verspec iron --survey Y1 --notqso y --maxr 0 --version $1


python $LSSCODE/LSS/scripts/plotting/densvs_improp_tpixw.py --version $1
python $LSSCODE/LSS/scripts/validation/validation_tsnr_zbin.py --version $1
python $LSSCODE/LSS/scripts/validation/validation_sky.py --version $1
python $LSSCODE/LSS/scripts/validation/validation_fiber.py --version $1
python $LSSCODE/LSS/scripts/validation/validation_focal.py --version $1
python $LSSCODE/LSS/scripts/validation/validation_cl.py --version $1
