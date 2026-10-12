#!/bin/bash

#pipeline for just BGS_BRIGHT, intended for DA2 v2.3

set -e

source /global/common/software/desi/desi_environment.sh main
#module load LSS/main
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
#export LSSCODE=$HOME ; do this before script, e.g., export LSSCODE=$HOME/LSScode for desica
#PYTHONPATH=$PYTHONPATH:$LSS/LSS/py

verspec=loa-v1
survey=DA2
LSS=$HOME/LSScode
PYTHONPATH=$LSS/LSS/py:$PYTHONPATH

#full_noveto; add --mkgtl for first dark/bright/dark1b instance, total 1hr 16 minutes (two tracers could be run at once?)
##10 minutes
#srun -N 1 -C cpu -t 04:00:00 -q interactive  python $LSS/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec $verspec --survey $survey --version $1 --mkgtl 



#randoms with spec info
#only should be done for each spec release 
#~1hr 12 minutes
#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec $verspec --type dark --combwspec y --fullr n --survey $survey --maxr 18 --version $1 --nproc 9
#x minutes
#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec $verspec --type dark1b --combwspec y --fullr n --survey $survey --maxr 18 --version $1 --nproc 9

#x minutes
#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec $verspec --type bright --combwspec y --fullr n --survey $survey --maxr 18 --version $1 --nproc 9


#full randoms for each program (dark,bright,dark1b) run this for new version of LSS catalogs 63 minutes total
#24 minutes
#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec $verspec --type bright --combwspec n --fullr y --survey $survey --maxr 18 --version $1 --nproc 9

#fill randoms with properties and apply vetos 128 minutes total
#24 minutes
#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --fillran y --apply_veto y  --maxr 18 --version $1


#make healpix maps 199 minutes total
#47 minutes
#srun -N 1 -C cpu -t 04:00:00 -q interactive  python $LSS/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec $verspec --survey $survey  --mkHPmaps y --version $1

#apply vetos to data, get FRAC_TLOBS info 42 minutes total, one big hang (tracers could be run in parallel?)
#5 minutes
#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec $verspec --survey $survey --maxr 0 --version $1

#Add FRAC_TLOBS info to randoms 70 minutes total
#11 minutes
#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --add_tl y --maxr 18 --version $1


#Apply map veto 153 minutes total (some weird hangs)
#21 minutes
#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_map_veto y --verspec $verspec --survey $survey --version $1

#redshift failure modeling; 106 minutes total
#13 minutes, most of time because of multiple writes
#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec $verspec --add_weight_zfail y --survey $survey  --use_map_veto _HPmapcut --version $1

#clustering catalogs; 85 minutes total; was reading not using dvs_ro, so that likely made it much longer than need be
#13 minutes
#srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSS/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT  --fulld n --survey $survey --verspec $verspec --clusd y --clusran y --splitGC y --nz y --par y  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --version $1 --extra_clus_dir 'nonKP/' --redo_fracz y  #--imsys_colname WEIGHT_IMLIN
#Make -21.35 subsample
srun -N 1 -C cpu -t 01:00:00 --qos interactive --account desi python scripts/mkCat_subsamp.py --input_tracer BGS_BRIGHT --ccut -21.35 --version $1 --survey $survey --verspec $verspec --outdir /global/cfs/cdirs/desi/survey/catalogs/$survey/LSS/$verspec/LSScats/$1/nonKP/  --mkfulldat y --clusd y --clusran y --splitGC y --nz y --imsys_clus y --imsys_clus_ran y --nran4imsys 10 
#Make it with zcmb correction
srun -N 1 -C cpu -t 01:00:00 --qos interactive --account desi python scripts/mkCat_subsamp.py --input_tracer BGS_BRIGHT --ccut zcmb-21.35 --version $1 --survey $survey --verspec $verspec --outdir /global/cfs/cdirs/desi/survey/catalogs/$survey/LSS/$verspec/LSScats/$1/nonKP/  --mkfulldat y --clusd y --clusran y --splitGC y --nz y --imsys_clus y --imsys_clus_ran y --nran4imsys 10 