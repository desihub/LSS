#!/bin/bash

#pipeline for all tracers starting from comb files; 16.3 hours total

set -e

source /global/common/software/desi/desi_environment.sh main
#module load LSS/main
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
#export LSSCODE=$HOME ; do this before script, e.g., export LSSCODE=$HOME/LSScode for desica
#PYTHONPATH=$PYTHONPATH:$LSS/LSS/py

verspec=nevis-v2test
survey=DA3
LSS=$HOME/LSScode
PYTHONPATH=$LSS/LSS/py:$PYTHONPATH

#full_noveto; add --mkgtl for first dark/bright/dark1b instance, total 1hr 16 minutes (two tracers could be run at once?)
#11 minutes
#srun -N 1 -C cpu -t 04:00:00 -q interactive  python $LSS/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec $verspec --survey $survey --version $1 --mkgtl 
#6.5 minutes
#srun -N 1 -C cpu -t 04:00:00 -q interactive  python $LSS/LSS/scripts/main/mkCat_main.py --type BGS_FAINT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec $verspec --survey $survey --version $1 



#randoms with spec info
#only should be done for each spec release 
#~1hr 12 minutes
#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec $verspec --type bright --combwspec y --fullr n --survey $survey --maxr 18 --version $1 --nproc 9


#full randoms for each program (dark,bright,dark1b) run this for new version of LSS catalogs 63 minutes total
#29 minutes
#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec $verspec --type bright --combwspec n --fullr y --survey $survey --maxr 18 --version $1 --nproc 9


#fill randoms with properties and apply vetos 128 minutes total
#26 minutes
#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --fillran y --apply_veto y  --maxr 18 --version $1
#13 minutes
#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --type BGS_FAINT --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --apply_veto y --maxr 18 --version $1
#11 minutes
#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --type ELG_LOP --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --apply_veto y --maxr 18 --version $1
#11 minutes
#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --type ELG_VLO --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --apply_veto y --maxr 18 --version $1


#make healpix maps 199 minutes total
#44 minutes
#srun -N 1 -C cpu -t 04:00:00 -q interactive  python $LSS/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec $verspec --survey $survey  --mkHPmaps y --version $1
#46 minutes
#srun -N 1 -C cpu -t 04:00:00 -q interactive  python $LSS/LSS/scripts/main/mkCat_main.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec $verspec --survey $survey  --mkHPmaps y --version $1
#only one of the ELGs and BGS should be necessary
#44 minutes
#srun -N 1 -C cpu -t 04:00:00 -q interactive  python $LSS/LSS/scripts/main/mkCat_main.py --type ELG_LOP --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec $verspec --survey $survey  --mkHPmaps y --version $1

#apply vetos to data, get FRAC_TLOBS info 42 minutes total, one big hang (tracers could be run in parallel?)
#5 minutes
#python $LSS/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec $verspec --survey $survey --maxr 0 --version $1
#3 minutes
#python $LSS/LSS/scripts/main/mkCat_main.py --type BGS_FAINT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec $verspec --survey $survey --maxr 0 --version $1
#22 minutes, big hang on initial read
#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type ELG_LOP --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec $verspec --survey $survey --maxr 0 --version $1
#3 minutes
#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type ELG_VLO --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec $verspec --survey $survey --maxr 0 --version $1

#Add FRAC_TLOBS info to randoms 70 minutes total
#11 minutes
#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --add_tl y --maxr 18 --version $1
#10 minutes
#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --type BGS_FAINT --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --add_tl y --maxr 18 --version $1
#11 minutes
#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --type ELG_LOP --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --add_tl y --maxr 18 --version $1
#11 minutes
#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --type ELG_VLO --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --add_tl y --maxr 18 --version $1


#Apply map veto 153 minutes total (some weird hangs)
#17 minutes
#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_map_veto y --verspec $verspec --survey $survey --version $1
#15 minutes
#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type BGS_FAINT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_map_veto y --verspec $verspec --survey $survey --version $1
#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_map_veto y --verspec $verspec --survey $survey --version $1 --maxr 0
#52 minutes ... 16 minutes to read data file; something seems wrong
#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type ELG_LOP --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_map_veto y --verspec $verspec --survey $survey --version $1
#15 minutes
#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type ELG_VLO --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_map_veto y --verspec $verspec --survey $survey --version $1

#redshift failure modeling; 106 minutes total
#32 minutes
#srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSS/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec $verspec --add_weight_zfail y --survey $survey  --use_map_veto _HPmapcut --version $1
#16 minutes
#srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSS/LSS/scripts/main/mkCat_main.py --type BGS_FAINT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec $verspec --add_weight_zfail y --survey $survey  --use_map_veto _HPmapcut --version $1

outdir=/global/cfs/cdirs/desi/survey/catalogs/$survey/LSS/$verspec/LSScats/$1/nonKP/
#srun -N 1 -C cpu -t 01:00:00 --qos interactive --account desi python $LSS/LSS/scripts/mkCat_subsamp.py --input_tracer BGS_BRIGHT --ccut -21.35 --version $1 --survey $survey --verspec $verspec --outdir $outdir --compmd NN --mkfulldat y --clusd y --clusran y --splitGC y --nz y --imsys_clus y --imsys_clus_ran y --nran4imsys 10 --imsys_zbin y

#srun -N 1 -C cpu -t 01:00:00 --qos interactive --account desi python $LSS/LSS/scripts/mkCat_subsamp.py --input_tracer BGS_FAINT --ccut -21.35 --version $1 --survey $survey --verspec $verspec --outdir $outdir --compmd NN --mkfulldat y --clusd y --clusran y --splitGC y --nz y --imsys_clus y --imsys_clus_ran y --nran4imsys 10 --imsys_zbin hi

srun -N 1 -C cpu -t 01:00:00 --qos interactive --account desi python $LSS/LSS/scripts/mkCat_subsamp.py --input_tracer BGS_BRIGHT --ccut -20.7 --version $1 --survey $survey --verspec $verspec --outdir $outdir --compmd NN --mkfulldat y --clusd y --clusran y --splitGC y --nz y --imsys_clus y --imsys_clus_ran y --nran4imsys 10 --imsys_zbin y

#combine BGS -21.35
srun -N 1 -C cpu -t 01:00:00 --qos interactive --account desi python $LSS/LSS/scripts/combine_BGSorELG.py --tracer BGS --ccut -21.35 --version $1 --survey $survey --verspec $verspec --outdir $outdir



#clustering catalogs; 85 minutes total; was reading not using dvs_ro, so that likely made it much longer than need be
#11 minutes
#srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSS/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT  --fulld n --survey $survey --verspec $verspec --clusd y --clusran y --splitGC y --nz y --par y  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --version $1 --extra_clus_dir 'nonKP/' --redo_fracz y --nearestneighbor y --des_resamp #--imsys_colname WEIGHT_IMLIN
#10 minutes
#srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSS/LSS/scripts/main/mkCat_main.py --type BGS_FAINT  --fulld n --survey $survey --verspec $verspec --clusd y --clusran y --splitGC y --nz y --par y  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --version $1 --extra_clus_dir 'nonKP/' --redo_fracz y --nearestneighbor y --des_resamp #--imsys_colname WEIGHT_IMLIN
#15 minutes
#srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSS/LSS/scripts/main/mkCat_main.py --type ELG_LOP --notqso y  --fulld n --survey $survey --verspec $verspec --clusd y --clusran y --splitGC y --nz y --par y  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --version $1 --extra_clus_dir 'nonKP/' --redo_fracz y --nearestneighbor y --des_resamp #--imsys_colname WEIGHT_IMLIN
#19 minutes (why so long?)
#srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSS/LSS/scripts/main/mkCat_main.py --type ELG_VLO --notqso y  --fulld n --survey $survey --verspec $verspec --clusd y --clusran y --splitGC y --nz y --par y  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --version $1 --extra_clus_dir 'nonKP/' --redo_fracz y --nearestneighbor y --des_resamp #--imsys_colname WEIGHT_IMLIN

#imaging systematics for non-ELG samples; 47 minutes total
#13 minutes
#srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSS/LSS/scripts/addsys2clus.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/ --version $1 --survey DA3 --extra_clus_dir nonKP --verspec $verspec --doimlin y --imsys_zbin fine --replace_syscol --splitDES
#10 minutes
#srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSS/LSS/scripts/addsys2clus.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/ --version $1 --survey DA3 --extra_clus_dir nonKP --verspec $verspec --doimlin y --imsys_zbin split --replace_syscol --splitDES
