#!/bin/bash

#pipeline for all tracers starting from comb files; 16.3 hours total

set -e

source /global/common/software/desi/desi_environment.sh main
#module load LSS/main
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
#export LSSCODE=$HOME ; do this before script, e.g., export LSSCODE=$HOME/LSScode for desica
#PYTHONPATH=$PYTHONPATH:$LSS/LSS/py

verspec=matterhorn-v2
survey=DA3
LSS=$HOME/LSScode
PYTHONPATH=$LSS/LSS/py:$PYTHONPATH

#full_noveto; add --mkgtl for first dark/bright/dark1b instance, total 1hr 16 minutes (two tracers could be run at once?)
#11 minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive  python $LSS/LSS/scripts/main/mkCat_main.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec $verspec --survey $survey --version $1 --mkgtl 
#2 minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive  python $LSS/LSS/scripts/main/mkCat_main.py --type LGE --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec $verspec --survey $survey --version $1 --mkgtl 
#13 minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive  python $LSS/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec $verspec --survey $survey --version $1 --mkgtl 
#6.5 minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive  python $LSS/LSS/scripts/main/mkCat_main.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec $verspec --survey $survey --version $1 
#25 minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type ELG_LOP --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec $verspec --survey $survey --version $1 
#10 minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type ELG_VLO --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec $verspec --survey $survey --version $1 
#7.5 minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type BGS_FAINT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec $verspec --survey $survey --version $1 



#randoms with spec info
#only should be done for each spec release 
#~1hr 12 minutes
#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec $verspec --type dark --combwspec y --fullr n --survey $survey --maxr 18 --version $1 --nproc 9
#x minutes
#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec $verspec --type dark1b --combwspec y --fullr n --survey $survey --maxr 18 --version $1 --nproc 9

#x minutes
#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec $verspec --type bright --combwspec y --fullr n --survey $survey --maxr 18 --version $1 --nproc 9


#full randoms for each program (dark,bright,dark1b) run this for new version of LSS catalogs 63 minutes total
#29 minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec $verspec --type dark --combwspec n --fullr y --survey $survey --maxr 18 --version $1 --nproc 9
#10 minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec $verspec --type dark1b --combwspec n --fullr y --survey $survey --maxr 18 --version $1 --nproc 9
#24 minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec $verspec --type bright --combwspec n --fullr y --survey $survey --maxr 18 --version $1 --nproc 9

#add LRG veto column 11 minutes total
#8 minutes 
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --add_veto y --verspec $verspec --survey $survey --maxr 18 --version $1
#3 minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type LGE --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --add_veto y --verspec $verspec --survey $survey --maxr 18 --version $1 --par y

#fill randoms with properties and apply vetos 128 minutes total
#26 minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --fillran y --apply_veto y  --maxr 18 --version $1
#15 minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --type LGE --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --fillran y --apply_veto y  --maxr 18 --version $1
#13 minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --apply_veto y --maxr 18 --version $1
#11 minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --type ELG_LOP --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --apply_veto y --maxr 18 --version $1
#11 minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --type ELG_VLO --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --apply_veto y --maxr 18 --version $1
#24 minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --fillran y --apply_veto y  --maxr 18 --version $1
#28 minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --type BGS_FAINT --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey  --apply_veto y  --maxr 18 --version $1


#make healpix maps 199 minutes total
#44 minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive  python $LSS/LSS/scripts/main/mkCat_main.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec $verspec --survey $survey  --mkHPmaps y --version $1
#18 minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive  python $LSS/LSS/scripts/main/mkCat_main.py --type LGE --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec $verspec --survey $survey  --mkHPmaps y --version $1
#46 minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive  python $LSS/LSS/scripts/main/mkCat_main.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec $verspec --survey $survey  --mkHPmaps y --version $1
#only one of the ELGs and BGS should be necessary
#44 minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive  python $LSS/LSS/scripts/main/mkCat_main.py --type ELG_LOP --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec $verspec --survey $survey  --mkHPmaps y --version $1
#47 minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive  python $LSS/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec $verspec --survey $survey  --mkHPmaps y --version $1

#apply vetos to data, get FRAC_TLOBS info 42 minutes total, one big hang (tracers could be run in parallel?)
#5 minutes
python $LSS/LSS/scripts/main/mkCat_main.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec $verspec --survey $survey --maxr 0 --version $1
#1 minutes
python $LSS/LSS/scripts/main/mkCat_main.py --type LGE --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec $verspec --survey $survey --maxr 0 --version $1
#3 minutes
python $LSS/LSS/scripts/main/mkCat_main.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec $verspec --survey $survey --maxr 0 --version $1
#22 minutes, big hang on initial read
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type ELG_LOP --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec $verspec --survey $survey --maxr 0 --version $1
#3 minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type ELG_VLO --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec $verspec --survey $survey --maxr 0 --version $1
#5 minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec $verspec --survey $survey --maxr 0 --version $1
#3 minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type BGS_FAINT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec $verspec --survey $survey --maxr 0 --version $1

#Add FRAC_TLOBS info to randoms 70 minutes total
#11 minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --add_tl y --maxr 18 --version $1
#10 minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --add_tl y --maxr 18 --version $1
#11 minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --add_tl y --maxr 18 --version $1
#11 minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --type BGS_FAINT --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --add_tl y --maxr 18 --version $1
#4 minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --type LGE --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --add_tl y --maxr 18 --version $1
#11 minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --type ELG_LOP --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --add_tl y --maxr 18 --version $1
#11 minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --type ELG_VLO --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --add_tl y --maxr 18 --version $1


#Apply map veto 153 minutes total (some weird hangs)
#17 minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_map_veto y --verspec $verspec --survey $survey --version $1
#21 minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_map_veto y --verspec $verspec --survey $survey --version $1
#17 minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type BGS_FAINT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_map_veto y --verspec $verspec --survey $survey --version $1
#16 minutes; seems way too long; 4 minute hang between file output and function exit...
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type LGE --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_map_veto y --verspec $verspec --survey $survey --version $1
#15 minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_map_veto y --verspec $verspec --survey $survey --version $1
#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_map_veto y --verspec $verspec --survey $survey --version $1 --maxr 0
#52 minutes ... 16 minutes to read data file; something seems wrong
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type ELG_LOP --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_map_veto y --verspec $verspec --survey $survey --version $1
#15 minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type ELG_VLO --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_map_veto y --verspec $verspec --survey $survey --version $1

#redshift failure modeling; 106 minutes total
#32 minutes
python $LSS/LSS/scripts/main/mkCat_main.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec $verspec --add_weight_zfail y --survey $survey  --use_map_veto _HPmapcut --version $1
#3 minutes
python $LSS/LSS/scripts/main/mkCat_main.py --type LGE --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec $verspec --add_weight_zfail y --survey $survey  --use_map_veto _HPmapcut --version $1
#16 minutes
python $LSS/LSS/scripts/main/mkCat_main.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec $verspec --add_weight_zfail y --survey $survey  --use_map_veto _HPmapcut --version $1
#25 minutes, most of time because of multiple writes
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type ELG_LOP --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec $verspec --add_weight_zfail y --survey $survey  --use_map_veto _HPmapcut --version $1
#7 minutes, most of time because of multiple writes
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type ELG_VLO --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec $verspec --add_weight_zfail y --survey $survey  --use_map_veto _HPmapcut --version $1
#13 minutes, most of time because of multiple writes
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec $verspec --add_weight_zfail y --survey $survey  --use_map_veto _HPmapcut --version $1
#10 minutes, most of time because of multiple writes
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type BGS_FAINT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec $verspec --add_weight_zfail y --survey $survey  --use_map_veto _HPmapcut --version $1

#clustering catalogs; 85 minutes total
#11 minutes
srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSS/LSS/scripts/main/mkCat_main.py --type QSO  --fulld n --survey $survey --verspec $verspec --clusd y --clusran y --splitGC y --nz y --par y  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --version $1 --extra_clus_dir 'nonKP/' --redo_fracz y --nearestneighbor y --des_resamp #--imsys_colname WEIGHT_IMLIN
#10 minutes
srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSS/LSS/scripts/main/mkCat_main.py --type LRG  --fulld n --survey $survey --verspec $verspec --clusd y --clusran y --splitGC y --nz y --par y  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --version $1 --extra_clus_dir 'nonKP/' --redo_fracz y --nearestneighbor y --des_resamp #--imsys_colname WEIGHT_IMLIN
#4 minutes
srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSS/LSS/scripts/main/mkCat_main.py --type LGE  --fulld n --survey $survey --verspec $verspec --clusd y --clusran y --splitGC y --nz y --par y  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --version $1 --extra_clus_dir 'nonKP/' --redo_fracz y --nearestneighbor y --des_resamp #--imsys_colname WEIGHT_IMLIN
#13 minutes
srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSS/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT  --fulld n --survey $survey --verspec $verspec --clusd y --clusran y --splitGC y --nz y --par y  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --version $1 --extra_clus_dir 'nonKP/' --redo_fracz y --nearestneighbor y --des_resamp #--imsys_colname WEIGHT_IMLIN
#13 minutes
srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSS/LSS/scripts/main/mkCat_main.py --type BGS_FAINT  --fulld n --survey $survey --verspec $verspec --clusd y --clusran y --splitGC y --nz y --par y  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --version $1 --extra_clus_dir 'nonKP/' --redo_fracz y --nearestneighbor y --des_resamp #--imsys_colname WEIGHT_IMLIN
#15 minutes
srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSS/LSS/scripts/main/mkCat_main.py --type ELG_LOP --notqso y  --fulld n --survey $survey --verspec $verspec --clusd y --clusran y --splitGC y --nz y --par y  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --version $1 --extra_clus_dir 'nonKP/' --redo_fracz y --nearestneighbor y --des_resamp #--imsys_colname WEIGHT_IMLIN
#19 minutes (why so long?)
srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSS/LSS/scripts/main/mkCat_main.py --type ELG_VLO --notqso y  --fulld n --survey $survey --verspec $verspec --clusd y --clusran y --splitGC y --nz y --par y  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --version $1 --extra_clus_dir 'nonKP/' --redo_fracz y --nearestneighbor y --des_resamp #--imsys_colname WEIGHT_IMLIN

#imaging systematics for non-ELG samples; 47 minutes total
#6 minutes
srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSS/LSS/scripts/addsys2clus.py --type LGE --basedir /global/cfs/cdirs/desi/survey/catalogs/ --version $1 --survey DA3 --extra_clus_dir nonKP --verspec matterhorn-v2 --doimlin y --imsys_zbin fine --replace_syscol --splitDES
#13 minutes
srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSS/LSS/scripts/addsys2clus.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/ --version $1 --survey DA3 --extra_clus_dir nonKP --verspec matterhorn-v2 --doimlin y --imsys_zbin fine --replace_syscol --splitDES
#10 minutes
srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSS/LSS/scripts/addsys2clus.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/ --version $1 --survey DA3 --extra_clus_dir nonKP --verspec matterhorn-v2 --doimlin y --imsys_zbin split --replace_syscol --splitDES
#10 minutes
srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSS/LSS/scripts/addsys2clus.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/ --version $1 --survey DA3 --extra_clus_dir nonKP --verspec matterhorn-v2 --doimlin y --imsys_zbin split --replace_syscol --splitDES
#8 minutes
srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSS/LSS/scripts/addsys2clus.py --type BGS_FAINT --basedir /global/cfs/cdirs/desi/survey/catalogs/ --version $1 --survey DA3 --extra_clus_dir nonKP --verspec matterhorn-v2 --doimlin y --imsys_zbin split --replace_syscol --splitDES