#!/bin/bash

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
#srun -N 1 -C cpu -t 04:00:00 -q interactive  python $LSS/LSS/scripts/main/mkCat_main.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec $verspec --survey $survey --version $1 --mkgtl
#3 minutes
#srun -N 1 -C cpu -t 04:00:00 -q interactive  python $LSS/LSS/scripts/main/mkCat_main.py --type LGE --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec $verspec --survey $survey --version $1 --mkgtl
#13 minutes
#srun -N 1 -C cpu -t 04:00:00 -q interactive  python $LSS/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec $verspec --survey $survey --version $1 --mkgtl
#6.5 minutes
#srun -N 1 -C cpu -t 04:00:00 -q interactive  python $LSS/LSS/scripts/main/mkCat_main.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec $verspec --survey $survey --version $1
#25 minutes
#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type ELG_LOP --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec $verspec --survey $survey --version $1
#10 minutes
#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type ELG_VLO --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec $verspec --survey $survey --version $1
#7.5 minutes
#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type BGS_FAINT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec $verspec --survey $survey --version $1

#randoms with spec info
#only should be done for each spec release 
#~1hr 12 minutes
#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec $verspec --type dark --combwspec y --fullr n --survey $survey --maxr 18 --version $1 --nproc 9
#x minutes
#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec $verspec --type dark1b --combwspec y --fullr n --survey $survey --maxr 18 --version $1 --nproc 9

#x minutes
#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec $verspec --type bright --combwspec y --fullr n --survey $survey --maxr 18 --version $1 --nproc 9


#run this for new version of LSS catalogs
#59 minutes
#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec $verspec --type dark --combwspec n --fullr y --survey $survey --maxr 18 --version $1 --nproc 9

#15 minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec $verspec --type dark1b --combwspec n --fullr y --survey $survey --maxr 18 --version $1 --nproc 9


srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec $verspec --type bright --combwspec n --fullr y --survey $survey --maxr 18 --version $1 --nproc 9

#add LRG veto column
#x minutes 
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --add_veto y --verspec $verspec --survey $survey --maxr 18 --version $1

srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type LGE --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --add_veto y --verspec $verspec --survey $survey --maxr 18 --version $1 --par y

#fill randoms with properties and apply vetos
#x minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --fillran y --apply_veto y  --maxr 18 --version $1
#x minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --type LGE --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --fillran y --apply_veto y  --maxr 18 --version $1

#x minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --apply_veto y --maxr 18 --version $1

#x minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --type ELG_LOP --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --apply_veto y --maxr 18 --version $1
#x minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --type ELG_VLO --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --apply_veto y --maxr 18 --version $1

#x minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --fillran y --apply_veto y  --maxr 18 --version $1

#x minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --type BGS_FAINT --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --fillran y --apply_veto y  --maxr 18 --version $1


#make healpix maps
#52 minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive  python $LSS/LSS/scripts/main/mkCat_main.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec $verspec --survey $survey  --mkHPmaps y --version $1
#x minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive  python $LSS/LSS/scripts/main/mkCat_main.py --type LGE --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec $verspec --survey $survey  --mkHPmaps y --version $1

#49 minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive  python $LSS/LSS/scripts/main/mkCat_main.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec $verspec --survey $survey  --mkHPmaps y --version $1

#only one of the ELGs and BGS should be necessary
#47 minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive  python $LSS/LSS/scripts/main/mkCat_main.py --type ELG_LOP --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec $verspec --survey $survey  --mkHPmaps y --version $1

#x minutes
srun -N 1 -C cpu -t 04:00:00 -q interactive  python $LSS/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec $verspec --survey $survey  --mkHPmaps y --version $1

#apply vetos to data, get FRAC_TLOBS info

#5 minutes
python $LSS/LSS/scripts/main/mkCat_main.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec $verspec --survey $survey --maxr 0 --version $1
#5 minutes
python $LSS/LSS/scripts/main/mkCat_main.py --type LGE --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec $verspec --survey $survey --maxr 0 --version $1

#3 minutes
python $LSS/LSS/scripts/main/mkCat_main.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec $verspec --survey $survey --maxr 0 --version $1

srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type ELG_LOP --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec $verspec --survey $survey --maxr 0 --version $1

srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type ELG_VLO --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec $verspec --survey $survey --maxr 0 --version $1

srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec $verspec --survey $survey --maxr 0 --version $1

srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type BGS_FAINT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec $verspec --survey $survey --maxr 0 --version $1

#Add FRAC_TLOBS info to randoms

srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --add_tl y --maxr 18 --version $1

srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --add_tl y --maxr 18 --version $1

srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --add_tl y --maxr 18 --version $1

srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --type BGS_FAINT --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --add_tl y --maxr 18 --version $1


srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --type LGE --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --add_tl y --maxr 18 --version $1

srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --type ELG_LOP --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --add_tl y --maxr 18 --version $1

srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main_ran.py --type ELG_VLO --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/   --verspec $verspec --survey $survey --add_tl y --maxr 18 --version $1


#Apply map veto

srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_map_veto y --verspec $verspec --survey $survey --version $1

srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_map_veto y --verspec $verspec --survey $survey --version $1

srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type BGS_FAINT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_map_veto y --verspec $verspec --survey $survey --version $1


srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type LGE --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_map_veto y --verspec $verspec --survey $survey --version $1

srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_map_veto y --verspec $verspec --survey $survey --version $1
#srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_map_veto y --verspec $verspec --survey $survey --version $1 --maxr 0

srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type ELG_LOP --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_map_veto y --verspec $verspec --survey $survey --version $1

srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type ELG_VLO --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_map_veto y --verspec $verspec --survey $survey --version $1

#redshift failure modeling
python $LSS/LSS/scripts/main/mkCat_main.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec $verspec --add_weight_zfail y --survey $survey  --use_map_veto _HPmapcut --version $1

python $LSS/LSS/scripts/main/mkCat_main.py --type LGE --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec $verspec --add_weight_zfail y --survey $survey  --use_map_veto _HPmapcut --version $1

python $LSS/LSS/scripts/main/mkCat_main.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec $verspec --add_weight_zfail y --survey $survey  --use_map_veto _HPmapcut --version $1

#25 minutes, most of time because of multiple writes
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type ELG_LOP --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec $verspec --add_weight_zfail y --survey $survey  --use_map_veto _HPmapcut --version $1

#x minutes, most of time because of multiple writes
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type ELG_VLO --notqso y --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec $verspec --add_weight_zfail y --survey $survey  --use_map_veto _HPmapcut --version $1

#x minutes, most of time because of multiple writes
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec $verspec --add_weight_zfail y --survey $survey  --use_map_veto _HPmapcut --version $1

#x minutes, most of time because of multiple writes
srun -N 1 -C cpu -t 04:00:00 -q interactive python $LSS/LSS/scripts/main/mkCat_main.py --type BGS_FAINT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --verspec $verspec --add_weight_zfail y --survey $survey  --use_map_veto _HPmapcut --version $1

#clustering catalogs

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSS/LSS/scripts/main/mkCat_main.py --type QSO  --fulld n --survey $survey --verspec $verspec --clusd y --clusran y --splitGC y --nz y --par y  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --version $1 --extra_clus_dir 'nonKP/' --redo_fracz y --nearestneighbor y --des_resamp #--imsys_colname WEIGHT_IMLIN

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSS/LSS/scripts/main/mkCat_main.py --type LRG  --fulld n --survey $survey --verspec $verspec --clusd y --clusran y --splitGC y --nz y --par y  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --version $1 --extra_clus_dir 'nonKP/' --redo_fracz y --nearestneighbor y --des_resamp #--imsys_colname WEIGHT_IMLIN

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSS/LSS/scripts/main/mkCat_main.py --type LGE  --fulld n --survey $survey --verspec $verspec --clusd y --clusran y --splitGC y --nz y --par y  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --version $1 --extra_clus_dir 'nonKP/' --redo_fracz y --nearestneighbor y --des_resamp #--imsys_colname WEIGHT_IMLIN

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSS/LSS/scripts/main/mkCat_main.py --type BGS_BRIGHT  --fulld n --survey $survey --verspec $verspec --clusd y --clusran y --splitGC y --nz y --par y  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --version $1 --extra_clus_dir 'nonKP/' --redo_fracz y --nearestneighbor y --des_resamp #--imsys_colname WEIGHT_IMLIN

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSS/LSS/scripts/main/mkCat_main.py --type BGS_FAINT  --fulld n --survey $survey --verspec $verspec --clusd y --clusran y --splitGC y --nz y --par y  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --version $1 --extra_clus_dir 'nonKP/' --redo_fracz y --nearestneighbor y --des_resamp #--imsys_colname WEIGHT_IMLIN

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSS/LSS/scripts/main/mkCat_main.py --type ELG_LOP --notqso y  --fulld n --survey $survey --verspec $verspec --clusd y --clusran y --splitGC y --nz y --par y  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --version $1 --extra_clus_dir 'nonKP/' --redo_fracz y --nearestneighbor y --des_resamp #--imsys_colname WEIGHT_IMLIN

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSS/LSS/scripts/main/mkCat_main.py --type ELG_VLO --notqso y  --fulld n --survey $survey --verspec $verspec --clusd y --clusran y --splitGC y --nz y --par y  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --version $1 --extra_clus_dir 'nonKP/' --redo_fracz y --nearestneighbor y --des_resamp #--imsys_colname WEIGHT_IMLIN

#imaging systematics for non-ELG samples

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSS/LSS/scripts/addsys2clus.py --type LGE --basedir /global/cfs/cdirs/desi/survey/catalogs/ --version $1 --survey DA3 --extra_clus_dir nonKP --verspec matterhorn-v2 --doimlin y --imsys_zbin fine --replace_syscol --splitDES

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSS/LSS/scripts/addsys2clus.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/ --version $1 --survey DA3 --extra_clus_dir nonKP --verspec matterhorn-v2 --doimlin y --imsys_zbin fine --replace_syscol --splitDES

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSS/LSS/scripts/addsys2clus.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/ --version $1 --survey DA3 --extra_clus_dir nonKP --verspec matterhorn-v2 --doimlin y --imsys_zbin split --replace_syscol --splitDES

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSS/LSS/scripts/addsys2clus.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/ --version $1 --survey DA3 --extra_clus_dir nonKP --verspec matterhorn-v2 --doimlin y --imsys_zbin split --replace_syscol --splitDES

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSS/LSS/scripts/addsys2clus.py --type BGS_FAINT --basedir /global/cfs/cdirs/desi/survey/catalogs/ --version $1 --survey DA3 --extra_clus_dir nonKP --verspec matterhorn-v2 --doimlin y --imsys_zbin split --replace_syscol --splitDES