#!/bin/bash

set -e

source /global/common/software/desi/desi_environment.sh main
module load LSS/main
#module swap desitarget/3.0.0
#export LSSCODE=$HOME ; do this before script, e.g., export LSSCODE=$HOME/LSScode for desica
#PYTHONPATH=$PYTHONPATH:$LSSCODE/LSS/py
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
verspec=matterhorn-v2
survey=DA3
LSSCODE=$HOME/LSScode
#PYTHONPATH=$PYTHONPATH:$LSSCODE/LSS/py

#srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/addsys2clus.py --type QSO   --basedir /global/cfs/cdirs/desi/survey/catalogs/    --doimlin y --survey $survey --verspec $verspec --imsys_zbin split  --version test --extra_clus_dir nonKP --replace_syscol

#srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/addsys2clus.py --type LRG   --basedir /global/cfs/cdirs/desi/survey/catalogs/    --doimlin y --survey $survey --verspec $verspec --imsys_zbin fine  --version test --extra_clus_dir nonKP --replace_syscol

#srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/addsys2clus.py --type ELG_LOPnotqso   --basedir /global/cfs/cdirs/desi/survey/catalogs/    --prep4sysnet y --survey $survey --verspec $verspec --imsys_zbin split  --version test --extra_clus_dir nonKP

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/addsys2clus.py --type ELG_LOPnotqso   --basedir /global/cfs/cdirs/desi/survey/catalogs/    --addsysnet y --survey $survey --verspec $verspec --imsys_zbin split  --version test --extra_clus_dir nonKP --replace_syscol

#srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi $LSSCODE/LSS/scripts/Y3_sysnetELG_zbins_new.sh $1

#srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSSCODE/LSS/scripts/main/mkCat_main.py --type ELG_LOP --notqso y  --fulld n --survey $survey --verspec $verspec --clusd y --clusran y --splitGC y --nz y --par y  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --version $1 --extra_clus_dir 'nonKP/' #--imsys_colname WEIGHT_SN

#srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSSCODE/LSS/scripts/main/mkCat_main.py --type ELG --notqso y  --fulld n --survey $survey --verspec $verspec --clusd y --clusran y --splitGC y --nz y --par y  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --version $1 --extra_clus_dir 'nonKP/' #--imsys_colname WEIGHT_SN


#mkdir /global/cfs/cdirs/desi/survey/catalogs/$survey/LSS/$verspec/LSScats/$1/nonKP/
#mv /global/cfs/cdirs/desi/survey/catalogs/$survey/LSS/$verspec/LSScats/$1/*clustering* /global/cfs/cdirs/desi/survey/catalogs/$survey/LSS/$verspec/LSScats/$1/nonKP/
