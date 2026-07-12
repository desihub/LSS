#!/bin/bash

set -e

source /global/common/software/desi/desi_environment.sh main
module load LSS/main
#module swap desitarget/3.0.0
#export LSSCODE=$HOME ; do this before script, e.g., export LSSCODE=$HOME/LSScode for desica
#PYTHONPATH=$PYTHONPATH:$LSSCODE/LSS/py
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
verspec=loa-v1
survey=DA2
LSSCODE=$HOME/LSScode/LSS/
version=v2
edir=CMBLENS
scriptdir=$LSSCODE/scripts
bdir=/global/cfs/cdirs/desi/survey/catalogs/
#PYTHONPATH=$PYTHONPATH:$LSSCODE/LSS/py

#srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $scriptdir/addsys2clus.py --type ELG_LOPnotqso   --basedir  $bdir   --prep4sysnet y --survey $survey --verspec $verspec --imsys_zbin split  --version $version --extra_clus_dir $edir --compwtmd fraczNN

$scriptdir/sysnetELG_LOPnotqso_zbins.sh '' ELG_LOPnotqso $bdir/$survey/LSS/$verspec/LSScats/$version/$edir/

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $scriptdir/addsys2clus.py --type ELG_LOPnotqso   --basedir $bdir    --addsysnet y --survey $survey --verspec $verspec --imsys_zbin split  --version $version --extra_clus_dir edir --replace_syscol

