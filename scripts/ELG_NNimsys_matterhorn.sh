#!/bin/bash

#run this after a node has already been allocated
set -e

source /global/common/software/desi/desi_environment.sh main
module load LSS/main
#module swap desitarget/3.0.0
#export LSSCODE=$HOME ; do this before script, e.g., export LSSCODE=$HOME/LSScode for desica
#PYTHONPATH=$PYTHONPATH:$LSSCODE/LSS/py
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
verspec=matterhorn-v2
survey=DA3
export LSSCODE=$HOME/LSScode/LSS/ #necessary for sysnet code at present
version=$1
edir=nonKP
scriptdir=$LSSCODE/scripts
bdir=/global/cfs/cdirs/desi/survey/catalogs/
#PYTHONPATH=$PYTHONPATH:$LSSCODE/LSS/py

python $scriptdir/addsys2clus.py --type ELG_LOPnotqso   --basedir  $bdir   --prep4sysnet y --survey $survey --verspec $verspec --imsys_zbin split  --version $version --extra_clus_dir $edir --compwtmd fraczNN --splitDES

python $scriptdir/addsys2clus.py --type ELG_LOPnotqso   --basedir  $bdir   --prep4sysnet y --survey $survey --verspec $verspec --imsys_zbin split  --version $version --extra_clus_dir $edir --compwtmd fraczNN --splitDES


$scriptdir/sysnet_splitDES_zbins.sh '' ELG_LOPnotqso $bdir/$survey/LSS/$verspec/LSScats/$version/$edir/

$scriptdir/sysnet_splitDES_zbins.sh '' ELG_VLOnotqso $bdir/$survey/LSS/$verspec/LSScats/$version/$edir/


python $scriptdir/addsys2clus.py --type ELG_LOPnotqso   --basedir $bdir    --addsysnet y --survey $survey --verspec $verspec --imsys_zbin split  --version $version --extra_clus_dir $edir --replace_syscol --splitDES

python $scriptdir/addsys2clus.py --type ELG_VLOnotqso   --basedir $bdir    --addsysnet y --survey $survey --verspec $verspec --imsys_zbin split  --version $version --extra_clus_dir $edir --replace_syscol --splitDES


