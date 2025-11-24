#!/bin/bash
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh test

#inputs should be catalog version (could be '') and the directory up to the version

#export LSSDIR=$HOME
#export SYSNETDIR=$HOME/desicode
#export LSSBASE=/global/cfs/cdirs/desi/survey/catalogs/
#PYTHONPATH=$PYTHONPATH:$LSSDIR/LSS/py

#verspec=loa-v1
#survey=DA2


$LSSCODE/LSS/scripts/run_sysnet_cd.sh N ELG_LOPnotqso0.8_1.1 true false 1024 0.003 dnnp pnll $1 $2
$LSSCODE/LSS/scripts/run_sysnet_cd.sh S ELG_LOPnotqso0.8_1.1 true false 1024 0.003 dnnp pnll $1 $2

$LSSCODE/LSS/scripts/run_sysnet_cd.sh N ELG_LOPnotqso0.8_1.1 false true 1024 0.004 dnnp pnll $1 $2
$LSSCODE/LSS/scripts/run_sysnet_cd.sh S ELG_LOPnotqso0.8_1.1 false true 1024 0.004 dnnp pnll $1 $2

$LSSCODE/LSS/scripts/run_sysnet_cd.sh N ELG_LOPnotqso1.1_1.6 true false 1024 0.003 dnnp pnll $1 $2
$LSSCODE/LSS/scripts/run_sysnet_cd.sh S ELG_LOPnotqso1.1_1.6 true false 1024 0.003 dnnp pnll $1 $2

$LSSCODE/LSS/scripts/run_sysnet_cd.sh N ELG_LOPnotqso1.1_1.6 false true 1024 0.004 dnnp pnll $1 $2
$LSSCODE/LSS/scripts/run_sysnet_cd.sh S ELG_LOPnotqso1.1_1.6 false true 1024 0.004 dnnp pnll $1 $2

