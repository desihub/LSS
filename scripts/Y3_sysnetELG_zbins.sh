#!/bin/bash
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main

export LSSDIR=$HOME
#export SYSNETDIR=$HOME/desicode
export LSSBASE=/global/cfs/cdirs/desi/survey/catalogs/
PYTHONPATH=$PYTHONPATH:$LSSDIR/LSS

verspec=loa-v1
survey=DA2

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/main/mkCat_main.py --basedir $LSSBASE --type ELG_LOP --notqso y --prepsysnet y --imsys_zbin y --fulld n --survey $survey --verspec $verspec --version $1

$LSSDIR/LSS/scripts/run_sysnet_cd.sh N ELG_LOPnotqso0.8_1.1 true false 1024 0.003 dnnp pnll $1 $LSSBASE/$survey/LSS/$verspec/LSScats/
$LSSDIR/LSS/scripts/run_sysnet_cd.sh S ELG_LOPnotqso0.8_1.1 true false 1024 0.003 dnnp pnll $1 $LSSBASE/$survey/LSS/$verspec/LSScats/

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi $LSSDIR/LSS/scripts/run_sysnet_cd.sh N ELG_LOPnotqso0.8_1.1 false true 1024 0.004 dnnp pnll $1 $LSSBASE/$survey/LSS/$verspec/LSScats/
srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi $LSSDIR/LSS/scripts/run_sysnet_cd.sh S ELG_LOPnotqso0.8_1.1 false true 1024 0.004 dnnp pnll $1 $LSSBASE/$survey/LSS/$verspec/LSScats/

$LSSDIR/LSS/scripts/run_sysnet_cd.sh N ELG_LOPnotqso1.1_1.6 true false 1024 0.003 dnnp pnll $1 $LSSBASE/$survey/LSS/$verspec/LSScats/
$LSSDIR/LSS/scripts/run_sysnet_cd.sh S ELG_LOPnotqso1.1_1.6 true false 1024 0.003 dnnp pnll $1 $LSSBASE/$survey/LSS/$verspec/LSScats/

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi $LSSDIR/LSS/scripts/run_sysnet_cd.sh N ELG_LOPnotqso1.1_1.6 false true 1024 0.004 dnnp pnll $1 $LSSBASE/$survey/LSS/$verspec/LSScats/
srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi $LSSDIR/LSS/scripts/run_sysnet_cd.sh S ELG_LOPnotqso1.1_1.6 false true 1024 0.004 dnnp pnll $1 $LSSBASE/$survey/LSS/$verspec/LSScats/


python scripts/main/mkCat_main.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --type ELG_LOP --notqso y --add_sysnet y --imsys_zbin y --fulld n --survey $survey --verspec $verspec --version $1

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/validation/validation_improp_full.py --tracers ELG_LOPnotqso --version $1 --verspec $verspec --survey $survey --weight_col WEIGHT_SN
