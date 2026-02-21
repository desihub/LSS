#!/bin/bash
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main

export LSSDIR=$HOME
#export SYSNETDIR=$HOME/desicode
export LSSBASE=/global/cfs/cdirs/desi/survey/catalogs/
PYTHONPATH=$PYTHONPATH:$LSSDIR/LSS

python scripts/main/mkCat_main.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --type ELG_LOP --notqso y --prepsysnet y --imsys_zbin y --fulld n --survey Y1 --verspec iron --version $1

$LSSDIR/LSS/scripts/run_sysnet_cd.sh N ELG_LOPnotqso0.8_1.1 true false 1024 0.003 dnnp pnll $1 $LSSBASE/Y1/LSS/iron/LSScats/
$LSSDIR/LSS/scripts/run_sysnet_cd.sh S ELG_LOPnotqso0.8_1.1 true false 1024 0.003 dnnp pnll $1 $LSSBASE/Y1/LSS/iron/LSScats/

$LSSDIR/LSS/scripts/run_sysnet_cd.sh N ELG_LOPnotqso0.8_1.1 false true 1024 0.004 dnnp pnll $1 $LSSBASE/Y1/LSS/iron/LSScats/
$LSSDIR/LSS/scripts/run_sysnet_cd.sh S ELG_LOPnotqso0.8_1.1 false true 1024 0.004 dnnp pnll $1 $LSSBASE/Y1/LSS/iron/LSScats/

$LSSDIR/LSS/scripts/run_sysnet_cd.sh N ELG_LOPnotqso1.1_1.6 true false 1024 0.003 dnnp pnll $1 $LSSBASE/Y1/LSS/iron/LSScats/
$LSSDIR/LSS/scripts/run_sysnet_cd.sh S ELG_LOPnotqso1.1_1.6 true false 1024 0.003 dnnp pnll $1 $LSSBASE/Y1/LSS/iron/LSScats/

$LSSDIR/LSS/scripts/run_sysnet_cd.sh N ELG_LOPnotqso1.1_1.6 false true 1024 0.004 dnnp pnll $1 $LSSBASE/Y1/LSS/iron/LSScats/
$LSSDIR/LSS/scripts/run_sysnet_cd.sh S ELG_LOPnotqso1.1_1.6 false true 1024 0.004 dnnp pnll $1 $LSSBASE/Y1/LSS/iron/LSScats/


python scripts/main/mkCat_main.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --type ELG_LOP --notqso y --add_sysnet y --imsys_zbin y --fulld n --survey Y1 --verspec iron --version $1

python scripts/validation/validation_improp_full.py --tracer ELG_LOPnotqso --version $1

