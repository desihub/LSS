source /global/common/software/desi/desi_environment.sh main

export LSSDIR=$HOME
export SYSNETDIR=$HOME/desicode
export LSSBASE=/global/cfs/cdirs/desi/survey/catalogs/
PYTHONPATH=$PYTHONPATH:$LSSDIR/LSS/py

python scripts/main/mkCat_main.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --type ELG_LOP --notqso y --prepsysnet y --fulld n --survey Y1 --verspec iron --version $1

$LSSDIR/LSS/scripts/run_sysnet.sh N $2$3 true false 1024 0.003 dnnp pnll $1 $LSSBASE/Y1/LSS/iron/LSScats/
$LSSDIR/LSS/scripts/run_sysnet.sh S $2$3 true false 1024 0.003 dnnp pnll $1 $LSSBASE/Y1/LSS/iron/LSScats/

$LSSDIR/LSS/scripts/run_sysnet.sh N $2$3 false true 1024 0.004 dnnp pnll $1 $LSSBASE/Y1/LSS/iron/LSScats/
$LSSDIR/LSS/scripts/run_sysnet.sh S $2$3 false true 1024 0.004 dnnp pnll $1 $LSSBASE/Y1/LSS/iron/LSScats/
