
source /global/common/software/desi/desi_environment.sh main

export LSSDIR=$HOME
export SYSNETDIR=$HOME/desicode
export LSSBASE=/global/cfs/cdirs/desi/survey/catalogs/
PYTHONPATH=$PYTHONPATH:$LSSDIR/LSS

python scripts/main/mkCat_main.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --type ELG_LOP --notqso y --prepsysnet y --fulld n --survey Y1 --verspec iron --version $1

$LSSDIR/LSS/scripts/run_sysnet.sh N ELG_LOPnotqso true false 1024 0.003 dnnp pnll $1 $LSSBASE/Y1/LSS/iron/LSScats/
$LSSDIR/LSS/scripts/run_sysnet.sh S ELG_LOPnotqso true false 1024 0.003 dnnp pnll $1 $LSSBASE/Y1/LSS/iron/LSScats/

$LSSDIR/LSS/scripts/run_sysnet.sh N ELG_LOPnotqso false true 1024 0.004 dnnp pnll $1 $LSSBASE/Y1/LSS/iron/LSScats/
$LSSDIR/LSS/scripts/run_sysnet.sh S ELG_LOPnotqso false true 1024 0.004 dnnp pnll $1 $LSSBASE/Y1/LSS/iron/LSScats/

python scripts/main/mkCat_main.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --type ELG_LOP --notqso y --add_sysnet y --fulld n --survey Y1 --verspec iron --version $1

python scripts/validation/validation_sky.py --tracer ELG_LOPnotqso --version $1

python scripts/main/mkCat_main.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --type LRG --prepsysnet y --fulld n --survey Y1 --verspec iron --version $1

$LSSDIR/LSS/scripts/run_sysnet.sh N LRG true false 1024 0.003 dnnp pnll $1 $LSSBASE/Y1/LSS/iron/LSScats/
$LSSDIR/LSS/scripts/run_sysnet.sh S LRG true false 1024 0.003 dnnp pnll $1 $LSSBASE/Y1/LSS/iron/LSScats/

$LSSDIR/LSS/scripts/run_sysnet.sh N LRG false true 1024 0.01 dnnp pnll $1 $LSSBASE/Y1/LSS/iron/LSScats/
$LSSDIR/LSS/scripts/run_sysnet.sh S LRG false true 1024 0.01 dnnp pnll $1 $LSSBASE/Y1/LSS/iron/LSScats/

python scripts/main/mkCat_main.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --type LRG --add_sysnet y --fulld n --survey Y1 --verspec iron --version $1

python scripts/main/mkCat_main.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --type QSO --prepsysnet y --fulld n --survey Y1 --verspec iron --version $1

$LSSDIR/LSS/scripts/run_sysnet.sh N QSO true false 1024 0.003 dnnp pnll $1 $LSSBASE/Y1/LSS/iron/LSScats/
$LSSDIR/LSS/scripts/run_sysnet.sh S QSO true false 1024 0.003 dnnp pnll $1 $LSSBASE/Y1/LSS/iron/LSScats/

$LSSDIR/LSS/scripts/run_sysnet.sh N QSO false true 1024 0.01 dnnp pnll $1 $LSSBASE/Y1/LSS/iron/LSScats/
$LSSDIR/LSS/scripts/run_sysnet.sh S QSO false true 1024 0.01 dnnp pnll $1 $LSSBASE/Y1/LSS/iron/LSScats/


python scripts/main/mkCat_main.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --type BGS_BRIGHT-21.5 --prepsysnet y --fulld n --survey Y1 --verspec iron --version $1

$LSSDIR/LSS/scripts/run_sysnet.sh N BGS_BRIGHT-21.5 true false 1024 0.003 dnnp pnll $1 $LSSBASE/Y1/LSS/iron/LSScats/
$LSSDIR/LSS/scripts/run_sysnet.sh S BGS_BRIGHT-21.5 true false 1024 0.003 dnnp pnll $1 $LSSBASE/Y1/LSS/iron/LSScats/

$LSSDIR/LSS/scripts/run_sysnet.sh N BGS_BRIGHT-21.5 false true 1024 0.01 dnnp pnll $1 $LSSBASE/Y1/LSS/iron/LSScats/
$LSSDIR/LSS/scripts/run_sysnet.sh S BGS_BRIGHT-21.5 false true 1024 0.01 dnnp pnll $1 $LSSBASE/Y1/LSS/iron/LSScats/
