source /global/common/software/desi/desi_environment.sh main

export LSSDIR=$HOME

$LSSDIR/LSS/scripts/run_sysnet_ab2ndgen.sh N $2$3 true false 1024 0.003 dnnp pnll $1 
$LSSDIR/LSS/scripts/run_sysnet_ab2ndgen.sh S $2$3 true false 1024 0.003 dnnp pnll $1 

$LSSDIR/LSS/scripts/run_sysnet_ab2ndgen.sh N $2$3 false true 1024 0.004 dnnp pnll $1 
$LSSDIR/LSS/scripts/run_sysnet_ab2ndgen.sh S $2$3 false true 1024 0.004 dnnp pnll $1 
