#!/bin/bash
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main

export LSSDIR=$HOME
#export SYSNETDIR=$HOME/desicode
export LSSBASE=/global/cfs/cdirs/desi/survey/catalogs/
PYTHONPATH=$PYTHONPATH:$LSSDIR/LSS

VERSION=$1

# Some NN parameters for North
LR_N=0.02    # learning rate
NBATCH_N=256  # Powers of 2
NCHAIN_N=5    # chains
NEPOCH_N=100  # number of epochs
NNS_N=(3 10)  # NN structure (# layers, # units)

# Some NN parameters for South
LR_S=0.01    # learning rate
NBATCH_S=1024 # Powers of 2
NCHAIN_S=5    # chains
NEPOCH_S=100  # number of epochs
NNS_S=(4 20)  # NN structure (# layers, # units)

BASEDIR=$LSSBASE/Y1/LSS/iron/LSScats/
RUN_SYSNET=$LSSDIR/LSS/scripts/run_sysnetELG_cd_new.sh

# If using the allsky randoms option when preparing for sysnet mkCat_main.py might not work without salloc or job; is run if supplied argument 
if [ $2 == 1 ]; then
  python scripts/main/mkCat_main.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --type BGS_BRIGHT-21.5 --prepsysnet y --imsys_zbin y --fulld n --survey Y1 --verspec iron --version $VERSION --use_allsky_rands y
fi

# Find learning rate for North
$RUN_SYSNET N BGS_BRIGHT-21.50.1_0.4 true false $NBATCH_N 0.003 dnnp pnll $VERSION $BASEDIR $NCHAIN_N $NEPOCH_N ${NNS_N[@]}
# Find learning rate for South
$RUN_SYSNET S BGS_BRIGHT-21.50.1_0.4 true false $NBATCH_S 0.003 dnnp pnll $VERSION $BASEDIR $NCHAIN_S $NEPOCH_S ${NNS_S[@]}

# Run SYSNet for North
$RUN_SYSNET N BGS_BRIGHT-21.50.1_0.4 false true $NBATCH_N $LR_N dnnp pnll $VERSION $BASEDIR $NCHAIN_N $NEPOCH_N ${NNS_N[@]}
# Run SYSNet for South
$RUN_SYSNET S BGS_BRIGHT-21.50.1_0.4 false true $NBATCH_S $LR_S dnnp pnll $VERSION $BASEDIR $NCHAIN_S $NEPOCH_S ${NNS_S[@]}

python scripts/main/mkCat_main.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --type BGS_BRIGHT-21.5 --add_sysnet y --imsys_zbin y --fulld n --survey Y1 --verspec iron --version $VERSION

python scripts/validation/validation_improp_full.py --tracer BGS_BRIGHT-21.5 --version $VERSION --weight_col WEIGHT_SN

