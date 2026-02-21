#!/bin/bash
# bash Y1_sysnetLRG_zbins_v1.4.1.sh
set -e
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main

export LSSDIR=$HOME
#export SYSNETDIR=$HOME/desicode
export LSSBASE=/global/cfs/cdirs/desi/survey/catalogs/
PYTHONPATH=$PYTHONPATH:$LSSDIR/LSS

VERSION=v1.4.1
TRACER=LRG
do_prep=false
do_lr=false
do_SN=false
add_weights=false
do_validation=true

# Some NN parameters for North
LR_N1=0.01     # learning rate for LRG1 N
LR_N2=0.008    # learning rate for LRG2 N
LR_N3=0.007    # learning rate for LRG3 N
NBATCH_N=256   # Powers of 2
NCHAIN_N=4     # chains
NEPOCH_N=70    # number of epochs
NNS_N=(3 20)   # NN structure (# layers, # units)

# Some NN parameters for South
LR_S1=0.008    # learning rate for LRG1 S
LR_S2=0.007    # learning rate for LRG2 S
LR_S3=0.007    # learning rate for LRG3 S
NBATCH_S=1024  # Powers of 2
NCHAIN_S=4     # chains
NEPOCH_S=70    # number of epochs
NNS_S=(3 20)   # NN structure (# layers, # units)

BASEDIR=$LSSBASE/Y1/LSS/iron/LSScats/
RUN_SYSNET=$LSSDIR/LSS/scripts/run_sysnetELG_cd_new.sh

if [ $do_prep == true ]
then
    # If using the allsky randoms option when preparing for sysnet mkCat_main.py might not work without salloc or job
    python $LSSDIR/LSS/scripts/main/mkCat_main.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --type LRG --prepsysnet y --imsys_zbin y --fulld n --survey Y1 --verspec iron --version $VERSION --use_allsky_rands y --par y
fi

if [ $do_lr == true ]
then
    # Find learning rate for North
    $RUN_SYSNET N LRG0.8_1.1 true false $NBATCH_N 0.003 dnnp pnll $VERSION $BASEDIR $NCHAIN_N $NEPOCH_N ${NNS_N[@]}
    $RUN_SYSNET N LRG0.6_0.8 true false $NBATCH_N 0.003 dnnp pnll $VERSION $BASEDIR $NCHAIN_N $NEPOCH_N ${NNS_N[@]}
    $RUN_SYSNET N LRG0.4_0.6 true false $NBATCH_N 0.003 dnnp pnll $VERSION $BASEDIR $NCHAIN_N $NEPOCH_N ${NNS_N[@]}
    # Find learning rate for South
    $RUN_SYSNET S LRG0.8_1.1 true false $NBATCH_S 0.003 dnnp pnll $VERSION $BASEDIR $NCHAIN_S $NEPOCH_S ${NNS_N[@]}
    $RUN_SYSNET S LRG0.6_0.8 true false $NBATCH_S 0.003 dnnp pnll $VERSION $BASEDIR $NCHAIN_S $NEPOCH_S ${NNS_N[@]}
    $RUN_SYSNET S LRG0.4_0.6 true false $NBATCH_S 0.003 dnnp pnll $VERSION $BASEDIR $NCHAIN_S $NEPOCH_S ${NNS_N[@]}
fi

if [ $do_SN == true ]
then
    # Run SYSNet for North
    srun -N 1 $RUN_SYSNET N LRG0.8_1.1 false true $NBATCH_N $LR_N3 dnnp pnll $VERSION $BASEDIR $NCHAIN_N $NEPOCH_N ${NNS_N[@]} &
    srun -N 1 $RUN_SYSNET N LRG0.6_0.8 false true $NBATCH_N $LR_N2 dnnp pnll $VERSION $BASEDIR $NCHAIN_N $NEPOCH_N ${NNS_N[@]} &
    srun -N 1 $RUN_SYSNET N LRG0.4_0.6 false true $NBATCH_N $LR_N1 dnnp pnll $VERSION $BASEDIR $NCHAIN_N $NEPOCH_N ${NNS_N[@]} &
    wait
    # Run SYSNet for South
    srun -N 1 $RUN_SYSNET S LRG0.8_1.1 false true $NBATCH_S $LR_S3 dnnp pnll $VERSION $BASEDIR $NCHAIN_S $NEPOCH_S ${NNS_S[@]} &
    srun -N 1 $RUN_SYSNET S LRG0.6_0.8 false true $NBATCH_S $LR_S2 dnnp pnll $VERSION $BASEDIR $NCHAIN_S $NEPOCH_S ${NNS_S[@]} &
    srun -N 1 $RUN_SYSNET S LRG0.4_0.6 false true $NBATCH_S $LR_S1 dnnp pnll $VERSION $BASEDIR $NCHAIN_S $NEPOCH_S ${NNS_S[@]} &
    wait
fi

if [ $add_weights == true ]
then
    python $LSSDIR/LSS/scripts/main/mkCat_main.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --type $TRACER --add_sysnet y --imsys_zbin y --fulld n --survey Y1 --verspec iron --version $VERSION --imsys_colname WEIGHT_SN
fi

if [ $do_validation == true ]
then
    python $LSSDIR/LSS/scripts/validation/validation_improp_full.py --tracer $TRACER --version $VERSION --weight_col WEIGHT_SN
fi
