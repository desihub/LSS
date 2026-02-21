#!/bin/bash
#SBATCH -N 1
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -J sys_Y1ab2ndgen
#SBATCH -t 1:00:00
#SBATCH -L SCRATCH
#SBATCH --output=/pscratch/sd/a/arosado/slurm/Y1ab2ndgen_regressis_%A_%a.out
#SBATCH --array=1-24
start_time=$(date +%s)
set -e
# This script can be run as a sbatch job, or an interactive job
# if run in an interactive node then manually set mock realization 
echo $SLURM_ARRAY_TASK_ID

source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main 
export LSSDIR=$HOME/LSS
PYTHONPATH=$PYTHONPATH:$LSSDIR

RUN_SYSNET=$LSSDIR/scripts/run_sysnetELG_ab2ndgen.sh 
#/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_v4_1//altmtl0/mock0/LSScats
# REGION only useful for ELGs right now
DATA_VERSION='v0.6'
REGION='all'
TRACER=ELG_LOPnotqso
#TRACER=LRG
MOCKVERSION=/SecondGenMocks/AbacusSummit_v4_1fixran/
BASEDIR=/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/$MOCKVERSION
MOCKCATVER='v0'
ALTMTL='y'

#REAL=$SLURM_ARRAY_TASK_ID
REAL=0
do_RF=false
do_prep=false
do_SN_ELG=false
do_SN_LRG=false
add_SN=true # if true make sure all weights are available, if not this will fail
do_validation_ffa=false
do_validation_altmtl=false
do_RF_validation_altmtl=false

#DIR=$BASEDIR/mock$REAL
DIR=$BASEDIR/altmtl$REAL/mock$REAL/LSScats

SUBSTRING=${TRACER:0:3}
echo $SUBSTRING

if [ $do_RF == true ]
then
    echo "doing RF regression and adding weights"
    # Compute and add regressis weights to mocks
    if [ $MOCKCATVER == 'v0' ]
    then
        srun -N 1 -n 1 python $LSSDIR/scripts/mock_tools/addsys.py --realization $REAL --tracer $TRACER --regressis y --add_regressis y --add_regressis_ran y --par y --base_dir $BASEDIR --use_altmtl $ALTMTL --data_version $DATA_VERSION
    else
        srun -N 1 -n 1 python $LSSDIR/scripts/mock_tools/addsys.py --realization $REAL --tracer $TRACER --regressis y --add_regressis y --add_regressis_ran y --par y --base_dir $BASEDIR --mockcatver $MOCKCATVER --use_altmtl $ALTMTL --data_version $DATA_VERSION
    fi
    echo "finished doing RF regression and adding RF weights"
fi

if [ $do_prep == true ]
then
    echo "preparing data for SYSNet"
    # prepare data for SYSNet
    if [ $MOCKCATVER == 'v0' ]
    then
        srun -N 1 -n 1 python $LSSDIR/scripts/mock_tools/addsys.py --tracer $TRACER --prepsysnet y --base_dir $BASEDIR --realization $REAL --use_allsky_rands 'y' --use_altmtl $ALTMTL --data_version $DATA_VERSION
    else
        srun -N 1 -n 1 python $LSSDIR/scripts/mock_tools/addsys.py --tracer $TRACER --prepsysnet y --base_dir $BASEDIR --realization $REAL --mockcatver $MOCKCATVER --use_allsky_rands 'y' --use_altmtl $ALTMTL  --data_version $DATA_VERSION
    fi
fi

if [ $do_SN_ELG == true ]
then
    # Run SYSNET different options for different zbins
    # Must have data prepared
    echo "running SYSNet on ${TRACER}"
    zbin1=0.8_1.1
    zbin2=1.1_1.6
    if [ $REGION == 'all' ]
    then 
        # Some NN parameters for North
        LR_N=0.009    # learning rate
        NBATCH_N=256  # Powers of 2
        NCHAIN_N=5    # chains
        NEPOCH_N=100  # number of epochs
        NNS_N=(3 10)  # NN structure (# layers, # units)

        # Some NN parameters for South
        LR_S=0.007    # learning rate
        NBATCH_S=1024 # Powers of 2
        NCHAIN_S=5    # chains
        NEPOCH_S=100  # number of epochs
        NNS_S=(4 20)  # NN structure (# layers, # units)

        # Find learning rate for North
        $RUN_SYSNET N $TRACER$zbin1 true false $NBATCH_N 0.003 dnnp pnll $DIR $NCHAIN_N $NEPOCH_N ${NNS_N[@]}
        $RUN_SYSNET N $TRACER$zbin2 true false $NBATCH_N 0.003 dnnp pnll $DIR $NCHAIN_N $NEPOCH_N ${NNS_N[@]}
        # Find learning rate for South
        $RUN_SYSNET S $TRACER$zbin1 true false $NBATCH_S 0.003 dnnp pnll $DIR $NCHAIN_S $NEPOCH_S ${NNS_S[@]}
        $RUN_SYSNET S $TRACER$zbin2 true false $NBATCH_S 0.003 dnnp pnll $DIR $NCHAIN_S $NEPOCH_S ${NNS_S[@]}

        # Run SYSNet for North
        srun -N 1 -n 1 $RUN_SYSNET N $TRACER$zbin1 false true $NBATCH_N $LR_N dnnp pnll $DIR $NCHAIN_N $NEPOCH_N ${NNS_N[@]} &
        srun -N 1 -n 1 $RUN_SYSNET N $TRACER$zbin2 false true $NBATCH_N $LR_N dnnp pnll $DIR $NCHAIN_N $NEPOCH_N ${NNS_N[@]} &
        wait

        # Run SYSNet for South
        srun -N 1 -n 1 $RUN_SYSNET S $TRACER$zbin1 false true $NBATCH_S $LR_S dnnp pnll $DIR $NCHAIN_S $NEPOCH_S ${NNS_S[@]} &
        srun -N 1 -n 1 $RUN_SYSNET S $TRACER$zbin2 false true $NBATCH_S $LR_S dnnp pnll $DIR $NCHAIN_S $NEPOCH_S ${NNS_S[@]} &
        wait

    else
        if [ $REGION == 'N' ]
        then
            # Some NN parameters for North
            LR=0.009    # learning rate
            NBATCH=256  # Powers of 2
            NCHAIN=5    # chains
            NEPOCH=100  # number of epochs
            NNS=(3 10)  # NN structure (# layers, # units)
        fi

        if [ $REGION == 'S' ]
        then
            # Some NN parameters for South
            LR=0.007    # learning rate
            NBATCH=1024 # Powers of 2
            NCHAIN=5    # chains
            NEPOCH=100  # number of epochs
            NNS=(4 20)  # NN structure (# layers, # units)
        fi

        $RUN_SYSNET $REGION $TRACER$zbin1 true false $NBATCH 0.003 dnnp pnll $DIR $NCHAIN $NEPOCH ${NNS[@]}
        $RUN_SYSNET $REGION $TRACER$zbin2 true false $NBATCH 0.003 dnnp pnll $DIR $NCHAIN $NEPOCH ${NNS[@]}

        srun -N 1 -n 1 $RUN_SYSNET $REGION $TRACER$zbin1 false true $NBATCH $LR dnnp pnll $DIR $NCHAIN $NEPOCH ${NNS[@]} &
        srun -N 1 -n 1 $RUN_SYSNET $REGION $TRACER$zbin2 false true $NBATCH $LR dnnp pnll $DIR $NCHAIN $NEPOCH ${NNS[@]} &
        wait
    fi
fi


if [ $do_SN_LRG == true ]
then
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
    
    # Find learning rate for North
    $RUN_SYSNET N LRG0.8_1.1 true false $NBATCH_N 0.003 dnnp pnll $DIR $NCHAIN_N $NEPOCH_N ${NNS_N[@]}
    $RUN_SYSNET N LRG0.6_0.8 true false $NBATCH_N 0.003 dnnp pnll $DIR $NCHAIN_N $NEPOCH_N ${NNS_N[@]}
    $RUN_SYSNET N LRG0.4_0.6 true false $NBATCH_N 0.003 dnnp pnll $DIR $NCHAIN_N $NEPOCH_N ${NNS_N[@]}
    # Find learning rate for South
    $RUN_SYSNET S LRG0.8_1.1 true false $NBATCH_S 0.003 dnnp pnll $DIR $NCHAIN_S $NEPOCH_S ${NNS_N[@]}
    $RUN_SYSNET S LRG0.6_0.8 true false $NBATCH_S 0.003 dnnp pnll $DIR $NCHAIN_S $NEPOCH_S ${NNS_N[@]}
    $RUN_SYSNET S LRG0.4_0.6 true false $NBATCH_S 0.003 dnnp pnll $DIR $NCHAIN_S $NEPOCH_S ${NNS_N[@]}

    # Run SYSNet for North
    srun -N 1 $RUN_SYSNET N LRG0.8_1.1 false true $NBATCH_N $LR_N3 dnnp pnll $DIR $NCHAIN_N $NEPOCH_N ${NNS_N[@]} &
    srun -N 1 $RUN_SYSNET N LRG0.6_0.8 false true $NBATCH_N $LR_N2 dnnp pnll $DIR $NCHAIN_N $NEPOCH_N ${NNS_N[@]} &
    srun -N 1 $RUN_SYSNET N LRG0.4_0.6 false true $NBATCH_N $LR_N1 dnnp pnll $DIR $NCHAIN_N $NEPOCH_N ${NNS_N[@]} &
    wait
    # Run SYSNet for South
    srun -N 1 $RUN_SYSNET S LRG0.8_1.1 false true $NBATCH_S $LR_S3 dnnp pnll $DIR $NCHAIN_S $NEPOCH_S ${NNS_S[@]} &
    srun -N 1 $RUN_SYSNET S LRG0.6_0.8 false true $NBATCH_S $LR_S2 dnnp pnll $DIR $NCHAIN_S $NEPOCH_S ${NNS_S[@]} &
    srun -N 1 $RUN_SYSNET S LRG0.4_0.6 false true $NBATCH_S $LR_S1 dnnp pnll $DIR $NCHAIN_S $NEPOCH_S ${NNS_S[@]} &
    wait
fi


if [ $add_SN == true ]
then
    echo "add SYSNet weights to catalogs"
    # add SYSNet weights
    if [ $MOCKCATVER == 'v0' ]
    then
        srun -N 1 -n 1 python $LSSDIR/scripts/mock_tools/addsys.py --realization $REAL --tracer $TRACER  --add_sysnet y --add_sysnet_ran y --par y --base_dir $BASEDIR --use_altmtl $ALTMTL --data_version $DATA_VERSION
    else
        srun -N 1 -n 1 python $LSSDIR/scripts/mock_tools/addsys.py --realization $REAL --tracer $TRACER  --add_sysnet y --add_sysnet_ran y --par y --base_dir $BASEDIR --use_altmtl $ALTMTL --mockcatver $MOCKCATVER --data_version $DATA_VERSION
    fi
    echo "finished adding SYSNet weights"
fi

if [ $do_validation_ffa == true ]
then
    VALIDATION=$LSSDIR/scripts/validation/validation_improp_mock_FFA.py
    echo "running validation on mocks"
    if [ $MOCKCATVER == 'v0' ]
    then
        python $VALIDATION --tracer ELG_LOP --mockn $REAL --dataver $DATA_VERSION --mockversion $MOCKVERSION
        python $VALIDATION --tracer ELG_LOP --mockn $REAL --dataver $DATA_VERSION --mockversion $MOCKVERSION --weight_col WEIGHT_SN 
    else
        python $VALIDATION --tracer ELG_LOP --mockn $REAL --dataver $DATA_VERSION --mockversion $MOCKVERSION --mockcatver $MOCKCATVER
        python $VALIDATION --tracer ELG_LOP --mockn $REAL --dataver $DATA_VERSION --mockversion $MOCKVERSION --mockcatver $MOCKCATVER --weight_col WEIGHT_SN
    fi
    echo "finished validations"
fi


if [ $do_validation_altmtl == true ]
then
    WEIGHT=WEIGHT_SN
    VALIDATION=$LSSDIR/scripts/validation/validation_improp_mock_altmtl.py
    echo "running validation on mocks"
    if [ $MOCKCATVER == 'v0' ]
    then
        srun -N 1 -n 1 python $VALIDATION --tracer $TRACER --mockn $REAL --dataver $DATA_VERSION --mockversion $MOCKVERSION --weight_col $WEIGHT &
        srun -N 1 -n 1 python $VALIDATION --tracer $TRACER --mockn $REAL --dataver $DATA_VERSION --mockversion $MOCKVERSION &
        wait
    else
        srun -N 1 -n 1 python $VALIDATION --tracer $TRACER --mockn $REAL --dataver $DATA_VERSION --mockversion $MOCKVERSION --mockcatver $MOCKCATVER --weight_col $WEIGHT &
        srun -N 1 -n 1 python $VALIDATION --tracer $TRACER --mockn $REAL --dataver $DATA_VERSION --mockversion $MOCKVERSION --mockcatver $MOCKCATVER &
        wait
    fi
    echo "finished validations"
fi

if [ $do_RF_validation_altmtl == true ]
then
    WEIGHT=WEIGHT_RF
    VALIDATION=$LSSDIR/scripts/validation/validation_improp_mock_altmtl.py
    echo "running validation on mocks"
    if [ $MOCKCATVER == 'v0' ]
    then
        srun -N 1 -n 1 python $VALIDATION --tracer $TRACER --mockn $REAL --dataver $DATA_VERSION --mockversion $MOCKVERSION --weight_col $WEIGHT 
        srun -N 1 -n 1 python $VALIDATION --tracer $TRACER --mockn $REAL --dataver $DATA_VERSION --mockversion $MOCKVERSION
    else
        srun -N 1 -n 1 python $VALIDATION --tracer $TRACER --mockn $REAL --dataver $DATA_VERSION --mockversion $MOCKVERSION --mockcatver $MOCKCATVER --weight_col $WEIGHT
        srun -N 1 -n 1 python $VALIDATION --tracer $TRACER --mockn $REAL --dataver $DATA_VERSION --mockversion $MOCKVERSION --mockcatver $MOCKCATVER
    fi
    echo "finished validations"
fi

# Record the end time
end_time=$(date +%s)

# Calculate the difference in seconds
time_diff=$((end_time - start_time))

# Convert seconds to hours, minutes, and seconds
hours=$((time_diff / 3600))
minutes=$(( (time_diff % 3600) / 60 ))
seconds=$((time_diff % 60))

echo "Script execution time: $hours hours $minutes minutes $seconds seconds"