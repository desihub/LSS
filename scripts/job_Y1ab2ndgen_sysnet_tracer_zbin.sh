#!/bin/bash
#SBATCH -N 2
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -J sys_Y1ab2ndgen
#SBATCH -t 4:00:00
#SBATCH -L SCRATCH
#SBATCH --output=/path/to/slurm_outputs/Y1ab2ndgen_sysnet_%A_%a.out
#SBATCH --array=0-24

set -e
# This script can be run as a sbatch job, or an interactive job
# if run in an interactive node then manually set mock realization 
echo $SLURM_ARRAY_TASK_ID

source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main 
export LSSDIR=$HOME/LSS
PYTHONPATH=$PYTHONPATH:$LSSDIR/py

RUN_SYSNET=$LSSDIR/scripts/run_sysnetELG_ab2ndgen.sh 

# REGION only useful for ELGs right now
DATA_VERSION='v1.3'
REGION='all'
TRACER=ELG_LOP_ffa
BASEDIR='/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_v4_1/'
MOCKVERSION='v0'
#REAL=$SLURM_ARRAY_TASK_ID
REAL=0
do_RF=false
do_prep=false
do_SN=false
add_SN=true # if true make sure all weights are available, if not this will fail
do_validation=false

SUBSTRING=${TRACER:0:3}
echo $SUBSTRING

if [ $do_RF == true ]
then
    echo "doing RF regression and adding weights"
    # Compute and add regressis weights to mocks
    srun -N 1 -n 1 python $LSSDIR/scripts/mock_tools/addsys.py --realization $REAL --tracer $TRACER --regressis y --add_regressis y --add_regressis_ran y --par y --base_dir $BASEDIR --mock_version $MOCKVERSION --data_version $DATA_VERSION
    echo "finished doing RF regression and adding RF weights"
fi

if [ $do_prep == true ]
then
    echo "preparing data for SYSNet"
    # prepare data for SYSNet
    if [ $MOCKVERSION == 'v0' ]
    then
        srun -N 1 -n 1 python $LSSDIR/scripts/mock_tools/addsys.py --tracer $TRACER --prepsysnet y --base_dir $BASEDIR --realization $REAL --use_allsky_rands 'y' --data_version $DATA_VERSION
    else
        srun -N 1 -n 1 python $LSSDIR/scripts/mock_tools/addsys.py --tracer $TRACER --prepsysnet y --base_dir $BASEDIR --realization $REAL --mockcatver $MOCKVERSION --use_allsky_rands 'y' --data_version $DATA_VERSION
    fi
fi

if [ $do_SN == true ]
then
    # Run SYSNET different options for different zbins
    # Must have data prepared
    echo "running SYSNet on ${TRACER}"
    zbin1=0.8_1.1
    zbin2=1.1_1.6
    DIR=$BASEDIR/mock$REAL/
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

if [ $add_SN == true ]
then
    echo "add SYSNet weights to catalogs"
    # add SYSNet weights
    if [ $MOCKVERSION == 'v0' ]
    then
        srun -N 1 -n 1 python $LSSDIR/scripts/mock_tools/addsys.py --realization $REAL --tracer $TRACER  --add_sysnet y --add_sysnet_ran y --par y --base_dir $BASEDIR --data_version $DATA_VERSION
    else
        srun -N 1 -n 1 python $LSSDIR/scripts/mock_tools/addsys.py --realization $REAL --tracer $TRACER  --add_sysnet y --add_sysnet_ran y --par y --base_dir $BASEDIR --mockcatver $MOCKVERSION --data_version $DATA_VERSION
    fi
    echo "finished adding SYSNet weights"
fi

if [ $do_validation == true ]
then
    echo "running validation on mocks"
    if [ $MOCKVERSION == 'v0' ]
    then
        python $LSSDIR/scripts/validation/validation_improp_mock_FFA.py --tracer ELG_LOP --mockn $REAL --base_dir $BASEDIR --data_version $DATA_VERSION
        python $LSSDIR/scripts/validation/validation_improp_mock_FFA.py --tracer ELG_LOP --mockn $REAL --base_dir $BASEDIR --weight_col WEIGHT_SN --data_version $DATA_VERSION
    else
        python $LSSDIR/scripts/validation/validation_improp_mock_FFA.py --tracer ELG_LOP --mockn $REAL --base_dir $BASEDIR --mockcatver $MOCKVERSION --data_version $DATA_VERSION
        python $LSSDIR/scripts/validation/validation_improp_mock_FFA.py --tracer ELG_LOP --mockn $REAL --base_dir $BASEDIR --mockcatver $MOCKVERSION --weight_col WEIGHT_SN --data_version $DATA_VERSION
    fi
    echo "finished validations"
fi