#!/bin/bash
#SBATCH -N 2
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -J sys_Y1ab2ndgen
#SBATCH -t 4:00:00
#SBATCH -L SCRATCH
#SBATCH --output=/path/to/slurm_outputs/Y1ab2ndgen_sysnet_%A_%a.out
#SBATCH --array=0-24

# This script can be run as a sbatch job, or an interactive job
# if run in an interactive node then manually set mock realization 
echo $SLURM_ARRAY_TASK_ID

source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main 
export LSSDIR=$HOME/LSS
PYTHONPATH=$PYTHONPATH:$LSSDIR/py

RUN_SYSNET=$LSSDIR/scripts/run_sysnetELG_ab2ndgen.sh 

# REGION only useful for ELGs right now
REGION='all'
TRACER=ELG_LOP_ffa
MOCKVERSION='v2'
REAL=$SLURM_ARRAY_TASK_ID
do_RF=false
do_prep=true
do_SN=true
add_SN=true # if true make sure all weights are available, if not this will fail
do_validation=true

SUBSTRING=${TRACER:0:3}
echo $SUBSTRING

if [ $do_RF == true ]
then
    echo "doing RF regression and adding weights"
    # Compute and add regressis weights to mocks
    srun -N 1 -n 1 python $LSSDIR/scripts/mock_tools/addsys.py --realization $REAL --tracer $TRACER --regressis y --add_regressis y --add_regressis_ran y --par y --mock_version $MOCKVERSION
    echo "finished doing RF regression and adding RF weights"
fi

if [ $do_prep == true ]
then
    echo "preparing data for SYSNet"
    # prepare data for SYSNet
    srun -N 1 -n 1 python $LSSDIR/scripts/mock_tools/addsys.py --tracer $TRACER --prepsysnet y --realization $REAL --mockcatver $MOCKVERSION --use_allsky_rands 'y'
fi

if [ $do_SN == true ]
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
        NEPOCH_N=200  # number of epochs
        NNS_N=(3 10)  # NN structure (# layers, # units)

        # Some NN parameters for South
        LR_S=0.007    # learning rate
        NBATCH_S=1024 # Powers of 2
        NCHAIN_S=5    # chains
        NEPOCH_S=200  # number of epochs
        NNS_S=(4 20)  # NN structure (# layers, # units)

        # Find learning rate for North
        $RUN_SYSNET N $TRACER$zbin1 true false $NBATCH_N 0.003 dnnp pnll $REAL $NCHAIN_N $NEPOCH_N ${NNS_N[@]} $MOCKVERSION
        $RUN_SYSNET N $TRACER$zbin2 true false $NBATCH_N 0.003 dnnp pnll $REAL $NCHAIN_N $NEPOCH_N ${NNS_N[@]} $MOCKVERSION
        # Find learning rate for South
        $RUN_SYSNET S $TRACER$zbin1 true false $NBATCH_S 0.003 dnnp pnll $REAL $NCHAIN_S $NEPOCH_S ${NNS_S[@]} $MOCKVERSION
        $RUN_SYSNET S $TRACER$zbin2 true false $NBATCH_S 0.003 dnnp pnll $REAL $NCHAIN_S $NEPOCH_S ${NNS_S[@]} $MOCKVERSION

        # Run SYSNet for North
        srun -N 1 -n 1 $RUN_SYSNET N $TRACER$zbin1 false true $NBATCH_N $LR_N dnnp pnll $REAL $NCHAIN_N $NEPOCH_N ${NNS_N[@]} $MOCKVERSION &
        srun -N 1 -n 1 $RUN_SYSNET N $TRACER$zbin2 false true $NBATCH_N $LR_N dnnp pnll $REAL $NCHAIN_N $NEPOCH_N ${NNS_N[@]} $MOCKVERSION &
        wait

        # Run SYSNet for South
        srun -N 1 -n 1 $RUN_SYSNET S $TRACER$zbin1 false true $NBATCH_S $LR_S dnnp pnll $REAL $NCHAIN_S $NEPOCH_S ${NNS_S[@]} $MOCKVERSION &
        srun -N 1 -n 1 $RUN_SYSNET S $TRACER$zbin2 false true $NBATCH_S $LR_S dnnp pnll $REAL $NCHAIN_S $NEPOCH_S ${NNS_S[@]} $MOCKVERSION &
        wait

    else
        if [ $REGION == 'N' ]
        then
            # Some NN parameters for North
            LR=0.009    # learning rate
            NBATCH=256  # Powers of 2
            NCHAIN=5    # chains
            NEPOCH=200  # number of epochs
            NNS=(3 10)  # NN structure (# layers, # units)
        fi

        if [ $REGION == 'S' ]
        then
            # Some NN parameters for South
            LR=0.007    # learning rate
            NBATCH=1024 # Powers of 2
            NCHAIN=5    # chains
            NEPOCH=200  # number of epochs
            NNS=(4 20)  # NN structure (# layers, # units)
        fi

        $RUN_SYSNET $REGION $TRACER$zbin1 true false $NBATCH 0.003 dnnp pnll $REAL $NCHAIN $NEPOCH ${NNS[@]} $MOCKVERSION
        $RUN_SYSNET $REGION $TRACER$zbin2 true false $NBATCH 0.003 dnnp pnll $REAL $NCHAIN $NEPOCH ${NNS[@]} $MOCKVERSION

        srun -N 1 -n 1 $RUN_SYSNET $REGION $TRACER$zbin1 false true $NBATCH $LR dnnp pnll $REAL $NCHAIN $NEPOCH ${NNS[@]} $MOCKVERSION &
        srun -N 1 -n 1 $RUN_SYSNET $REGION $TRACER$zbin2 false true $NBATCH $LR dnnp pnll $REAL $NCHAIN $NEPOCH ${NNS[@]} $MOCKVERSION &
        wait
    fi
fi

if [ $add_SN == true ]
then
    echo "add SYSNet weights to catalogs"
    # add SYSNet weights
    srun -N 1 -n 1 python $LSSDIR/scripts/mock_tools/addsys.py --realization $REAL --tracer $TRACER --add_sysnet y --add_sysnet_ran y --par y --mockcatver $MOCKVERSION
    echo "finished adding SYSNet weights"
fi

if [ $do_validation == true ]
then
    echo "running validation on mocks"
    python $LSSDIR/scripts/validation/validation_improp_mock_FFA.py --tracer ELG_LOP --mockn $REAL --mockcatver $MOCKVERSION
    python $LSSDIR/scripts/validation/validation_improp_mock_FFA.py --tracer ELG_LOP --mockn $REAL --mockcatver $MOCKVERSION --weight_col WEIGHT_SN
    echo "finished validations"
fi