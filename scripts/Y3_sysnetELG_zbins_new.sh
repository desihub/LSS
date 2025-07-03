#!/bin/bash
# To run this script in an interactive node:
# srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi Y3_sysnetELG_zbins_new.sh test
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh test

export LSSDIR=$HOME
export LSSBASE=/global/cfs/cdirs/desi/survey/catalogs/
PYTHONPATH=$PYTHONPATH:$LSSDIR/LSS/py

version=$1
verspec=loa-v1
survey=DA2

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

BASEDIR=$LSSBASE/$survey/LSS/$verspec/LSScats/
RUN_SYSNET=$LSSDIR/LSS/scripts/run_sysnetELG_cd_mpi.sh

# If using the allsky randoms option when preparing for sysnet mkCat_main.py might not work without salloc or job
python scripts/main/mkCat_main.py --basedir $LSSBASE --type ELG_LOP --notqso y --prepsysnet y --imsys_zbin y --fulld n --survey $survey --verspec $verspec --version $version --use_allsky_rands y

# Find learning rate for North
$RUN_SYSNET N ELG_LOPnotqso0.8_1.1 true false $NBATCH_N 0.003 dnnp pnll $VERSION $BASEDIR $NCHAIN_N $NEPOCH_N ${NNS_N[@]}
$RUN_SYSNET N ELG_LOPnotqso1.1_1.6 true false $NBATCH_N 0.003 dnnp pnll $VERSION $BASEDIR $NCHAIN_N $NEPOCH_N ${NNS_N[@]}
# Find learning rate for South
$RUN_SYSNET S ELG_LOPnotqso0.8_1.1 true false $NBATCH_S 0.003 dnnp pnll $VERSION $BASEDIR $NCHAIN_S $NEPOCH_S ${NNS_S[@]}
$RUN_SYSNET S ELG_LOPnotqso1.1_1.6 true false $NBATCH_S 0.003 dnnp pnll $VERSION $BASEDIR $NCHAIN_S $NEPOCH_S ${NNS_S[@]}

# Run SYSNet for North
$RUN_SYSNET N ELG_LOPnotqso0.8_1.1 false true $NBATCH_N $LR_N dnnp pnll $VERSION $BASEDIR $NCHAIN_N $NEPOCH_N ${NNS_N[@]}
$RUN_SYSNET N ELG_LOPnotqso1.1_1.6 false true $NBATCH_N $LR_N dnnp pnll $VERSION $BASEDIR $NCHAIN_N $NEPOCH_N ${NNS_N[@]} 
# Run SYSNet for South
$RUN_SYSNET S ELG_LOPnotqso0.8_1.1 false true $NBATCH_S $LR_S dnnp pnll $VERSION $BASEDIR $NCHAIN_S $NEPOCH_S ${NNS_S[@]}
$RUN_SYSNET S ELG_LOPnotqso1.1_1.6 false true $NBATCH_S $LR_S dnnp pnll $VERSION $BASEDIR $NCHAIN_S $NEPOCH_S ${NNS_S[@]}

python scripts/main/mkCat_main.py --basedir $LSSBASE --type ELG_LOP --notqso y --add_sysnet y --imsys_zbin y --fulld n --survey $survey --verspec $verspec --version $version

python scripts/validation/validation_improp_full.py --tracers ELG_LOPnotqso --version $version --verspec $verspec --survey $survey --weight_col WEIGHT_SN
