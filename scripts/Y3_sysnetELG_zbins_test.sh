#!/bin/bash
# bash $LSSCODE/LSS/scripts/Y3_sysnetELG_zbins_test.sh v2 ELG
# bash $LSSCODE/LSS/scripts/Y3_sysnetELG_zbins_test.sh v2 ELG_LOP
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh test
export OMP_NUM_THREADS=2
PYTHONPATH=$PYTHONPATH:$LSSCODE/LSS/py

version=$1
verspec=loa-v1
survey=DA2
type=$2

if [ $type = 'ELG_LOP' ] 
then
    tracer=ELG_LOPnotqso
fi

if [ $type = 'ELG' ] 
then
    tracer=ELGnotqso
fi

echo $tracer

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

LSSBASE=/global/cfs/cdirs/desi/survey/catalogs/
BASEDIR=$LSSBASE/$survey/LSS/$verspec/LSScats/
echo $BASEDIR

# Flags and directories used by SYSNet 
sysnet_app=$LSSCODE/LSS/py/LSS/imaging/sysnet_app_mpi.py
sysnet_dir=$BASEDIR/$version/sysnet/
srun_flags="-N 1 -C cpu --qos interactive --account desi "
lrfinder_flags="-ax all --model dnnp --loss pnll -fl"
train_flags="-ax all --model dnnp --loss pnll --eta_min 0.00001 -k"
north_flags="-lr $LR_N -bs $NBATCH_N --nn_structure ${NNS_N[@]} -ne $NEPOCH_N -nc $NCHAIN_N"
south_flags="-lr $LR_S -bs $NBATCH_S --nn_structure ${NNS_S[@]} -ne $NEPOCH_S -nc $NCHAIN_S"

# If using the allsky randoms option when preparing for sysnet mkCat_main.py might not work without salloc or job
srun -n 1 -t 20 $srun_flags python scripts/main/mkCat_main.py --basedir $LSSBASE --type $type --notqso y --prepsysnet y --imsys_zbin y --fulld n --survey $survey --verspec $verspec --version $version --use_allsky_rands y

# Get learning rates
    # for North
srun -n 1 -t 5 $srun_flags python $sysnet_app $lrfinder_flags $north_flags -i $sysnet_dir/prep_${tracer}0.8_1.1_N.fits -o $sysnet_dir/${tracer}0.8_1.1_N
srun -n 1 -t 5 $srun_flags python $sysnet_app $lrfinder_flags $north_flags -i $sysnet_dir/prep_${tracer}1.1_1.6_N.fits -o $sysnet_dir${tracer}1.1_1.6_N
    # for South
srun -n 1 -t 5 $srun_flags python $sysnet_app $lrfinder_flags $south_flags -i $sysnet_dir/prep_${tracer}0.8_1.1_S.fits -o $sysnet_dir/${tracer}0.8_1.1_S
srun -n 1 -t 5 $srun_flags python $sysnet_app $lrfinder_flags $south_flags -i $sysnet_dir/prep_${tracer}1.1_1.6_S.fits -o $sysnet_dir/${tracer}1.1_1.6_S

# Train NN
    # for North
srun -n 25 -t 10 $srun_flags python $sysnet_app $train_flags $north_flags -i $sysnet_dir/prep_${tracer}0.8_1.1_N.fits -o $sysnet_dir/${tracer}0.8_1.1_N
srun -n 25 -t 10 $srun_flags python $sysnet_app $train_flags $north_flags -i $sysnet_dir/prep_${tracer}1.1_1.6_N.fits -o $sysnet_dir/${tracer}1.1_1.6_N
    # for South
srun -n 25 -t 10 $srun_flags python $sysnet_app $train_flags $south_flags -i $sysnet_dir/prep_${tracer}0.8_1.1_S.fits -o $sysnet_dir/${tracer}0.8_1.1_S
srun -n 25 -t 10 $srun_flags python $sysnet_app $train_flags $south_flags -i $sysnet_dir/prep_${tracer}1.1_1.6_S.fits -o $sysnet_dir/${tracer}1.1_1.6_S

srun -n 1 -t 20 $srun_flags python scripts/main/mkCat_main.py --basedir $LSSBASE --type $type --notqso y --add_sysnet y --imsys_zbin y --fulld n --survey $survey --verspec $verspec --version $version

#python scripts/validation/validation_improp_full.py --tracers $type --version $version --verspec $verspec --survey $survey --weight_col WEIGHT_SN
