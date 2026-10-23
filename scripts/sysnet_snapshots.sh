#!/bin/bash
# This script generates and ensemble of SYSNet models from SYSNet snapshots.
# 1. First train the NN model, saving snapshots every T_0 epochs.
# 2. Pull the snapshots, forward model and save predictions in HEALPix maps with Nside=256.
# Note: This script automatically runs in an interactive node
# export LSSCODE=$HOME # Change to your LSS code directory.
# bash $LSSCODE/LSS/scripts/sysnet_snapshots.sh ELG
# bash $LSSCODE/LSS/scripts/sysnet_snapshots.sh ELG_LOP
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
export OMP_NUM_THREADS=2
PYTHONPATH=$PYTHONPATH:$LSSCODE/LSS/py
# export PYTHONPATH=$HOME/sysnetdev:$PYTHONPATH

# /global/cfs/cdirs/desi/survey/catalogs/DA2/LSS/loa-v1/LSScats/v2/
survey=DA2
verspec=loa-v1
version=v2
type=$1

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
NCHAIN_N=25    # chains
NEPOCH_N=400  # number of epochs
NNS_N=(3 10)  # NN structure (# layers, # units)

# Some NN parameters for South
LR_S=0.007    # learning rate
NBATCH_S=1024 # Powers of 2
NCHAIN_S=25    # chains
NEPOCH_S=400  # number of epochs
NNS_S=(4 20)  # NN structure (# layers, # units)

LSSBASE=/global/cfs/cdirs/desi/survey/catalogs/
BASEDIR=$LSSBASE/$survey/LSS/$verspec/LSScats/
echo $BASEDIR

# Flags and directories used by SYSNet 
sysnet_app=$LSSCODE/LSS/py/LSS/imaging/sysnet_appensemble_mpi.py
sysnet_dir=$BASEDIR/$version/sysnet/ # This directory should contain the prepared data tables
sysnet_snap_dir=$BASEDIR/$version/sysnet_snapshots/
# sysnet_snap_dir=$SCRATCH/$version/sysnet_snapshots/ # output directory
srun_flags="-N 1 -C cpu --qos interactive --account desi "
train_flags="-ax all --model dnnp --loss pnll --eta_min 0.00001 --snapshot_ensemble -k --no_eval"
north_flags="-lr $LR_N -bs $NBATCH_N --nn_structure ${NNS_N[@]} -ne $NEPOCH_N -nc $NCHAIN_N"
south_flags="-lr $LR_S -bs $NBATCH_S --nn_structure ${NNS_S[@]} -ne $NEPOCH_S -nc $NCHAIN_S"

# Train NN
    # for North
srun -n 125 -t 30 $srun_flags python $sysnet_app $train_flags $north_flags -i $sysnet_dir/prep_${tracer}0.8_1.1_N.fits -o $sysnet_snap_dir/${tracer}0.8_1.1_N
srun -n 125 -t 30 $srun_flags python $sysnet_app $train_flags $north_flags -i $sysnet_dir/prep_${tracer}1.1_1.6_N.fits -o $sysnet_snap_dir/${tracer}1.1_1.6_N
    # for South
srun -n 125 -t 30 $srun_flags python $sysnet_app $train_flags $south_flags -i $sysnet_dir/prep_${tracer}0.8_1.1_S.fits -o $sysnet_snap_dir/${tracer}0.8_1.1_S
srun -n 125 -t 30 $srun_flags python $sysnet_app $train_flags $south_flags -i $sysnet_dir/prep_${tracer}1.1_1.6_S.fits -o $sysnet_snap_dir/${tracer}1.1_1.6_S

# Pull and combine SYSNet snapshots
pull_snaps=$LSSCODE/LSS/py/LSS/imaging/pull_sysnet_snapshot_and_combine_mpi.py
snap_flags="--type $tracer --north_nn_structure ${NNS_N[@]} --south_nn_structure ${NNS_S[@]} --model dnnp"
srun -n 125 -t 30 $srun_flags python $pull_snaps $snap_flags -i $sysnet_dir -sd $sysnet_snap_dir -o $sysnet_snap_dir