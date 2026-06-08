#!/bin/bash
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh test
export OMP_NUM_THREADS=2
PYTHONPATH=$PYTHONPATH:$LSSCODE/LSS/py

tracer=$1
BASEDIR=$2

sysnet_app=$LSSCODE/LSS/py/LSS/imaging/sysnet_app_mpi.py
sysnet_dir=$BASEDIR/sysnet/

echo $sysnet_dir $tracer

if [[ "$tracer" == *"ELG"* ]]; then
    regions=("N" "S")
    zlims=(0.8 1.1 1.1 1.6) # each pair is a redshift range
    # Some NN parameters for North
    LR=(0.009 0.007)    # learning rate
    NBATCH=(256 1024)   # batchsize, powers of 2
    NNS=(3 10 4 20)     # each pair is a NN structure (# layers, # units)
    NCHAIN=5            # number of chains
    NEPOCH=100          # number of epochs
    
elif [[ "$tracer" == *"LRG"* ]]; then
    regions=("N" "S")
    zlims=(0.4 0.6 0.6 0.8 0.8 1.1)
    LR=(0.009 0.007)
    NBATCH=(256 1024)
    NNS=(3 10 4 20)
    NCHAIN=5
    NEPOCH=100
    
elif [[ "$tracer" == *"QSO"* ]]; then
    regions=("N" "SnotDES" "DES")
    zlims=(0.8 1.3 1.3 2.1 2.1 3.5)
    LR=(0.01 0.01 0.02)
    NBATCH=(256 1024 256)
    NNS=(3 10 4 20 3 10)
    NCHAIN=5
    NEPOCH=100
    
else
    echo "ERROR: Tracer $tracer not supported"
    exit 1
fi

echo ${regions[@]}
lrfinder_flags="-ax all --model dnnp --loss pnll -fl"
train_flags="-ax all --model dnnp --loss pnll --eta_min 0.00001 -k"

for ((i=0; i<${#regions[@]}; i++)); do
    #iterate over regions
    NNlayers="${NNS[$((i*2))]}"
    NNunits="${NNS[$((i*2+1))]}"
    for ((j=0; j<$(( ${#zlims[@]} / 2 )); j++)); do
        zmin="${zlims[$((j*2))]}"
        zmax="${zlims[$((j*2 + 1))]}"
        region_flags="-lr ${LR[$i]} -bs ${NBATCH[$i]} --nn_structure $NNlayers $NNunits -ne $NEPOCH -nc $NCHAIN -i $sysnet_dir/prep_${tracer}${zmin}_${zmax}_${regions[$i]}.fits -o $sysnet_dir/${tracer}${zmin}_${zmax}_${regions[$i]}"
        echo $region_flags
        # find learning rate 
        srun -n 1  -t 5  $srun_flags python $sysnet_app $lrfinder_flags $region_flags
        # start training
        srun -n 25 -t 10 $srun_flags python $sysnet_app $train_flags $region_flags
    done
done