#!/usr/bin/env bash

echo "$(date), Step 4-7:concatenate, prepar AltMLT"

#
# script parameters
#
LSS_DIR=$1
DS_DIR=$2           # root directory of mock with version
FIRST_ID_RANK=$3    # first id seed in rank job array

#
# SLURM parameters
#
# task rank
PROCID=${SLURM_PROCID:-0}
# ID seed mock to process in task rank
IDS=$((FIRST_ID_RANK + PROCID ))

NCPU=${SLURM_CPUS_PER_TASK:-1}
# NCPU_M2=$((NCPU-2))

#
# internal variables
#
HOLI_DIR=$LSS_DIR/scripts/mock_tools/holi_pipeline
export PATH=$HOLI_DIR:$PATH
cd $HOLI_DIR
seed=$(printf "seed%04d" "$IDS")
# dataflow name input
in_4=imforFA0_Y3_noimagingmask_applied.fits
in_5=forFA0.fits
in_6=forFA0_withcontaminants.fits
in_7=forFA0_concat.fits
out6=$DS_DIR/$(printf "forFA%04d.fits" "$IDS")


#
# step 4 5 6
#

(
    #
    # step 4 mask
    #
    echo "================= step 4"
    ## Env
    source /global/common/software/desi/users/adematti/cosmodesi_environment.sh dr1
    # use local package LSS, refresh after source env
    export PYTHONPATH=$LSS_DIR/py:$PYTHONPATH
    export PATH=$LSS_DIR/bin:$PATH
    # # login nodes have no GPU; avoid jax/tensorflow CUDA plugin init (silences
    # # duplicate cuDNN/cuFFT/cuBLAS factory registration warnings)
    # export JAX_PLATFORMS=cpu
    # export CUDA_VISIBLE_DEVICES=-1
    # export TF_CPP_MIN_LOG_LEVEL=3
    # # install
    # #python -m pip install -U scalene

    ## in4=out3, out4
    out3_ELG=$DS_DIR/$seed/ELG/$in_4
    out3_LRG=$DS_DIR/$seed/LRG/$in_4
    out3_QSO=$DS_DIR/$seed/QSO/$in_4
    out4_ELG=$DS_DIR/$seed/ELG/$in_5
    out4_LRG=$DS_DIR/$seed/LRG/$in_5
    out4_QSO=$DS_DIR/$seed/QSO/$in_5
    ## join tracers
    time ./join_imaging_mask_stdpars.py --inputs $out3_ELG --outputs $out4_ELG
    
    ## profile version
    # rm $out4_ELG
    # python -m scalene run -o profile1.json ./join_imaging_mask_stdpars.py --inputs $out3_ELG --outputs $out4_ELG
    # python -m scalene view --standalone profile1.json
    # exit 0

    time ./join_imaging_mask_stdpars.py --inputs $out3_LRG --outputs $out4_LRG 
    time ./join_imaging_mask_stdpars.py --inputs $out3_QSO --outputs $out4_QSO 

    # wait

    #
    # step 5 contaminant, **not for LRG**
    #
    echo "================= step 5"
    date
    out5_ELG=$DS_DIR/$seed/ELG/$in_6
    out5_QSO=$DS_DIR/$seed/QSO/$in_6
    time ./add_contaminants_to_mock_stdpars.py --inputs $out4_ELG --outputs $out5_ELG 
    time ./add_contaminants_to_mock_stdpars.py --inputs $out4_QSO --outputs $out5_QSO 
    #wait

    #
    # step 6 concatenate tracers
    #
    echo "================= step 6"
    date
    in6_ELG=$out5_ELG
    in6_QSO=$out5_QSO
    in6_LRG=$out4_LRG
    
    time ./concatenate_tracers_to_fba_stdpars.py --inputs $in6_ELG $in6_LRG $in6_QSO --outputs $out6 --id_seed $IDS
)


#
# step 7 initialize amtl mocks
#
(
    echo "================= step 7"
    date
    ## Env
    source /global/common/software/desi/desi_environment.sh 26.3
    #module swap desitarget/3.0.0
    # local desitarget ecsv to fits
    target_dir=$(get_pars.py $HOLI_PARS TARGET_dir)
    export PYTHONPATH=$target_dir/py:$PYTHONPATH
    export PATH=$target_dir/bin:$PATH
    # use local package LSS, refresh after source env
    export PYTHONPATH=$LSS_DIR/py:$PYTHONPATH
    export PATH=$LSS_DIR/bin:$HOLI_DIR:$PATH

    altmtlxxxx=$(printf "altmtl%04d" "$IDS")
    out7=$DS_DIR/$altmtlxxxx
    echo $out7

    # NOTE: do NOT use $NCPU here. NTASKS seeds already run this step in
    # parallel on the node, so reusing cpus-per-task as the internal
    # make_ledger pool size oversubscribes memory (NTASKS*nproc workers
    # at once) and triggers OOM kills. Use a small, dedicated value.
    NPROC_LEDGER=$(get_pars.py $HOLI_PARS init_amtl.nproc)
    if (( NPROC_LEDGER > NCPU )); then
        NPROC_LEDGER=$NCPU
    fi

    time ./initialize_amtl_mocks_da2_stdpars.py --inputs $out6 --outputs $out7 --obscon DARK --nproc $NPROC_LEDGER
)

date
exit 0
