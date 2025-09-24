#!/bin/bash

echo "Please, choose the tracer: LRG, QSO or ELG"
read -p "Tracer: " TRACER

if [[ "$TRACER" != "LRG" && "$TRACER" != "QSO" && "$TRACER" != "ELG" ]]; then
    echo "Invalid tracer. Please, chose among LRG, QSO or ELG."
    exit 1
fi

source /global/common/software/desi/desi_environment.sh main
VARIABLE="/global/cfs/cdirs/desi/mocks/cai/Uchuu-SHAM/Y3-v1.5/0000/"
VARIABLE2="Uchuu-SHAM_${TRACER}_Y3-v1.5_0000_clustering.dat"

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python prepare_mocks_Y3.py \
    --survey DA2 --tracer "$TRACER" --input_mockpath "$VARIABLE" --input_mockfile "$VARIABLE2" \
    --output_fullpathfn output1.fits --nproc 128



export PATH=/global/cfs/cdirs/desi/users/raichoor/fiberassign-rerun-main/fiberassign_main_godesi23.10/bin:$PATH
export PYTHONPATH=/global/cfs/cdirs/desi/users/raichoor/fiberassign-rerun-main/fiberassign_main_godesi23.10/py:$PYTHONPATH
export SKYHEALPIXS_DIR=$DESI_ROOT/target/skyhealpixs/v1

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python getpota_input_Y3.py \
    --input output1.fits --output output2.fits --realization 1 --prog 'DARK'

if [[ "$TRACER" == "LRG" ]]; then
    source /global/common/software/desi/desi_environment.sh main
    srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python readwrite_pixel_bitmask_ttype_Y3.py \
        -t LRG -i /pscratch/sd/e/efdez/Uchuu/LSS/scripts/mock_tools/output1
    srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python pota2clus_fast_Y3.py \
        --pota_fn output2.fits --tracer LRG --in_tarf output1_LRGimask.fits --mockver 'Uchuu'

elif [[ "$TRACER" == "QSO" ]]; then
    source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
    srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python pota2clus_fast_Y3.py \
        --pota_fn output2.fits --tracer QSO --mockver 'Uchuu' -mk_inputran y

elif [[ "$TRACER" == "ELG" ]]; then
    srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python pota2clus_fast_Y3.py \
        --pota_fn output2.fits --tracer ELG --mockver 'Uchuu'
fi

#optional
rm -f output1.fits output2.fits output2_LRGimask.fits
