#!/bin/bash
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
source /global/common/software/desi/desi_environment.sh main

# Change this with your mock name and path!!!
VARIABLE="path/to/mock" #/pscratch/sd/e/efdez/Uchuu/Y1/
VARIABLE2="mock_name" #uchuu-desi-y1_v2_0p65
VARIABLE3="output_file_name" 



############DO NOT CHANGE ANYTHING HERE############
srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python prepare_mocks_Y3_LRG.py --survey DA2 --tracer LRG --input_mockpath "$VARIABLE" --input_mockfile "$VARIABLE2" --output_fullpathfn output1.fits --nproc 128

export PATH=/global/cfs/cdirs/desi/users/raichoor/fiberassign-rerun-main/fiberassign_main_godesi23.10/bin:$PATH
export PYTHONPATH=/global/cfs/cdirs/desi/users/raichoor/fiberassign-rerun-main/fiberassign_main_godesi23.10/py:$PYTHONPATH
export SKYHEALPIXS_DIR=$DESI_ROOT/target/skyhealpixs/v1

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python getpota_input.py --input output1.fits --output output2.fits --realization 1

source /global/common/software/desi/desi_environment.sh main


srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python readwrite_pixel_bitmask_ttype.py -t LRG -i /pscratch/sd/e/efdez/Uchuu/LSS/scripts/mock_tools/output1

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python pota2clus_fast.py --pota_fn output2.fits --tracer LRG --in_tarf output1_LRGimask.fits --mockver 'Uchuu'

rm -f output1.fits output2.fits output2_LRGimask.fits ##you can remove this line if you want to keep any of the output files
####################################################
