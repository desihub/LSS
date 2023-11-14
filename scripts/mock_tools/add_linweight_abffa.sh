#!/bin/bash
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
PYTHONPATH=$PYTHONPATH:$HOME/LSS/py
for (( i=$1;i<=$2;i++ ))
do
 srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/mock_tools/addsys.py --tracer ${3}_ffa --imsys y --add_imsys_ran y --par y --realization $i
 python scripts/validation/validation_improp_mock.py --tracer $3 --mockn $i --weight_col WEIGHT_IMLIN

done

