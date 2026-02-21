#!/bin/bash
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
PYTHONPATH=$PYTHONPATH:$HOME/LSS
for (( i=$1;i<=$2;i++ ))
do
 srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/mock_tools/addsys.py --tracer ${3}_ffa --imsys n --add_imsys_ran y --par y --realization $i --mockcatver v2
 python scripts/validation/validation_improp_mock_FFA.py --tracer $3 --mockn $i --weight_col WEIGHT_IMLIN --mockcatver v2

done

