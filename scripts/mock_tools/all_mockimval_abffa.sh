#!/bin/bash
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
PYTHONPATH=$PYTHONPATH:$HOME/LSS/py
for (( i=$1;i<=$2;i++ ))
do
 python scripts/validation/validation_improp_mock_FFA.py --tracer $3 --mockn $i --weight_col $4 --mockcatver v2

done

