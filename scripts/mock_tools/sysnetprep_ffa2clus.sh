#!/bin/bash
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
PYTHONPATH=$PYTHONPATH:$HOME/LSS
for (( i=$1;i<=$2;i++ ))
do
 python scripts/mock_tools/addsys.py --tracer ${3}_ffa --prepsysnet y  --realization $i
done