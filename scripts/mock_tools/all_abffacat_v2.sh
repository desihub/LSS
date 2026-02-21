#!/bin/bash
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
PYTHONPATH=$PYTHONPATH:$HOME/LSS
for (( i=$1;i<=$2;i++ ))
do
 srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/mock_tools/ffa2clus_fast.py --mockcatver v2  --realization $i

done

