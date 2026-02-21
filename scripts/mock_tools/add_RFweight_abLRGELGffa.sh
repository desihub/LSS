#!/bin/bash
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
PYTHONPATH=$PYTHONPATH:$HOME/LSS
for (( i=$1;i<=$2;i++ ))
do
  srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/mock_tools/addsys.py --tracer LRG_ffa --regressis y --add_regressis y --add_regressis_ran y --par y --realization $i --mockcatver v2
  srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/mock_tools/addsys.py --tracer ELG_LOP_ffa --regressis y --add_regressis y --add_regressis_ran y --par y --realization $i --mockcatver v2

done

