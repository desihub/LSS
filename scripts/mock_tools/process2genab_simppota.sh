#!/bin/bash
source /global/common/software/desi/desi_environment.sh main
PYTHONPATH=$PYTHONPATH:$HOME/LSS
for (( i=$1;i<=$2;i++ ))
do
 python scripts/mock_tools/pota2clus_simp.py  --veto _gtlimaging --realization $i --maxr 4
done