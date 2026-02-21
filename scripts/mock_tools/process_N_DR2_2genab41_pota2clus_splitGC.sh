#!/bin/bash
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
export LSSCODE=$HOME/LSScode
PYTHONPATH=$PYTHONPATH:$LSSCODE/LSS

for (( i=$1;i<=$2;i++ ))
do
 python scripts/mock_tools/pota2clus_fast.py --realization $i --outloc prod --mkdat n --mkran n --nz n 
done
