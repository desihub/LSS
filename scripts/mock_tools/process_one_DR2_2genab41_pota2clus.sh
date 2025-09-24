#!/bin/bash
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
export LSSCODE=$HOME/LSScode
PYTHONPATH=$PYTHONPATH:$LSSCODE/LSS/py

python scripts/mock_tools/pota2clus_fast.py --realization $1 --outloc prod
