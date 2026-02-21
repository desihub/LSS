#!/bin/bash
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
PYTHONPATH=$PYTHONPATH:$HOME/LSS

ver=v4_1fixran

python scripts/mock_tools/pota2clus_fast.py --realization $1 --mockver AbacusSummit_$ver --outloc prod
