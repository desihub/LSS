#!/bin/bash
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
PYTHONPATH=$PYTHONPATH:$HOME/LSS
for (( i=$1;i<=$2;i++ ))
do
  python scripts/mock_tools/addsys.py --tracer BGS_BRIGHT-21.5_ffa --imsys y --add_imsys_ran y --par y --realization $i --base_dir /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummitBGS/FFA/ 

done

