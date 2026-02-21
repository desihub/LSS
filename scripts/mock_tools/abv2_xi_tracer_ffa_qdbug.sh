#!/bin/bash
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
PYTHONPATH=$PYTHONPATH:$HOME/LSS
for (( i=$1;i<=$2;i++ ))
do
 srun -N 1 -C gpu -t 00:20:00 --gpus 4 --qos debug --account desi_g python scripts/xirunpc.py --gpu --tracer ${3}_ffa --region NGC SGC --corr_type smu --weight_type default_FKP --njack 0 --basedir /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/mock$i/v2/ --outdir /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/mock$i/v2/xi/

done

