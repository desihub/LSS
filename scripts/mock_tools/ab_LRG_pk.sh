#!/bin/bash
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
PYTHONPATH=$PYTHONPATH:$HOME/LSS
for (( i=$1;i<=$2;i++ ))
do
  srun -N 1 -n 64 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/pkrun.py --tracer LRG --region NGC SGC --weight_type default_FKP  --rebinning y --basedir /global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummit_v4_1/altmtl$i/kibo-v1/mock$i/LSScats/ --version DA2_altmtl_$i
done

