#!/bin/bash
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
PYTHONPATH=$PYTHONPATH:$HOME/LSS
for ((i=$1;i<=$2;i++ ))
do
 srun -N 1 -n 128 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/pkrun.py --tracer LRG --survey Y1 --verspec iron --version mock$i --basedir /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_v4_2/altmtl$i/mock$i/LSScats/  --region NGC SGC --weight_type default_FKP --rebinning y --zlim 0.4 1.1 --option noNorth
 srun -N 1 -n 128 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/pkrun.py --tracer LRG --survey Y1 --verspec iron --version mock$i --basedir /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_v4_2/altmtl$i/mock$i/LSScats/  --region NGC SGC --weight_type default_FKP --rebinning y --zlim 0.4 1.1 --option noDES
 srun -N 1 -n 128 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/pkrun.py --tracer LRG --survey Y1 --verspec iron --version mock$i --basedir /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_v4_2/altmtl$i/mock$i/LSScats/  --region NGC SGC --weight_type default_FKP --rebinning y --zlim 0.4 1.1 
done

