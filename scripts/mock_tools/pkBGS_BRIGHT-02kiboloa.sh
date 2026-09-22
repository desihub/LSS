#!/bin/bash
source /global/common/software/desi/desi_environment.sh main

for ((i=$1;i<=$2;i++ ))
do
  echo $i
  srun -N 1 -n 4 clustering-stats --tracer BGS_BRIGHT-02 --zrange 0.1 0.4 --stats mesh2_spectrum --cat_dir /dvs_ro/cfs/cdirs/desi/survey/catalogs//DA2/mocks/SecondGenMocks/AbacusSummitBGS_v2/altmtl$i/kibo-v1/mock$i/LSScats/ --stats_dir $SCRATCH/mock$i/kibo-v1 --combine  --nran 1  
  srun -N 1 -n 4 clustering-stats --tracer BGS_BRIGHT-02 --zrange 0.1 0.4 --stats mesh2_spectrum --cat_dir /dvs_ro/cfs/cdirs/desi/survey/catalogs//DA2/mocks/SecondGenMocks/AbacusSummitBGS_v2/altmtl$i/loa-v1/mock$i/LSScats/ --stats_dir $SCRATCH/mock$i/kibo-v1 --combine  --nran 1 --expand_randoms data-dr2-v2
done