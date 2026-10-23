#!/bin/bash
source /global/common/software/desi/desi_environment.sh main
module load desitarget/3.0.0
module load LSS/main

for ((i=$1;i<=$2;i++ ))
do
  echo $i
  srun -N 1 -C cpu --cpus-per-task=64 -t 00:40:00 --ntasks=1 --qos interactive --account desi python $LSS/scripts/mock_tools/initialize_amtl_mocks_da2.py /pscratch/sd/d/desica/DA2/mocks/holi_v3/forFA$i.fits /pscratch/sd/d/desica/DA2/mocks/holi_v3/altmtl$i/ DARK
done
	 
