#!/bin/bash
#SBATCH --nodes=8
##SBATCH --ntasks=8
##SBATCH --ntasks-per-node=1
##SBATCH --cpus-per-task=128
#SBATCH --account=desi
#SBATCH -t 00:30:00
#SBATCH --constraint=cpu
#SBATCH -q debug
#test
scriptdir=/global/homes/d/desica/LSScode/LSS/scripts

for i in $(seq $1 $2); do
    $scriptdir/mock_tools/holi_darkv4_DA3DA2_lss_pipeline_run1_sysnet.sh "$i" &
done

wait

