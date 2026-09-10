#!/bin/bash
#SBATCH --nodes=300
#SBATCH --account=desi
#SBATCH -t 00:30:00
#SBATCH --constraint=cpu
#SBATCH -q regular
#test
scriptdir=/global/homes/d/desica/LSScode/LSS/scripts


for i in $(seq 0 99); do
    srun --exclusive -N 1 -n 1 --cpus-per-task=128 $scriptdir/mock_tools/holi_bgs_bmask_DA3DA2_lss_pipeline_run1.sh "$i" &
done

for i in $(seq 700 899); do
    srun --exclusive -N 1 -n 1 --cpus-per-task=128 $scriptdir/mock_tools/holi_bgs_bmask_DA3DA2_lss_pipeline_run1.sh "$i" &
done


wait

