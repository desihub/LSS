#!/bin/bash
#SBATCH --nodes=$1
##SBATCH --ntasks=8
##SBATCH --ntasks-per-node=1
##SBATCH --cpus-per-task=128
#SBATCH --account=desi
#SBATCH -t 00:30:00
#SBATCH --constraint=cpu
#SBATCH -q debug
#test
scriptdir=/global/homes/d/desica/LSScode/LSS/scripts

mmin=$2
mmax=$((mmin+$1-1))

for i in $(seq mmin mmax); do
    srun --exclusive -N 1 -n 1 --cpus-per-task=128 $scriptdir/mock_tools/holi_darkv4_DA3DA2_lss_snprep.sh "$i" &
done

wait

