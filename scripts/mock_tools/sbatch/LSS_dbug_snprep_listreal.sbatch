#!/bin/bash
#SBATCH --nodes=4
##SBATCH --ntasks=4
##SBATCH --ntasks-per-node=1
##SBATCH --cpus-per-task=128
#SBATCH --account=desi
#SBATCH -t 00:30:00
#SBATCH --constraint=cpu
#SBATCH -q debug
#test
scriptdir=/global/homes/d/desica/LSScode/LSS/scripts

srun --exclusive -N 1 -n 1 --cpus-per-task=128 $scriptdir/mock_tools/holi_darkv4_DA3DA2_lss_snprep.sh "1" &
srun --exclusive -N 1 -n 1 --cpus-per-task=128 $scriptdir/mock_tools/holi_darkv4_DA3DA2_lss_snprep.sh "324" &
srun --exclusive -N 1 -n 1 --cpus-per-task=128 $scriptdir/mock_tools/holi_darkv4_DA3DA2_lss_snprep.sh "328" &
srun --exclusive -N 1 -n 1 --cpus-per-task=128 $scriptdir/mock_tools/holi_darkv4_DA3DA2_lss_snprep.sh "778" &

wait

