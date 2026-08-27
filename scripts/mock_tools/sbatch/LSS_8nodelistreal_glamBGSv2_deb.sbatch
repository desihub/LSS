#!/bin/bash
#SBATCH --nodes=8
#SBATCH --account=desi
#SBATCH -t 00:30:00
#SBATCH --constraint=cpu
#SBATCH -q debug
#test
scriptdir=/global/homes/d/desica/LSScode/LSS/scripts


#for i in $(seq $1 $2); do
srun --exclusive -N 1 -n 1 --cpus-per-task=128 $scriptdir/mock_tools/glam_bgsv2_DA3DA2_lss_pipeline_run1.sh "119" &
srun --exclusive -N 1 -n 1 --cpus-per-task=128 $scriptdir/mock_tools/glam_bgsv2_DA3DA2_lss_pipeline_run1.sh "125" &
srun --exclusive -N 1 -n 1 --cpus-per-task=128 $scriptdir/mock_tools/glam_bgsv2_DA3DA2_lss_pipeline_run1.sh "133" &
srun --exclusive -N 1 -n 1 --cpus-per-task=128 $scriptdir/mock_tools/glam_bgsv2_DA3DA2_lss_pipeline_run1.sh "321" &
srun --exclusive -N 1 -n 1 --cpus-per-task=128 $scriptdir/mock_tools/glam_bgsv2_DA3DA2_lss_pipeline_run1.sh "346" &
srun --exclusive -N 1 -n 1 --cpus-per-task=128 $scriptdir/mock_tools/glam_bgsv2_DA3DA2_lss_pipeline_run1.sh "572" &
srun --exclusive -N 1 -n 1 --cpus-per-task=128 $scriptdir/mock_tools/glam_bgsv2_DA3DA2_lss_pipeline_run1.sh "577" &
srun --exclusive -N 1 -n 1 --cpus-per-task=128 $scriptdir/mock_tools/glam_bgsv2_DA3DA2_lss_pipeline_run1.sh "599" &


#done

wait

