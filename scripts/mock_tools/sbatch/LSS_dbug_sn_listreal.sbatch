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

$scriptdir/mock_tools/holi_darkv4_DA3DA2_lss_pipeline_run1_sysnet.sh "1" &
$scriptdir/mock_tools/holi_darkv4_DA3DA2_lss_pipeline_run1_sysnet.sh "324" &
$scriptdir/mock_tools/holi_darkv4_DA3DA2_lss_pipeline_run1_sysnet.sh "328" &
$scriptdir/mock_tools/holi_darkv4_DA3DA2_lss_pipeline_run1_sysnet.sh "778" &

wait

