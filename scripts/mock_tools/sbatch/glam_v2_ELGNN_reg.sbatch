#!/bin/bash
#SBATCH --nodes=41
#SBATCH --account=desi
#SBATCH -t 00:30:00
#SBATCH --constraint=cpu
#SBATCH -q regular
#test
scriptdir=/global/homes/d/desica/LSScode/LSS/scripts


for i in $(seq $1 $2); do
    $scriptdir/mock_tools/glam_lss_elgNN_pipeline_interactive.sh "$i" &
done

wait

