#!/bin/bash
#SBATCH --nodes=300
#SBATCH --account=desi
#SBATCH -t 00:30:00
#SBATCH --constraint=cpu
#SBATCH -q regular
#test
scriptdir=/global/homes/d/desica/LSScode/LSS/scripts


for i in $(seq $1 $2); do
    $scriptdir/mock_tools/holiv4_lss_elgNN_pipeline.sh "$i" &
done

wait

