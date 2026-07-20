#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --account=desi
#SBATCH --constraint=cpu
#SBATCH -q regular
#SBATCH -J Holi8
#SBATCH -t 36:00:00
#SBATCH --output=holi8_%j.out
#SBATCH --error=holi8_%j.err
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=jcolley@lpnhe.in2p3.fr

#
# script parameters
#
LSS_DIR=$1
DS_DIR=$2      # root directory of mock with version
FIRST_ID=$3    # id seed to process


cd $LSS_DIR/scripts/mock_tools/holi_pipeline

srun --exclusive -n $SLURM_NTASKS -c $SLURM_CPUS_PER_TASK ./step8.sh $LSS_DIR $DS_DIR $FIRST_ID