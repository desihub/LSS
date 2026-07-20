#!/usr/bin/env bash

# Description:
# ===========
# top level of Holi pipeline with slurm
# job array to process chunk of simulation = cpus-per-task

#SBATCH --ntasks=10
#SBATCH --cpus-per-task=20
#SBATCH --account=desi
#SBATCH --constraint=cpu
#SBATCH -q regular
#SBATCH -J Holi1-7
#SBATCH -t 9:10:00
#SBATCH --array=0-1
#SBATCH --output=holi1-7_%j.out
#SBATCH --error=holi1-7_%j.err
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=jcolley@lpnhe.in2p3.fr


FIRST_ID=$1
#
# Environment
#
#source /global/common/software/desi/desi_environment.sh main
#source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main

#TODO: use file parameters
LSS_DIR=/global/cfs/cdirs/desi/users/colley/LSS
#DS_DIR=/global/cfs/cdirs/desi/mocks/fa4acm/holi/webjax_v4.80
DS_DIR=/pscratch/sd/j/jcolley/holi/webjax_v4.80
# 
HOLI_DIR=$LSS_DIR/scripts/mock_tools/holi_pipeline
cd $HOLI_DIR


# run chunk
# chunk size is the number of CPU : cpus-per-task
./slurm2_chunk_Ncpu.sh $LSS_DIR $DS_DIR $FIRST_ID
