#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=desi
#SBATCH --constraint=cpu
#SBATCH -q debug
#SBATCH -t 00:10:00
#SBATCH -J kibo_comb

source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main

BASE_DIR="/dvs_ro/cfs/cdirs/desi/survey/catalogs/DA2/LSS/kibo-v1/LSScats/v1/BAO/blinded/"
SAVE_DIR="/pscratch/sd/s/sandersn/DA2/kibo-v1/v1/blinded/"
CAP="SGC"
VERBOSE="True"

srun python combined_catalog.py --base_dir $BASE_DIR --save_dir $SAVE_DIR --cap $CAP --verbose $VERBOSE
