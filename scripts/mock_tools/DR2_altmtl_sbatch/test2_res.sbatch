#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=desi
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=256
#SBATCH --constraint=cpu
##SBATCH -q debug
#SBATCH -q regular
#SBATCH -t 72:00:00
#SBATCH -o /pscratch/sd/d/desica/DA3/mocks/holi_v4/altmtl/%x.o%j
#SBATCH -e /pscratch/sd/d/desica/DA3/mocks/holi_v4/altmtl/%x.e%j
##SBATCH --dependency=afterany:46276251
#SBATCH --nodes=1
#SBATCH --reservation=_CAP_holi_dr3_altmtl

echo "about to load desi env"
source /global/common/software/desi/desi_environment.sh main
echo "loaded desi env"