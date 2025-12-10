#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=desi
#SBATCH --constraint=cpu
#SBATCH -q debug
#SBATCH -t 00:30:00

source /global/common/software/desi/desi_environment.sh main
module load LSS/DR2-mocks-v0
mocknum=1
mockver=holi/altMTL
scriptdir=/global/homes/d/desica/LSScode/LSS/scripts


source /global/common/software/desi/users/adematti/cosmodesi_environment.sh test

srun -N 1 -C cpu -t 01:00:00 --qos interactive --account desi python scripts/mock_tools/pota2clus_fast.py --realization $mocknum --base_dir /global/cfs/cdirs/desi/mocks/cai/ --data_dir /global/cfs/cdirs/desi/survey/catalogs/DA2/LSS/loa-v1/LSScats/v2/ --mockver $mockver --mockcatver v0  --specrel loa-v1