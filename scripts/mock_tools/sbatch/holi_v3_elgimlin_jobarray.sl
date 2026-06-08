#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=desi
#SBATCH --constraint=cpu
#SBATCH -q regular
#SBATCH -t 00:20:00
#SBATCH --array=0,12,13,14
#this can be used for any jobs that timed out
source /global/common/software/desi/desi_environment.sh main
module load LSS/main
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
mocknum=$SLURM_ARRAY_TASK_ID
scriptdir=/global/homes/d/desica/LSScode/LSS/scripts
sim=holi_v3

srun python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir /global/cfs/cdirs/desi/mocks/cai/LSS/ --outmd cfs --simName $sim --mocknum $mocknum --survey DA2 --specdata loa-v1 --tracer ELG --notqso y  --doimlin y --par y --imsys_zbin split
