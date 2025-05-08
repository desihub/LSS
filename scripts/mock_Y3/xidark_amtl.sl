#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --qos=debug
#SBATCH --nodes=1
#SBATCH --constraint=gpu
#SBATCH -J xi
#SBATCH --gpus-per-node=4
#SBATCH --array=7-7
#SBATCH --account=desi
#SBATCH -o ./stdout/%x_%A.o%a
#SBATCH -e ./stdout/%x_%A.e%a
##SBATCH --dependency=afterok:37216992

source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
PYTHONPATH=$PYTHONPATH:$HOME/installed_packages/desihub/LSS/py


code_path=$HOME/installed_packages/desihub/LSS/scripts


weight_type=default_FKP
verspec=loa-v1
version=v1.1
indir=/pscratch/sd/j/jerryou/DA2/mocks/EZmock/altmtl$SLURM_ARRAY_TASK_ID/$verspec/mock$SLURM_ARRAY_TASK_ID/LSScats/
outdir=/pscratch/sd/j/jerryou/DA2/mocks/EZmock/xi/$verspec/altmtl$SLURM_ARRAY_TASK_ID/

time srun python ${code_path}/xirunpc.py --gpu --tracer LRG --region NGC SGC --corr_type smu --weight_type ${weight_type} --njack 0 --basedir $indir --outdir $outdir
#time srun python ${code_path}/xirunpc.py --gpu --tracer ELG_LOPnotqso --region NGC SGC --corr_type smu --weight_type ${weight_type} --njack 0 --basedir $indir --outdir $outdir
#time srun python ${code_path}/xirunpc.py --gpu --tracer QSO --region NGC SGC --corr_type smu --weight_type ${weight_type} --njack 0 --basedir $indir --outdir $outdir
