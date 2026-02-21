#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --constraint=gpu
#SBATCH --gpus-per-node=4
#SBATCH --array=0-24
#SBATCH --account=desi

source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
PYTHONPATH=$PYTHONPATH:$HOME/LSS

indir=/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/altmtl$SLURM_ARRAY_TASK_ID/mock$SLURM_ARRAY_TASK_ID/LSScats

srun python scripts/xirunpc.py --gpu --tracer LRG --region NGC SGC --corr_type smu --weight_type default_FKP --njack 0 --basedir $indir --outdir $indir/xi/ --rpcut 2.5
srun python scripts/xirunpc.py --gpu --tracer ELG_LOPnotqso --region NGC SGC --corr_type smu --weight_type default_FKP --njack 0 --basedir $indir --outdir $indir/xi/ --rpcut 2.5
srun python scripts/xirunpc.py --gpu --tracer QSO --region NGC SGC --corr_type smu --weight_type default_FKP --njack 0 --basedir $indir --outdir $indir/xi/ --rpcut 2.5