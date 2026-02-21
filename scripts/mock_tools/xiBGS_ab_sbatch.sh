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

srun python scripts/xirunpc.py --gpu --tracer BGS_BRIGHT-21.5_ffa --region NGC SGC --corr_type smu --weight_type default_FKP --njack 0 --basedir /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummitBGS/SecondGenMocks/AbacusSummit/mock$SLURM_ARRAY_TASK_ID/ --outdir /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummitBGS/SecondGenMocks/AbacusSummit/mock$SLURM_ARRAY_TASK_ID/xi/