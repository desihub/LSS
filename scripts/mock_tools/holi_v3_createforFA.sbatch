#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=desi
#SBATCH --constraint=cpu
#SBATCH -q regular
#SBATCH -t 00:10:00
#SBATCH --array=450-999

bash scripts/mock_tools/run1_456.sh $SLURM_ARRAY_TASK_ID webjax_v4.81 webjax_v4.80 webjax_v4.80 /pscratch/sd/d/desica/DA2/mocks/holi_v3 holi