#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --constraint=cpu
#SBATCH --array=2,3,11,32,36,40,43,68,72,85,98,104,108,113,114,217,228,279,282,317,320,385,416,630,649,683,692,726,753,754,926,927,928,929,937,938,939,940,941,973,975,978,979,986,988,989,990,991,992,993,994,995,996,997,998,999,330,487,706
#SBATCH --account=desi

source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
PYTHONPATH=$PYTHONPATH:$HOME/LSS/py

srun scripts/mock_tools/run1_EZmock_LSS.sh $SLURM_ARRAY_TASK_ID