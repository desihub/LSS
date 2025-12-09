#!/bin/bash

#SBATCH -C cpu
#SBATCH -N 1
#SBATCH --qos regular
#SBATCH --mail-user=lucasnap1@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --account desi
#SBATCH --output /pscratch/sd/d/desica/AltMTL/runaltmtldark.log
#SBATCH --time=47:00:00

source /global/common/software/desi/desi_environment.sh main
export LSSCODE=$HOME/LSSCODE
PYTHONPATH=$PYTHONPATH:$LSSCODE/LSS/py

srun -N 1 python3 /global/homes/d/desica/LSScode/LSS/bin/SingleNode-runAltMTLParallel.py --altMTLBaseDir=/pscratch/sd/d/desica/AltMTL -obscon=DARK --survey=main --realMTLBaseDir='/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/' --multiDate --nproc 128