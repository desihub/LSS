#!/bin/bash

#SBATCH -A desi
#SBATCH -C cpu
#SBATCH -p debug
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --output=JOB_OUT_%x_%j.txt
#SBATCH --error=JOB_ERR_%x_%j.txt

# get environment
source /global/common/software/desi/desi_environment.sh main
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
module swap pyrecon/main pyrecon/mpi
PYTHONPATH=$PYTHONPATH:$HOME/LSS/py 

echo $1

srun -N 1 -C cpu -t 00:30:00 -q interactive python scripts/main/apply_blinding_main_fromfile_fcomp.py --type $1 --baoblind y --mkclusdat y --mkclusran y --minr 0 --maxr 18 --version v0.1 --useMPI n
srun -n 128 -N 1 -C cpu -t 01:00:00 -q interactive python scripts/main/apply_blinding_main_fromfile_fcomp.py --dorecon y --rsdblind y --fnlblind y --useMPI y --version v0.1 --minr 0 --maxr 18 --type $1
