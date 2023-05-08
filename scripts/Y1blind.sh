#!/bin/bash

# get environment
source /global/common/software/desi/desi_environment.sh master
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
module swap pyrecon/main pyrecon/mpi
PYTHONPATH=$PYTHONPATH:$HOME/LSS/py 

salloc -N 1 -C cpu -t 04:00:00 --qos interactive --account desi
python scripts/main/apply_blinding_main_fromfile_fcomp.py --type LRG --baoblind y --mkclusdat y --mkclusran y --minr 0 --maxr 18 --version v0.1 --useMPI n
srun -n 128 python scripts/main/apply_blinding_main_fromfile_fcomp.py --dorecon y --rsdblind y --fnlblind y --useMPI y --version v0.1 --minr 0 --maxr 18 --type LRG