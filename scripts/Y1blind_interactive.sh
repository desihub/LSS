#!/bin/bash

#run all after getting interactive node
# get environment
source /global/common/software/desi/desi_environment.sh main
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
module swap pyrecon/main pyrecon/mpi
PYTHONPATH=$PYTHONPATH:$HOME/LSS/py 


python scripts/main/apply_blinding_main_fromfile_fcomp.py --type ELG_LOPnotqso --baoblind y --mkclusdat y --mkclusran y --minr 0 --maxr 18 --version v0.1 --useMPI n
srun -n 128 python scripts/main/apply_blinding_main_fromfile_fcomp.py --dorecon y --rsdblind y --fnlblind y --useMPI y --version v0.1 --minr 0 --maxr 18 --type ELG_LOPnotqso

python scripts/main/apply_blinding_main_fromfile_fcomp.py --type BGS_BRIGHT-21.5 --baoblind y --mkclusdat y --mkclusran y --minr 0 --maxr 18 --version v0.1 --useMPI n
srun -n 128 python scripts/main/apply_blinding_main_fromfile_fcomp.py --dorecon y --rsdblind y --fnlblind y --useMPI y --version v0.1 --minr 0 --maxr 18 --type BGS_BRIGHT-21.5

python scripts/main/apply_blinding_main_fromfile_fcomp.py --type QSO --baoblind y --mkclusdat y --mkclusran y --minr 0 --maxr 18 --version v0.1 --useMPI n
srun -n 128 python scripts/main/apply_blinding_main_fromfile_fcomp.py --dorecon y --rsdblind y --fnlblind y --useMPI y --version v0.1 --minr 0 --maxr 18 --type QSO

python scripts/main/apply_blinding_main_fromfile_fcomp.py --type LRG --baoblind y --mkclusdat y --mkclusran y --minr 0 --maxr 18 --version v0.1 --useMPI n
srun -n 128 python scripts/main/apply_blinding_main_fromfile_fcomp.py --dorecon y --rsdblind y --fnlblind y --useMPI y --version v0.1 --minr 0 --maxr 18 --type LRG