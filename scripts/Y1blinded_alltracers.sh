#!/bin/bash

source /global/common/software/desi/desi_environment.sh main
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
PYTHONPATH=$PYTHONPATH:$HOME/LSS
module swap pyrecon/main pyrecon/mpi

srun -n 128 -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/main/apply_blinding_main_fromfile_fcomp.py --type ELG_LOPnotqso --basedir_out /global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron --version $1 --baoblind y --mkclusdat y --mkclusran y --maxr 18 --dorecon y --rsdblind y --fnlblind y --getFKP y

srun -n 128 -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/main/apply_blinding_main_fromfile_fcomp.py --type LRG --basedir_out /global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron --version $1 --baoblind y --mkclusdat y --mkclusran y --maxr 18 --dorecon y --rsdblind y --fnlblind y --getFKP y

srun -n 128 -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/main/apply_blinding_main_fromfile_fcomp.py --type QSO --basedir_out /global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron --version $1 --baoblind y --mkclusdat y --mkclusran y --maxr 18 --dorecon y --rsdblind y --fnlblind y --getFKP y

srun -n 128 -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/main/apply_blinding_main_fromfile_fcomp.py --type BGS_BRIGHT-21.5 --basedir_out /global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron --version $1 --baoblind y --mkclusdat y --mkclusran y --maxr 18 --dorecon y --rsdblind y --fnlblind y --getFKP y