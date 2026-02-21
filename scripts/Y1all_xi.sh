#!/bin/bash
#make sure to set $LSSCODE to wherever the LSS git repo is (e.g., $HOME)
#provide the version, e.g. Y1all_xi.sh v0.1/blinded
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
PYTHONPATH=$PYTHONPATH:$HOME/LSS


srun -N 1 -C gpu -t 00:30:00 --gpus 4 --qos interactive --account desi_g python scripts/xirunpc.py --tracer QSO --survey Y1 --nran 3 --verspec iron --version $1 --region NGC SGC --corr_type smu --weight_type default_FKP  --njack 0 --outdir /global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/$1/AJR/xi/

srun -N 1 -C gpu -t 00:30:00 --gpus 4 --qos interactive --account desi_g python scripts/xirunpc.py --tracer LRG --survey Y1 --nran 5 --verspec iron --version $1 --region NGC SGC --corr_type smu --weight_type default_FKP  --njack 0 --outdir /global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/$1/AJR/xi/

srun -N 1 -C gpu -t 00:30:00 --gpus 4 --qos interactive --account desi_g python scripts/xirunpc.py --tracer ELG_LOPnotqso --nran 6 --survey Y1 --verspec iron --version $1 --region NGC SGC --corr_type smu --weight_type default_FKP  --njack 0 --outdir /global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/$1/AJR/xi/

srun -N 1 -C gpu -t 00:30:00 --gpus 4 --qos interactive --account desi_g python scripts/xirunpc.py --tracer BGS_BRIGHT-21.5 --nran 1 --survey Y1 --verspec iron --version $1 --region NGC SGC --corr_type smu --weight_type default_FKP --njack 0 --outdir /global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/$1/AJR/xi/
