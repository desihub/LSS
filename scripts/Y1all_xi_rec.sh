#!/bin/bash
#make sure to set $LSSCODE to wherever the LSS git repo is (e.g., $HOME)
#provide the version, e.g. Y1all_xi.sh v0.1/blinded
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
PYTHONPATH=$PYTHONPATH:$HOME/LSS

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/xirunpc.py --tracer LRG --rec_type $2 --survey Y1 --verspec iron --version $1 --region NGC SGC --corr_type smu --weight_type default_FKP --nthreads 128 --njack 0 --recon_dir $3 --outdir /global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/$1/$3/xi/

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/xirunpc.py --tracer ELG_LOPnotqso --rec_type $2 --survey Y1 --verspec iron --version $1 --region NGC SGC --corr_type smu --weight_type default_FKP --nthreads 128 --njack 0 --recon_dir $3 --outdir /global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/$1/$3/xi/

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/xirunpc.py --tracer BGS_BRIGHT-21.5 --rec_type $2 --survey Y1 --verspec iron --version $1 --region NGC SGC --corr_type smu --weight_type default_FKP --nthreads 128 --njack 0 --recon_dir $3 --outdir /global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/$1/$3/xi/

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/xirunpc.py --tracer QSO --rec_type $2 --survey Y1 --verspec iron --version $1 --region NGC SGC --corr_type smu --weight_type default_FKP --nthreads 128 --njack 0 --recon_dir $3 --outdir /global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/$1/$3/xi/
