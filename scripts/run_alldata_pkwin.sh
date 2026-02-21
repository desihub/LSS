#!/bin/bash
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
PYTHONPATH=$PYTHONPATH:$HOME/LSS

srun -N 1 -n 128 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/pkrun.py --tracer ELG_LOPnotqso --survey Y1 --verspec iron --version $1 --region NGC SGC --weight_type default_FKP --outdir /global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/$1 --rebinning y  --calc_win y --rpcut 2.5

srun -N 1 -n 128 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/pkrun.py --tracer ELG_LOPnotqso --survey Y1 --verspec iron --version $1 --region NGC SGC --weight_type default_FKP --outdir /global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/$1 --rebinning y  --calc_win y 

srun -N 1 -n 128 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/pkrun.py --tracer LRG --survey Y1 --verspec iron --version $1 --region NGC SGC --weight_type default_FKP --outdir /global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/$1 --rebinning y  --calc_win y --rpcut 2.5

srun -N 1 -n 128 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/pkrun.py --tracer LRG --survey Y1 --verspec iron --version $1 --region NGC SGC --weight_type default_FKP --outdir /global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/$1 --rebinning y  --calc_win y 

srun -N 1 -n 128 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/pkrun.py --tracer QSO --survey Y1 --verspec iron --version $1 --region NGC SGC --weight_type default_FKP --outdir /global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/$1 --rebinning y  --calc_win y --rpcut 2.5

srun -N 1 -n 128 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/pkrun.py --tracer QSO --survey Y1 --verspec iron --version $1 --region NGC SGC --weight_type default_FKP --outdir /global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/$1 --rebinning y  --calc_win y 

srun -N 1 -n 128 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/pkrun.py --tracer BGS_BRIGHT-21.5 --survey Y1 --verspec iron --version $1 --region NGC SGC --weight_type default_FKP --outdir /global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/$1 --rebinning y  --calc_win y --rpcut 2.5

srun -N 1 -n 128 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/pkrun.py --tracer BGS_BRIGHT-21.5 --survey Y1 --verspec iron --version $1 --region NGC SGC --weight_type default_FKP --outdir /global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/$1 --rebinning y  --calc_win y 
