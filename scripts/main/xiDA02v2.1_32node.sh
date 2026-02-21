#!/bin/bash

VERSPEC='guadalupe'
VER='2.1'
WT='default_FKP'
OUTDIR='/global/cfs/cdirs/desi/survey/catalogs/DA02/LSS/guadalupe/LSScats/test/xi/'

source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
PYTHONPATH=$PYTHONPATH:$HOME/LSS


srun -N 32 -C haswell -n 32 -t 04:00:00 -q interactive   python xirunpc.py --tracer ELG_LOPnotqso  --outdir $OUTDIR --survey DA02 --verspec $VERSPEC --version $VER --weight_type $WT --corr_type smu

srun -N 32 -C haswell -n 32 -t 04:00:00 -q interactive   python xirunpc.py --tracer BGS_BRIGHT --outdir $OUTDIR --survey DA02 --verspec $VERSPEC --version $VER --weight_type $WT --corr_type smu

srun -N 32 -C haswell -n 32 -t 04:00:00 -q interactive   python xirunpc.py --tracer QSO --outdir $OUTDIR --survey DA02 --verspec $VERSPEC --version $VER --weight_type $WT --corr_type smu

srun -N 32 -C haswell -n 32 -t 04:00:00 -q interactive   python xirunpc.py --tracer QSO --outdir $OUTDIR --survey DA02 --verspec $VERSPEC --version $VER --weight_type $WT --corr_type smu --zlim lowz

srun -N 32 -C haswell -n 32 -t 04:00:00 -q interactive   python xirunpc.py --tracer LRG --outdir $OUTDIR --survey DA02 --verspec $VERSPEC --version $VER --weight_type $WT --corr_type smu

srun -N 32 -C haswell -n 32 -t 04:00:00 -q interactive   python xirunpc.py --tracer LRG QSO --outdir $OUTDIR --survey DA02 --verspec $VERSPEC --version $VER --weight_type $WT --corr_type smu

srun -N 32 -C haswell -n 32 -t 04:00:00 -q interactive   python xirunpc.py --tracer QSO ELG_LOPnotqso  --outdir $OUTDIR --survey DA02 --verspec $VERSPEC --version $VER --weight_type $WT --corr_type smu

srun -N 32 -C haswell -n 32 -t 04:00:00 -q interactive   python xirunpc.py --tracer LRG ELG_LOPnotqso --outdir $OUTDIR --survey DA02 --verspec $VERSPEC --version $VER --weight_type $WT --corr_type smu


