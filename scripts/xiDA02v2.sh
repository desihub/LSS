#!/bin/bash

source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
PYTHONPATH=$PYTHONPATH:$HOME/LSS

VERSPEC='guadalupe'
VER='2'
WT='default_FKP'
OUTDIR='/global/cfs/cdirs/desi/survey/catalogs/DA02/LSS/guadalupe/LSScats/2/xi/'


srun -N 1  python xirunpc.py --tracer ELG_LOPnotqso --survey DA02 --outdir $OUTDIR --verspec $VERSPEC --version $VER --weight_type RF_FKP &
srun -N 1  python xirunpc.py --tracer LRG --survey DA02 --outdir $OUTDIR --verspec $VERSPEC --version $VER --weight_type $WT  &
srun -N 1  python xirunpc.py --tracer QSO --survey DA02 --outdir $OUTDIR --verspec $VERSPEC --version $VER --weight_type $WT  &
srun -N 1  python xirunpc.py --tracer BGS_ANY --survey DA02 --outdir $OUTDIR --verspec $VERSPEC --version $VER --weight_type $WT  --nran 2 &
srun -N 1  python xirunpc.py --tracer BGS_BRIGHT --survey DA02 --outdir $OUTDIR --verspec $VERSPEC --version $VER --weight_type $WT --nran 2 &

srun -N 1  python xirunpc.py --tracer ELG_LOPnotqso LRG --survey DA02 --outdir $OUTDIR --verspec $VERSPEC --version $VER --weight_type RF_FKP  &
srun -N 1  python xirunpc.py --tracer ELG_LOPnotqso QSO --survey DA02 --outdir $OUTDIR --verspec $VERSPEC --version $VER --weight_type RF_FKP  &
srun -N 1  python xirunpc.py --tracer QSO LRG --survey DA02 --outdir $OUTDIR --verspec $VERSPEC --version $VER --weight_type $WT  &


wait