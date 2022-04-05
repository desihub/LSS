#!/bin/bash

source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
PYTHONPATH=$PYTHONPATH:$HOME/LSS/py

VERSPEC='guadalupe'
VER='2'
WT='default_FKP'
OUTDIR='/global/cfs/cdirs/desi/survey/catalogs/DA02/LSS/guadalupe/LSScats/2/xi/'


srun -N 1  python xirunpc.py --tracer ELG_LOPnotqso --outdir $OUTDIR --verspec $VERSPEC --version $VER --weight_type RF_FKP &
srun -N 1  python xirunpc.py --tracer LRG --outdir $OUTDIR --verspec $VERSPEC --version $VER --weight_type $WT  &
srun -N 1  python xirunpc.py --tracer QSO --outdir $OUTDIR --verspec $VERSPEC --version $VER --weight_type $WT  &
srun -N 1  python xirunpc.py --tracer BGS_ANY --outdir $OUTDIR --verspec $VERSPEC --version $VER --weight_type $WT  --nran 2 &
srun -N 1  python xirunpc.py --tracer BGS_BRIGHT --outdir $OUTDIR --verspec $VERSPEC --version $VER --weight_type $WT --nran 2 &

srun -N 1  python xirunpc.py --tracer ELG_LOPnotqso LRG --outdir $OUTDIR --verspec $VERSPEC --version $VER --weight_type RF_FKP  &
srun -N 1  python xirunpc.py --tracer ELG_LOPnotqso QSO --outdir $OUTDIR --verspec $VERSPEC --version $VER --weight_type RF_FKP  &
srun -N 1  python xirunpc.py --tracer QSO LRG --outdir $OUTDIR --verspec $VERSPEC --version $VER --weight_type $WT  &


wait