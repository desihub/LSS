#!/bin/bash

VERSPEC='guadalupe'
VER='2'
WT='RF_FKP_angular'

source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
PYTHONPATH=$PYTHONPATH:$HOME/LSS

srun -N 1  python xirunpc.py --tracer ELG_LOPnotqso  --survey DA02 --verspec $VERSPEC --version $VER --weight_type $WT &

srun -N 1  python xirunpc.py --tracer ELG_LOPnotqso  --survey DA02 --verspec $VERSPEC --version $VER --weight_type $WT --zlim extended_elgzmask &

srun -N 1  python xirunpc.py --tracer ELG_LOPnotqso  --survey DA02 --verspec $VERSPEC --version $VER --weight_type $WT --zlim extended &

srun -N 1  python xirunpc.py --tracer ELG_LOPnotqso  --survey DA02 --verspec $VERSPEC --version $VER --weight_type $WT --zlim smallshells_elgzmask &

srun -N 1  python xirunpc.py --tracer ELG_LOPnotqso  --survey DA02 --verspec $VERSPEC --version $VER --weight_type $WT --zlim smallshells &

wait