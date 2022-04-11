#!/bin/bash

VERSPEC='guadalupe'
VER='2'
WT='default_FKP'
WTe='RF_FKP'

source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
PYTHONPATH=$PYTHONPATH:$HOME/LSS/py

srun -N 32 -C haswell -n 1 -t 04:00:00 -q interactive   python xirunpc.py --tracer BGS_BRIGHT --survey DA02 --verspec $VERSPEC --version $VER --weight_type $WT 

srun -N 32 -C haswell -n 1 -t 04:00:00 -q interactive   python xirunpc.py --tracer BGS_ANY --survey DA02 --verspec $VERSPEC --version $VER --weight_type $WT 

srun -N 32 -C haswell -n 1 -t 04:00:00 -q interactive   python xirunpc.py --tracer QSO --survey DA02 --verspec $VERSPEC --version $VER --weight_type $WT 

srun -N 32 -C haswell -n 1 -t 04:00:00 -q interactive   python xirunpc.py --tracer LRG --survey DA02 --verspec $VERSPEC --version $VER --weight_type $WT 

srun -N 32 -C haswell -n 1 -t 04:00:00 -q interactive   python xirunpc.py --tracer ELG_LOPnotqso  --survey DA02 --verspec $VERSPEC --version $VER --weight_type $WTe 

srun -N 32 -C haswell -n 1 -t 04:00:00 -q interactive   python xirunpc.py --tracer LRG QSO --survey DA02 --verspec $VERSPEC --version $VER --weight_type $WT 

srun -N 32 -C haswell -n 1 -t 04:00:00 -q interactive   python xirunpc.py --tracer QSO ELG_LOPnotqso  --survey DA02 --verspec $VERSPEC --version $VER --weight_type $WTe 

srun -N 32 -C haswell -n 1 -t 04:00:00 -q interactive   python xirunpc.py --tracer LRG ELG_LOPnotqso  --survey DA02 --verspec $VERSPEC --version $VER --weight_type $WTe 

