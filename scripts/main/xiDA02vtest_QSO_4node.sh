#!/bin/bash

VERSPEC='guadalupe'
VER='test'
WT='default_FKP'
OUTDIR='/global/cfs/cdirs/desi/survey/catalogs/DA02/LSS/guadalupe/LSScats/test/xi/'

source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
PYTHONPATH=$PYTHONPATH:$HOME/LSS


srun -N 4 -C haswell -n 4 -t 04:00:00 -q interactive   python xirunpc.py --tracer QSO --outdir $OUTDIR --survey DA02 --verspec $VERSPEC --version $VER --weight_type $WT --corr_type smu

srun -N 4 -C haswell -n 4 -t 04:00:00 -q interactive   python xirunpc.py --tracer QSO --outdir $OUTDIR --survey DA02 --verspec $VERSPEC --version $VER --weight_type $WT --corr_type smu --zlim lowz


