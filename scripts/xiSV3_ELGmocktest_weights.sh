#!/bin/bash

source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
PYTHONPATH=$PYTHONPATH:$HOME/LSS


BDIR='/global/project/projectdirs/desi/users/acarnero/mtl_mock000_univ1/'
WT='default_angular_bitwise'
WT2='default_angular_bitwise_FKP'
WT3='default_angular'
WT4='default'
WT5='default_bitwise'
OUTDIR='/global/cfs/cdirs/desi/survey/catalogs/SV3/mock/test/xi/'


srun -N 1  python xirunpc.py --tracer ELG --basedir $BDIR --outdir $OUTDIR  --weight_type $WT --nran 8 &
srun -N 1  python xirunpc.py --tracer ELG --basedir $BDIR --outdir $OUTDIR  --weight_type $WT --nran 8 --bin_type log &

srun -N 1  python xirunpc.py --tracer ELG --basedir $BDIR --outdir $OUTDIR  --weight_type $WT3 --nran 8 &
srun -N 1  python xirunpc.py --tracer ELG --basedir $BDIR --outdir $OUTDIR  --weight_type $WT3 --nran 8 --bin_type log &

srun -N 1  python xirunpc.py --tracer ELG --basedir $BDIR --outdir $OUTDIR  --weight_type $WT4 --nran 8 &
srun -N 1  python xirunpc.py --tracer ELG --basedir $BDIR --outdir $OUTDIR  --weight_type $WT4 --nran 8 --bin_type log &

srun -N 1  python xirunpc.py --tracer ELG --basedir $BDIR --outdir $OUTDIR  --weight_type $WT5 --nran 8 &
srun -N 1  python xirunpc.py --tracer ELG --basedir $BDIR --outdir $OUTDIR  --weight_type $WT5 --nran 8 --bin_type log &


wait