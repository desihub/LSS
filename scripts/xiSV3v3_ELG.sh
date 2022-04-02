#!/bin/bash

VERSPEC='fuji'
VER='test'
WT='default_FKP_angular_bitwise'


srun -N 1  python xirunpc.py --tracer ELG  --verspec $VERSPEC --version $VER --weight_type $WT &
srun -N 1  python xirunpc.py --tracer ELG_HIP  --verspec $VERSPEC --version $VER --weight_type $WT &
srun -N 1  python xirunpc.py --tracer ELG_HIPnotqso  --verspec $VERSPEC --version $VER --weight_type $WT &

srun -N 1  python xirunpc.py --tracer ELG  --verspec $VERSPEC --version $VER --weight_type $WT --zlim extended_elgzmask &
srun -N 1  python xirunpc.py --tracer ELG_HIP  --verspec $VERSPEC --version $VER --weight_type $WT --zlim extended_elgzmask &
srun -N 1  python xirunpc.py --tracer ELG_HIPnotqso  --verspec $VERSPEC --version $VER --weight_type $WT --zlim extended_elgzmask &

srun -N 1  python xirunpc.py --tracer ELG  --verspec $VERSPEC --version $VER --weight_type $WT --zlim extended &
srun -N 1  python xirunpc.py --tracer ELG_HIP  --verspec $VERSPEC --version $VER --weight_type $WT --zlim extended &
srun -N 1  python xirunpc.py --tracer ELG_HIPnotqso  --verspec $VERSPEC --version $VER --weight_type $WT --zlim extended &

srun -N 1  python xirunpc.py --tracer ELG  --verspec $VERSPEC --version $VER --weight_type $WT --zlim smallshells_elgzmask &
srun -N 1  python xirunpc.py --tracer ELG_HIP  --verspec $VERSPEC --version $VER --weight_type $WT --zlim smallshells_elgzmask &
srun -N 1  python xirunpc.py --tracer ELG_HIPnotqso  --verspec $VERSPEC --version $VER --weight_type $WT --zlim smallshells_elgzmask &

srun -N 1  python xirunpc.py --tracer ELG  --verspec $VERSPEC --version $VER --weight_type $WT --zlim smallshells &
srun -N 1  python xirunpc.py --tracer ELG_HIP  --verspec $VERSPEC --version $VER --weight_type $WT --zlim smallshells &
srun -N 1  python xirunpc.py --tracer ELG_HIPnotqso  --verspec $VERSPEC --version $VER --weight_type $WT --zlim smallshells &

wait