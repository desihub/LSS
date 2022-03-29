#!/bin/bash

VERSPEC='fuji'
VER='3'
WT='default_angular_bitwise'
OUTDIR='/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/fuji/LSScats/3/xi/'


srun -N 1  python xirunpc.py --type ELG --outdir $OUTDIR --verspec $VERSPEC --version $VER --weight_type $WT --nran 8 &
srun -N 1  python xirunpc.py --type ELG --outdir $OUTDIR --verspec $VERSPEC --version $VER --weight_type $WT --nran 8 --bin_type log &
srun -N 1  python xirunpc.py --type ELG_HIP --outdir $OUTDIR --verspec $VERSPEC --version $VER --weight_type $WT --nran 8 &
srun -N 1  python xirunpc.py --type ELG_HIP --outdir $OUTDIR --verspec $VERSPEC --version $VER --weight_type $WT --nran 8 --bin_type log &
srun -N 1  python xirunpc.py --type ELG_HIPnotqso --outdir $OUTDIR --verspec $VERSPEC --version $VER --weight_type $WT --nran 8 &
srun -N 1  python xirunpc.py --type ELG_HIPnotqso --outdir $OUTDIR --verspec $VERSPEC --version $VER --weight_type $WT --nran 8 --bin_type log &
srun -N 1  python xirunpc.py --type LRG --outdir $OUTDIR --verspec $VERSPEC --version $VER --weight_type $WT  &
srun -N 1  python xirunpc.py --type LRG --outdir $OUTDIR --verspec $VERSPEC --version $VER --weight_type $WT  --bin_type log &
srun -N 1  python xirunpc.py --type LRG_main --outdir $OUTDIR --verspec $VERSPEC --version $VER --weight_type $WT  &
srun -N 1  python xirunpc.py --type LRG_main --outdir $OUTDIR --verspec $VERSPEC --version $VER --weight_type $WT  --bin_type log &
srun -N 1  python xirunpc.py --type QSO --outdir $OUTDIR --verspec $VERSPEC --version $VER --weight_type $WT  &
srun -N 1  python xirunpc.py --type QSO --outdir $OUTDIR --verspec $VERSPEC --version $VER --weight_type $WT  --bin_type log &
srun -N 1  python xirunpc.py --type BGS_ANY --outdir $OUTDIR --verspec $VERSPEC --version $VER --weight_type $WT  &
srun -N 1  python xirunpc.py --type BGS_ANY --outdir $OUTDIR --verspec $VERSPEC --version $VER --weight_type $WT  --bin_type log &
srun -N 1  python xirunpc.py --type BGS_BRIGHT --outdir $OUTDIR --verspec $VERSPEC --version $VER --weight_type $WT  &
srun -N 1  python xirunpc.py --type BGS_BRIGHT --outdir $OUTDIR --verspec $VERSPEC --version $VER --weight_type $WT  --bin_type log &
wait