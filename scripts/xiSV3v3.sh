#script to run all SV3 xismu
#!/bin/bash
#SBATCH --qos=interactive
#SBATCH --nodes=16
#SBATCH --time=4:00
#SBATCH --constraint=haswell

verspec = 'fuji'
ver = '3'
wt = 'default_angular_bitwise'
outdir = '/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/fuji/LSScats/3/xi/'


srun -N 1 -c 64 python xirunpc.py --type ELG --outdir $outdir --verspec $verspec --version $ver --weight_type $wt --nran 8 &
srun -N 1 -c 64 python xirunpc.py --type ELG --outdir $outdir --verspec $verspec --version $ver --weight_type $wt --nran 8 --bin_type log &
srun -N 1 -c 64 python xirunpc.py --type ELG_HIP --outdir $outdir --verspec $verspec --version $ver --weight_type $wt --nran 8 &
srun -N 1 -c 64 python xirunpc.py --type ELG_HIP --outdir $outdir --verspec $verspec --version $ver --weight_type $wt --nran 8 --bin_type log &
srun -N 1 -c 64 python xirunpc.py --type ELG_HIPnotqso --outdir $outdir --verspec $verspec --version $ver --weight_type $wt --nran 8 &
srun -N 1 -c 64 python xirunpc.py --type ELG_HIPnotqso --outdir $outdir --verspec $verspec --version $ver --weight_type $wt --nran 8 --bin_type log &
srun -N 1 -c 64 python xirunpc.py --type LRG --outdir $outdir --verspec $verspec --version $ver --weight_type $wt  &
srun -N 1 -c 64 python xirunpc.py --type LRG --outdir $outdir --verspec $verspec --version $ver --weight_type $wt  --bin_type log &
srun -N 1 -c 64 python xirunpc.py --type LRG_main --outdir $outdir --verspec $verspec --version $ver --weight_type $wt  &
srun -N 1 -c 64 python xirunpc.py --type LRG_main --outdir $outdir --verspec $verspec --version $ver --weight_type $wt  --bin_type log &
srun -N 1 -c 64 python xirunpc.py --type QSO --outdir $outdir --verspec $verspec --version $ver --weight_type $wt  &
srun -N 1 -c 64 python xirunpc.py --type QSO --outdir $outdir --verspec $verspec --version $ver --weight_type $wt  --bin_type log &
srun -N 1 -c 64 python xirunpc.py --type BGS_ANY --outdir $outdir --verspec $verspec --version $ver --weight_type $wt  &
srun -N 1 -c 64 python xirunpc.py --type BGS_ANY --outdir $outdir --verspec $verspec --version $ver --weight_type $wt  --bin_type log &
srun -N 1 -c 64 python xirunpc.py --type BGS_BRIGHT --outdir $outdir --verspec $verspec --version $ver --weight_type $wt  &
srun -N 1 -c 64 python xirunpc.py --type BGS_BRIGHT --outdir $outdir --verspec $verspec --version $ver --weight_type $wt  --bin_type log &
wait