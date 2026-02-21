#!/bin/bash
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
PYTHONPATH=$PYTHONPATH:$HOME/LSS

for (( i=$1;i<=$2;i++ ))
do
 indir=/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/altmtl$i/mock$i/LSScats
 srun -N 1 -C gpu -t 00:20:00 --gpus 4 --qos debug --account desi_g python scripts/xirunpc.py --gpu --tracer LRG --region NGC SGC --corr_type smu --weight_type default_FKP --njack 0 --basedir $indir --outdir $indir/xi/
 srun -N 1 -C gpu -t 00:20:00 --gpus 4 --qos debug --account desi_g python scripts/xirunpc.py --gpu --tracer ELG_LOPnotqso --region NGC SGC --corr_type smu --weight_type default_FKP --njack 0 --basedir $indir --outdir $indir/xi/
 srun -N 1 -C gpu -t 00:20:00 --gpus 4 --qos debug --account desi_g python scripts/xirunpc.py --gpu --tracer QSO --region NGC SGC --corr_type smu --weight_type default_FKP --njack 0 --basedir $indir --outdir $indir/xi/

done

