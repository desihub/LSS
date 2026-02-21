#!/bin/bash
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
PYTHONPATH=$PYTHONPATH:$HOME/LSS
for (( i=$1;i<=$2;i++ ))
do
 srun -N 1 -n 64 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/pkrun.py --tracer ELG_LOP_ffa --region NGC SGC --weight_type default_FKP --rebinning y --basedir /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/mock$i/ --outdir /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/mock$i/ 
  srun -N 1 -n 64 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/pkrun.py --tracer LRG_ffa --region NGC SGC --weight_type default_FKP --rebinning y --basedir /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/mock$i/ --outdir /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/mock$i/
  srun -N 1 -n 64 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/pkrun.py --tracer LRG_ffa --region NGC SGC --weight_type default_FKP_addRF --rebinning y --basedir /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/mock$i/ --outdir /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/mock$i/
  srun -N 1 -n 64 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/pkrun.py --tracer LRG_ffa --region NGC SGC --weight_type default_FKP_addSN --rebinning y --basedir /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/mock$i/ --outdir /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/mock$i/
 srun -N 1 -n 64 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/pkrun.py --tracer QSO_ffa --region NGC SGC --weight_type default_FKP --rebinning y --basedir /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/mock$i/ --outdir /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/mock$i/ 
done

