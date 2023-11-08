#!/bin/bash
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
PYTHONPATH=$PYTHONPATH:$HOME/LSS/py
for (( i=$1;i<=$2;i++ ))
do
 srun -N 1 -n 128 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/pkrun.py --tracer ELG_LOP_ffa --region NGC SGC --weight_type default --rebinning y --basedir /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/mock$i/ --outdir /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/mock$i/ --rpcut 2.5
 srun -N 1 -n 128 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/pkrun.py --tracer ELG_LOP_ffa --region NGC SGC --weight_type default --rebinning y --basedir /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/mock$i/ --outdir /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/mock$i/ 
 srun -N 1 -n 128 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/pkrun.py --tracer LRG_ffa --region NGC SGC --weight_type default --rebinning y --basedir /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/mock$i/ --outdir /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/mock$i/ --rpcut 2.5
  srun -N 1 -n 128 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/pkrun.py --tracer LRG_ffa --region NGC SGC --weight_type default --rebinning y --basedir /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/mock$i/ --outdir /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/mock$i/
 srun -N 1 -n 128 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/pkrun.py --tracer QSO_ffa --region NGC SGC --weight_type default --rebinning y --basedir /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/mock$i/ --outdir /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/mock$i/ --rpcut 2.5
 srun -N 1 -n 128 -C cpu -t 04:00:00 --qos interactive --account desi python scripts/pkrun.py --tracer QSO_ffa --region NGC SGC --weight_type default --rebinning y --basedir /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/mock$i/ --outdir /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/mock$i/ 
done

