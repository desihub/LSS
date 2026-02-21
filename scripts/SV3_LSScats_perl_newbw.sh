#!/bin/bash

source /global/common/software/desi/desi_environment.sh master
PYTHONPATH=$PYTHONPATH:$HOME/LSS

python SV3/mkCat_SV3.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec fuji --version $1
python SV3/mkCat_SV3.py --type ELG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec fuji --version $1
python SV3/mkCat_SV3.py --type ELG_HIP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec fuji --version $1
python SV3/mkCat_SV3.py --type ELG_HIP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec fuji --notqso y --version $1
python SV3/mkCat_SV3.py --type ELG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec fuji --notqso y --version $1
python SV3/mkCat_SV3.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec fuji --version $1
python SV3/mkCat_SV3.py --type BGS_ANY --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec fuji --version $1
python SV3/mkCat_SV3.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec fuji --version $1


srun -N 1 -C cpu -t 04:00:00 -q interactive python getmask_type.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --tracer LRG --survey SV3 --verspec fuji --maxr 0 --version $1 --mver 1.1
srun -N 1 -C cpu -t 04:00:00 -q interactive python getmask_type.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --tracer ELG --survey SV3 --verspec fuji --maxr 0 --version $1 --mver 1
srun -N 1 -C cpu -t 04:00:00 -q interactive python getmask_type.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --tracer ELG_HIP --survey SV3 --verspec fuji --maxr 0 --version $1 --mver 1
srun -N 1 -C cpu -t 04:00:00 -q interactive python getmask_type.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --tracer ELG_HIPnotqso --survey SV3 --verspec fuji --maxr 0 --version $1 --mver 1
srun -N 1 -C cpu -t 04:00:00 -q interactive python getmask_type.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --tracer ELGnotqso --survey SV3 --verspec fuji --maxr 0 --version $1 --mver 1

python SV3/mkCat_SV3.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec fuji --maxr 0 --version $1
python SV3/mkCat_SV3.py --type ELG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec fuji --maxr 0 --version $1
python SV3/mkCat_SV3.py --type ELG_HIP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec fuji --maxr 0 --version $1
python SV3/mkCat_SV3.py --type ELG_HIP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec fuji --maxr 0 --notqso y --version $1
python SV3/mkCat_SV3.py --type ELG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec fuji --maxr 0 --notqso y --version $1
python SV3/mkCat_SV3.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec fuji --maxr 0 --version $1
python SV3/mkCat_SV3.py --type BGS_ANY --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec fuji --maxr 0 --version $1
python SV3/mkCat_SV3.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec fuji --maxr 0 --version $1


python SV3/mkCat_SV3.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --clus y --clusran y  --nz y --maxr 18 --verspec fuji --version $1
python SV3/mkCat_SV3.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --clus y --clusran y  --nz y --maxr 18 --verspec fuji --ccut main --version $1
python SV3/mkCat_SV3.py --type ELG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --clus y --clusran y  --nz y --maxr 18 --verspec fuji --version $1
python SV3/mkCat_SV3.py --type ELG_HIP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --clus y --clusran y  --nz y --maxr 18 --verspec fuji --version $1
python SV3/mkCat_SV3.py --type ELG_HIP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --clus y --clusran y  --nz y --maxr 18 --verspec fuji --notqso y --version $1
python SV3/mkCat_SV3.py --type ELG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --clus y --clusran y  --nz y --maxr 18 --verspec fuji --notqso y --version $1
python SV3/mkCat_SV3.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --clus y --clusran y  --nz y --maxr 18 --verspec fuji --version $1
python SV3/mkCat_SV3.py --type BGS_ANY --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --clus y --clusran y --add_ke y --nz y --maxr 18 --verspec fuji --version $1
python SV3/mkCat_SV3.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --clus y --clusran y --add_ke y  --nz y --maxr 18 --verspec fuji --version $1

