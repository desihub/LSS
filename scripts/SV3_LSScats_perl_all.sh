#!/bin/bash

python SV3/mkCat_SV3.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec fuji --version $1
python SV3/mkCat_SV3.py --type ELG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec fuji --version $1
python SV3/mkCat_SV3.py --type ELG_HIP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec fuji --version $1
python SV3/mkCat_SV3.py --type ELG_HIP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec fuji --notqso y --version $1
python SV3/mkCat_SV3.py --type ELG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec fuji --notqso y --version $1
python SV3/mkCat_SV3.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec fuji --version $1
python SV3/mkCat_SV3.py --type BGS_ANY --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec fuji --version $1
python SV3/mkCat_SV3.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec fuji --version $1

srun -N 1 -C cpu -t 01:00:00 -q interactive python SV3/mkCat_SV3_ran.py  --type LRG  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec fuji --fullr y --version $1
srun -N 1 -C cpu -t 01:00:00 -q interactive python SV3/mkCat_SV3_ran.py  --type ELG  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec fuji --fullr y --version $1
srun -N 1 -C cpu -t 01:00:00 -q interactive python SV3/mkCat_SV3_ran.py  --type ELG_HIP  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec fuji --fullr y --version $1
srun -N 1 -C cpu -t 01:00:00 -q interactive python SV3/mkCat_SV3_ran.py  --type ELG_HIP  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec fuji --fullr y --notqso y --version $1
srun -N 1 -C cpu -t 01:00:00 -q interactive python SV3/mkCat_SV3_ran.py  --type ELG  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec fuji --fullr y --notqso y --version $1
srun -N 1 -C cpu -t 01:00:00 -q interactive python SV3/mkCat_SV3_ran.py  --type QSO  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec fuji --fullr y --version $1
srun -N 1 -C cpu -t 01:00:00 -q interactive python SV3/mkCat_SV3_ran.py  --type BGS_ANY  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec fuji --fullr y --version $1
srun -N 1 -C cpu -t 01:00:00 -q interactive python SV3/mkCat_SV3_ran.py  --type BGS_BRIGHT  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec fuji --fullr y --version $1

srun -N 1 -C cpu -t 04:00:00 -q interactive python getLRGmask.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --survey SV3 --verspec fuji --maxr 18 --version $1

python SV3/mkCat_SV3.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec fuji --maxr 0 --version $1
python SV3/mkCat_SV3.py --type ELG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec fuji --maxr 0 --version $1
python SV3/mkCat_SV3.py --type ELG_HIP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec fuji --maxr 0 --version $1
python SV3/mkCat_SV3.py --type ELG_HIP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec fuji --maxr 0 --notqso y --version $1
python SV3/mkCat_SV3.py --type ELG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec fuji --maxr 0 --notqso y --version $1
python SV3/mkCat_SV3.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec fuji --maxr 0 --version $1
python SV3/mkCat_SV3.py --type BGS_ANY --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec fuji --maxr 0 --version $1
python SV3/mkCat_SV3.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec fuji --maxr 0 --version $1

srun -N 1 -C cpu -t 01:00:00 -q interactive python SV3/mkCat_SV3_ran.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec fuji --version $1
srun -N 1 -C cpu -t 01:00:00 -q interactive python SV3/mkCat_SV3_ran.py --type ELG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec fuji --version $1
srun -N 1 -C cpu -t 01:00:00 -q interactive python SV3/mkCat_SV3_ran.py --type ELG_HIP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec fuji --version $1
srun -N 1 -C cpu -t 01:00:00 -q interactive python SV3/mkCat_SV3_ran.py --type ELG_HIP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec fuji --notqso y --version $1
srun -N 1 -C cpu -t 01:00:00 -q interactive python SV3/mkCat_SV3_ran.py --type ELG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec fuji --notqso y --version $1
srun -N 1 -C cpu -t 01:00:00 -q interactive python SV3/mkCat_SV3_ran.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec fuji --version $1
srun -N 1 -C cpu -t 01:00:00 -q interactive python SV3/mkCat_SV3_ran.py --type BGS_ANY --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec fuji --version $1
srun -N 1 -C cpu -t 01:00:00 -q interactive python SV3/mkCat_SV3_ran.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec fuji --version $1

python SV3/mkCat_SV3.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --clus y --clusran y  --nz y --maxr 18 --verspec fuji --version $1
python SV3/mkCat_SV3.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --clus y --clusran y  --nz y --maxr 18 --verspec fuji --ccut main --version $1
python SV3/mkCat_SV3.py --type ELG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --clus y --clusran y  --nz y --maxr 18 --verspec fuji --version $1
python SV3/mkCat_SV3.py --type ELG_HIP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --clus y --clusran y  --nz y --maxr 18 --verspec fuji --version $1
python SV3/mkCat_SV3.py --type ELG_HIP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --clus y --clusran y  --nz y --maxr 18 --verspec fuji --notqso y --version $1
python SV3/mkCat_SV3.py --type ELG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --clus y --clusran y  --nz y --maxr 18 --verspec fuji --notqso y --version $1
python SV3/mkCat_SV3.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --clus y --clusran y  --nz y --maxr 18 --verspec fuji --version $1
python SV3/mkCat_SV3.py --type BGS_ANY --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --clus y --clusran y  --nz y --maxr 18 --verspec fuji --version $1
python SV3/mkCat_SV3.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --clus y --clusran y  --nz y --maxr 18 --verspec fuji --version $1
