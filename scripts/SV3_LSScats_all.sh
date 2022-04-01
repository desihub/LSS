#!/bin/bash

python SV3/mkCat_SV3.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec fuji
python SV3/mkCat_SV3.py --type ELG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec fuji
python SV3/mkCat_SV3.py --type ELG_HIP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec fuji
python SV3/mkCat_SV3.py --type ELG_HIP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec fuji --notqso y
python SV3/mkCat_SV3.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec fuji
python SV3/mkCat_SV3.py --type BGS_ANY --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec fuji
python SV3/mkCat_SV3.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld y --verspec fuji

srun -N 1 -C haswell -c 64 -t 01:00:00 -q interactive python SV3/mkCat_SV3_ran.py  --type LRG  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec fuji --fullr y
srun -N 1 -C haswell -c 64 -t 01:00:00 -q interactive python SV3/mkCat_SV3_ran.py  --type ELG  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec fuji --fullr y
srun -N 1 -C haswell -c 64 -t 01:00:00 -q interactive python SV3/mkCat_SV3_ran.py  --type ELG_HIP  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec fuji --fullr y
srun -N 1 -C haswell -c 64 -t 01:00:00 -q interactive python SV3/mkCat_SV3_ran.py  --type ELG_HIP  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec fuji --fullr y --notqso y
srun -N 1 -C haswell -c 64 -t 01:00:00 -q interactive python SV3/mkCat_SV3_ran.py  --type QSO  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec fuji --fullr y 
srun -N 1 -C haswell -c 64 -t 01:00:00 -q interactive python SV3/mkCat_SV3_ran.py  --type BGS_ANY  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec fuji --fullr y 
srun -N 1 -C haswell -c 64 -t 01:00:00 -q interactive python SV3/mkCat_SV3_ran.py  --type BGS_BRIGHT  --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec fuji --fullr y 

srun -N 1 -C haswell -t 04:00:00 -q interactive python getLRGmask.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --survey SV3 --verspec fuji --maxr 18 

python SV3/mkCat_SV3.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec fuji --maxr 0
python SV3/mkCat_SV3.py --type ELG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec fuji --maxr 0
python SV3/mkCat_SV3.py --type ELG_HIP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec fuji --maxr 0
python SV3/mkCat_SV3.py --type ELG_HIP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec fuji --maxr 0 --notqso y
python SV3/mkCat_SV3.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec fuji --maxr 0
python SV3/mkCat_SV3.py --type BGS_ANY --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec fuji --maxr 0
python SV3/mkCat_SV3.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n --apply_veto y --verspec fuji --maxr 0

srun -N 1 -C haswell -c 64 -t 01:00:00 -q interactive python SV3/mkCat_SV3_ran.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec fuji 
srun -N 1 -C haswell -c 64 -t 01:00:00 -q interactive python SV3/mkCat_SV3_ran.py --type ELG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec fuji 
srun -N 1 -C haswell -c 64 -t 01:00:00 -q interactive python SV3/mkCat_SV3_ran.py --type ELG_HIP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec fuji 
srun -N 1 -C haswell -c 64 -t 01:00:00 -q interactive python SV3/mkCat_SV3_ran.py --type ELG_HIP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec fuji --notqso y
srun -N 1 -C haswell -c 64 -t 01:00:00 -q interactive python SV3/mkCat_SV3_ran.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec fuji 
srun -N 1 -C haswell -c 64 -t 01:00:00 -q interactive python SV3/mkCat_SV3_ran.py --type BGS_ANY --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec fuji 
srun -N 1 -C haswell -c 64 -t 01:00:00 -q interactive python SV3/mkCat_SV3_ran.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fullr n --apply_veto y --verspec fuji 

python SV3/mkCat_SV3.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --clus y --clusran y  --nz y --maxr 18 --verspec fuji
python SV3/mkCat_SV3.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --clus y --clusran y  --nz y --maxr 18 --verspec fuji --ccut main
python SV3/mkCat_SV3.py --type ELG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --clus y --clusran y  --nz y --maxr 18 --verspec fuji
python SV3/mkCat_SV3.py --type ELG_HIP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --clus y --clusran y  --nz y --maxr 18 --verspec fuji
python SV3/mkCat_SV3.py --type ELG_HIP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --clus y --clusran y  --nz y --maxr 18 --verspec fuji --notqso y
python SV3/mkCat_SV3.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --clus y --clusran y  --nz y --maxr 18 --verspec fuji
python SV3/mkCat_SV3.py --type BGS_ANY --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --clus y --clusran y  --nz y --maxr 18 --verspec fuji
python SV3/mkCat_SV3.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --clus y --clusran y  --nz y --maxr 18 --verspec fuji
