#!/bin/bash

python SV3/mkCat_SV3.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --clus y --clusran y  --nz y --maxr 18
python SV3/mkCat_SV3.py --type LRG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --clus y --clusran y  --nz y --maxr 18 --ccut main
python SV3/mkCat_SV3.py --type ELG --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --clus y --clusran y  --nz y --maxr 18
python SV3/mkCat_SV3.py --type ELG_HIP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --clus y --clusran y  --nz y --maxr 18
python SV3/mkCat_SV3.py --type ELG_HIP --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --clus y --clusran y  --nz y --maxr 18 --notqso y
python SV3/mkCat_SV3.py --type QSO --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --clus y --clusran y  --nz y --maxr 18
python SV3/mkCat_SV3.py --type BGS_ANY --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --clus y --clusran y  --nz y --maxr 18
python SV3/mkCat_SV3.py --type BGS_BRIGHT --basedir /global/cfs/cdirs/desi/survey/catalogs/  --fulld n  --clus y --clusran y  --nz y --maxr 18
