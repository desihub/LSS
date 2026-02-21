#!/bin/bash

source /global/common/software/desi/desi_environment.sh master
PYTHONPATH=$PYTHONPATH:$HOME/LSS

srun -N 1 -C haswell -t 04:00:00 -q interactive python run_mocks_multipass.py --realmin $1 --realmax $2 --footprint Y1 --nproc 64

python mkCat_mock.py --mockmin $1 --mockmax $2 --survey Y1 --combd y --combr y --combdr y --countran y --tracer dark --add_gtl y

python mkCat_mock.py --tracer LRG --mockmin $1 --mockmax $2 --survey Y1 --fulld y --fullr y --apply_veto y --mkclusran y --mkclusdat y --mkclusran_allpot y --mkclusdat_allpot y --nz y

python mkCat_mock.py --tracer ELG --mockmin $1 --mockmax $2 --survey Y1 --fulld y --fullr y --apply_veto y --mkclusran y --mkclusdat y --mkclusran_allpot y --mkclusdat_allpot y --nz y

python mkCat_mock.py --tracer QSO --mockmin $1 --mockmax $2 --survey Y1 --fulld y --fullr y --apply_veto y --mkclusran y --mkclusdat y --mkclusran_allpot y --mkclusdat_allpot y --nz y