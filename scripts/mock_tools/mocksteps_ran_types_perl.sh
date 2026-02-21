#!/bin/bash

source /global/common/software/desi/desi_environment.sh main
PYTHONPATH=$PYTHONPATH:$HOME/LSS

python mkCat_mock.py --tracer LRG --mockmin $1 --mockmax $2 --survey Y1 --fulld n --fullr y --apply_veto_ran y --mkclusran y --mkclusdat n --mkclusran_allpot y --mkclusdat_allpot n --nz y

python mkCat_mock.py --tracer ELG --mockmin $1 --mockmax $2 --survey Y1 --fulld n --fullr y --apply_veto_ran y --mkclusran y --mkclusdat n --mkclusran_allpot y --mkclusdat_allpot n --nz y

python mkCat_mock.py --tracer QSO --mockmin $1 --mockmax $2 --survey Y1 --fulld n --fullr y --apply_veto_ran y --mkclusran y --mkclusdat n --mkclusran_allpot y --mkclusdat_allpot n --nz y