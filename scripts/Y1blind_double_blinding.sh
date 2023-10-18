#!/bin/bash

# get environment
source /global/common/software/desi/desi_environment.sh main
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
module swap pyrecon/main pyrecon/mpi
PYTHONPATH=$PYTHONPATH:$HOME/LSS/py 


TRACER='BGS_BRIGHT-21.5' # 'LRG' 'ELG_LOPnotqso' 'QSO'
SURVEY='Y1'
VERSPEC='iron'
VER='v0.6/blinded'
# OUTDIR='/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/v0.1/blinded/xi/uendert/'
INDIR='/global/cfs/cdirs/desi/survey/catalogs/'
OUTDIR='/global/cfs/cdirs/desi/users/uendert/desi_blinding/double_blinded/'

# Double blinded cosmology
w0=-0.90
wa=0.26
fnl=20.0

srun -n 128 -N 1 -C cpu -t 04:00:00 -q interactive python scripts/main/apply_blinding_main_fromfile_fcomp.py --type $TRACER --survey $SURVEY --basedir_in $INDIR --basedir_out $OUTDIR  --verspec $VERSPEC --version $VER --baoblind y --mkclusdat n --mkclusran n --resamp y --maxr 18 --dorecon y --rsdblind y --fnlblind y --getFKP y --specified_w0 $w0 --specified_wa $wa --specified_fnl $fnl --get_par_mode specified > out.log
# srun -n 128 -N 1 -C cpu -t 04:00:00 -q interactive python scripts/main/apply_blinding_main_fromfile_fcomp.py --type $TRACER --survey $SURVEY --basedir_in $INDIR --basedir_out $OUTDIR  --verspec $VERSPEC --version $VER --baoblind y --mkclusdat n --mkclusran n --resamp y --maxr 18 --dorecon y --rsdblind y --fnlblind n --getFKP y --specified_w0 $w0 --specified_wa $wa --get_par_mode specified > out_nofnl.log
