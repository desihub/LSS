#!/bin/bash
set -e
source /global/common/software/desi/desi_environment.sh main
module load LSS/main
export LSSCODE=$LSS
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
mocknum=$1
scriptdir=/global/homes/d/desica/LSScode/LSS/scripts
sim=holi_v4/altmtl
survey=DA3
surveycat=DA2 #this will make it process the DR2 footprint
#PYTHONPATH=/global/homes/d/desica/LSScode/LSS/py:$PYTHONPATH

python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey $survey --surveycat $surveycat --specdata loa-v1 --tracer ELG_LOP --notqso y --par y --addsysnet y --replace_syscol

python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey $survey --surveycat $surveycat --specdata loa-v1 --tracer ELG --notqso y --par y --addsysnet y --replace_syscol

python $scriptdir/mock_tools/mkCat_amtl.py --mocknum $mocknum --tracer LRG --simName $sim --transfer_cfs
