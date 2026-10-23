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
#test
cfsdir='/dvs_ro/cfs/cdirs/desi/mocks/cai/LSS/DA2/mocks/'
#holi_v4/altmtl/altmtl'+str(i)+'/loa-v1/mock'+str(i)
#echo $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey $survey --surveycat $surveycat --add_gtl y --specdata loa-v1 --tracer dark --targDir $SCRATCH/$survey/mocks/$sim --combd y --joindspec y --par y --usepota y

python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey $survey --surveycat $surveycat --add_gtl y --specdata loa-v1 --tracer dark --targDir $cfsdir/$sim --combd y --joindspec y --par y --usepota y --pota $cfsdir/$sim/altmtl$mocknum/fba$mocknum/pota-DARK.fits

