#!/bin/bash
source /global/common/software/desi/desi_environment.sh main
module load LSS/main
export LSSCODE=$LSS
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
mocknum=0
scriptdir=/global/homes/d/desica/LSScode/LSS/scripts
sim=Uchuu-SHAM_BGS
survey=DA2
surveycat=DA2 #this will make it process the DR2 footprint
#PYTHONPATH=/global/homes/d/desica/LSScode/LSS/py:$PYTHONPATH
#test

#python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey $survey --surveycat $surveycat --add_gtl y --specdata loa-v1 --tracer bright --targDir $SCRATCH/$survey/mocks/$sim --combd y --joindspec y --par y --usepota y

#python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey $survey --surveycat $surveycat --specdata loa-v1 --targDir $SCRATCH/$survey/mocks/$sim --tracer BGS_BRIGHT --fulld y --apply_veto y --par y 

python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir /global/cfs/projectdirs/desi/mocks/cai/LSS/ --simName $sim --mocknum $mocknum --survey $survey --surveycat $surveycat --specdata loa-v1 --tracer BGS_BRIGHT-21.35  --mkclusdat y --mkclusran y --splitGC y --nz y --par y --outmd cfs

#python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir /global/cfs/projectdirs/desi/mocks/cai/LSS/ --simName $sim --mocknum $mocknum --survey $survey --surveycat $surveycat --specdata loa-v1 --tracer BGS_BRIGHT-21.35 --doimlin y --replace_syscol --par y --imsys_zbin split --outmd cfs

