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

#echo $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey $survey --surveycat $surveycat --add_gtl y --specdata loa-v1 --tracer dark --targDir $SCRATCH/$survey/mocks/$sim --combd y --joindspec y --par y --usepota y

python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey $survey --surveycat $surveycat --add_gtl y --specdata loa-v1 --tracer dark --targDir $SCRATCH/$survey/mocks/$sim --combd y --joindspec y --par y --usepota y --pota $SCRATCH/$survey/mocks/$sim/mock$mocknum/$survey/pota-DARK.fits

python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey $survey --surveycat $surveycat --specdata loa-v1 --targDir $SCRATCH/$survey/mocks/$sim --tracer LRG --fulld y --apply_veto y --par y &
python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey $survey --surveycat $surveycat --specdata loa-v1 --targDir $SCRATCH/$survey/mocks/$sim --tracer QSO --fulld y --apply_veto y --par y &
python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey $survey --surveycat $surveycat --specdata loa-v1 --targDir $SCRATCH/$survey/mocks/$sim --tracer ELG_LOP --notqso y --fulld y --apply_veto y --par y &
python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey $survey --surveycat $surveycat --specdata loa-v1 --targDir $SCRATCH/$survey/mocks/$sim --tracer ELG --notqso y --fulld y --apply_veto y --par y &
wait


python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey $survey --surveycat $surveycat --specdata loa-v1 --tracer QSO  --mkclusdat y --mkclusran y --splitGC y --nz y --par y
python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey $survey --surveycat $surveycat --specdata loa-v1 --tracer LRG  --mkclusdat y --mkclusran y --splitGC y --nz y --par y
python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey $survey --surveycat $surveycat --specdata loa-v1 --tracer ELG_LOP --notqso y  --mkclusdat y --mkclusran y --splitGC y --nz y --par y
python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey $survey --surveycat $surveycat --specdata loa-v1 --tracer ELG --notqso y  --mkclusdat y --mkclusran y --splitGC y --nz y --par y


python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey $survey --surveycat $surveycat --specdata loa-v1 --tracer QSO  --doimlin y --replace_syscol --par y --imsys_zbin split
python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey $survey --surveycat $surveycat --specdata loa-v1 --tracer LRG  --doimlin y --replace_syscol --par y --imsys_zbin fine
python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey $survey --surveycat $surveycat --specdata loa-v1 --tracer ELG_LOP --notqso y --par y --prep4sysnet y --nran4imsys 18

python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey $survey --surveycat $surveycat --specdata loa-v1 --tracer ELG --notqso y --par y --prep4sysnet y --nran4imsys 18

