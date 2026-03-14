scriptdir=/global/homes/d/desica/LSScode/LSS/scripts
source /global/common/software/desi/desi_environment.sh main
module load LSS/main

mocknum=$1

elgV=$2
lrgV=$3
qsoV=$4
output_path=$5
mock=$6

python $scriptdir/mock_tools/join_imaging_mask_z1.py --realization $mocknum --tracer LRG --mock $mock --mock_version $lrgV
python $scriptdir/mock_tools/join_imaging_mask_z1.py --realization $mocknum --tracer ELG --mock $mock --mock_version $elgV
python $scriptdir/mock_tools/join_imaging_mask_z1.py --realization $mocknum --tracer QSO --mock $mock --mock_version $qsoV

python $scriptdir/mock_tools/add_contaminants_to_mock_z1.py --realization $mocknum --tracer QSO --mock $mock --mock_version $qsoV
python $scriptdir/mock_tools/add_contaminants_to_mock_z1.py --realization $mocknum --tracer ELG --mock $mock --mock_version $elgV

python $scriptdir/mock_tools/concatenate_tracers_to_fba_z1.py --realization $mocknum --mock_version_forLRG $lrgV --mock_version_forELG $elgV --mock_version_forQSO $qsoV --mock holi --output_path $output_path
