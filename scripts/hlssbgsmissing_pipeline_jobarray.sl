#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=desi
#SBATCH --constraint=cpu
#SBATCH -q regular
#SBATCH -t 02:00:00
#SBATCH --array=215, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 362, 363, 364, 365, 366, 367, 368, 369, 370, 500, 513, 552, 568, 600, 601, 602, 603, 604, 605, 606, 607, 608, 609, 610, 611, 612, 613, 614, 615, 616, 617, 618, 619, 620, 621, 622, 623, 624, 625, 626, 627, 628, 629, 630, 631, 632, 633, 634, 635, 636, 637, 638, 639, 640, 641, 642, 643, 644, 645, 646, 647, 648, 649, 665, 674, 729, 753, 779, 780, 781, 782, 783, 784, 785, 786, 787, 788, 789, 891, 892, 893, 894, 895, 896, 897, 898, 899, 950, 998
#SBATCH --reservation=finish_dr2_mocks
##SBATCH --dependency=afterany:52044899
#test
source /global/common/software/desi/desi_environment.sh main
module load LSS/main
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
mocknum=$SLURM_ARRAY_TASK_ID
scriptdir=/global/homes/d/desica/LSScode/LSS/scripts
sim=holi_bgs
#PYTHONPATH=/global/homes/d/desica/LSScode/LSS/py:$PYTHONPATH

srun python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey DA2 --add_gtl y --specdata loa-v1 --tracer bright --targDir $SCRATCH/DA2/mocks/$sim --combd y --joindspec y --par y --usepota y

python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey DA2 --specdata loa-v1 --targDir $SCRATCH/DA2/mocks/$sim --tracer BGS_BRIGHT --fulld y --apply_veto y --par y 

srun python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey DA2 --specdata loa-v1 --tracer BGS_BRIGHT-21.35  --mkclusdat y --mkclusran y --splitGC y --nz y --par y 

srun python $scriptdir/mock_tools/mkCat_amtl.py --base_altmtl_dir $SCRATCH --simName $sim --mocknum $mocknum --survey DA2 --specdata loa-v1 --tracer BGS_BRIGHT-21.35 --doimlin y --replace_syscol --par y --imsys_zbin split

srun python $scriptdir/mock_tools/mkCat_amtl.py --mocknum $mocknum --tracer LRG --simName $sim --transfer_cfs