#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=desi
#SBATCH --constraint=cpu
#SBATCH -q regular
##SBATCH -q debug
#SBATCH -t 04:00:00
##SBATCH --array=601-650
##SBATCH --array=1500-1510
##SBATCH --array=0,1,2,5,7,8,9,10,11,12,13,14,23,24,25,28,29,30,31,33,34,35,38,39,43,44,45,48,49,52,53,54,56,57,58,59,60,61,62,64,65,67,68,69,91,95,100,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223
##SBATCH --array=224,228,229,230,231,232,234,235,237,238,239,240,241,242,243,245,246,247,253,255,256,257,258,261,264,266,267,268,271,272,274,275,276,277,278,280,282,283,284,287,288,289,290,291,294,296,297,300,304,306,308,310,312,313,314,316,317,318,319,320
##SBATCH --array=321,322,323,324,325,330,332,334,335,336,338,339,341,342,344,345,347,348,350,351,352,358,360,361,362,364,366,368,369,370,372,374,375,376,378,382,383,386,388,389,390,392,393,394,395,397,398,399,400,401,402,403,407,408,409,412,413,415,416,417,418,420,423,424,426,427,431,432,433,434
#SBATCH --array=435,436,438,442,445,447,448,450,471,472,473,474,475,476,477,478,479,480,481,482,483,484,485,486,487,488,489,490,491,492,493,494,495,496,497,498,499,500,601,602,604,605,611,612,613,614,615,616,617,620,621,623,624,625,628,630,633,635,636,637,638,640,641,642,644,648,649,650
#SBATCH --output=combine456test_%j.out
#SBATCH --error=combine456test_%j.err
#test

scriptdir=/global/homes/d/desica/LSScode/LSS/scripts
source /global/common/software/desi/desi_environment.sh main
module load LSS/main

mocknum=$SLURM_ARRAY_TASK_ID

elgV=v4.00
lrgV=v4.00
qsoV=v4.00
output_path=/pscratch/sd/d/desica/DA2/mocks/holi_v2

srun python $scriptdir/mock_tools/join_imaging_mask_z1.py --realization $mocknum --tracer LRG --mock holi --mock_version $lrgV
srun python $scriptdir/mock_tools/join_imaging_mask_z1.py --realization $mocknum --tracer ELG --mock holi --mock_version $elgV
srun python $scriptdir/mock_tools/join_imaging_mask_z1.py --realization $mocknum --tracer QSO --mock holi --mock_version $qsoV

srun python $scriptdir/mock_tools/add_contaminants_to_mock_z1.py --realization $mocknum --tracer QSO --mock holi --mock_version $qsoV
srun python $scriptdir/mock_tools/add_contaminants_to_mock_z1.py --realization $mocknum --tracer ELG --mock holi --mock_version $elgV

srun python $scriptdir/mock_tools/concatenate_tracers_to_fba_z1.py --realization $mocknum --mock_version_forLRG $lrgV --mock_version_forELG $elgV --mock_version_forQSO $qsoV --mock holi --output_path $output_path
