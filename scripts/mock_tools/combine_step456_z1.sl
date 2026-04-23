#!/bin/bash
#SBATCH --cpus-per-task=32
#SBATCH --account=desi
#SBATCH --constraint=cpu
#SBATCH -q regular
#SBATCH -t 01:00:00
#SBATCH --array=50-199
#SBATCH --output=c456_Hv3_%j.out
#SBATCH --error=c456_Hv3_%j.err

scriptdir=/global/homes/d/desica/LSScode/LSS/scripts
source /global/common/software/desi/desi_environment.sh main
module load LSS/main

mocknum=$SLURM_ARRAY_TASK_ID

elgV=webjax_v4.81
lrgV=webjax_v4.80
qsoV=webjax_v4.80
output_path=/pscratch/sd/d/desica/DA2/mocks/holi_v3

srun python $scriptdir/mock_tools/join_imaging_mask_z1.py --realization $mocknum --tracer LRG --mock holi --mock_version $lrgV
srun python $scriptdir/mock_tools/join_imaging_mask_z1.py --realization $mocknum --tracer ELG --mock holi --mock_version $elgV
srun python $scriptdir/mock_tools/join_imaging_mask_z1.py --realization $mocknum --tracer QSO --mock holi --mock_version $qsoV

srun python $scriptdir/mock_tools/add_contaminants_to_mock_z1.py --realization $mocknum --tracer QSO --mock holi --mock_version $qsoV
srun python $scriptdir/mock_tools/add_contaminants_to_mock_z1.py --realization $mocknum --tracer ELG --mock holi --mock_version $elgV

srun python $scriptdir/mock_tools/concatenate_tracers_to_fba_z1.py --realization $mocknum --mock_version_forLRG $lrgV --mock_version_forELG $elgV --mock_version_forQSO $qsoV --mock holi --output_path $output_path
