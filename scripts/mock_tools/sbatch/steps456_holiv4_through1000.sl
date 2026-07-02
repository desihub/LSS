#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=desi
#SBATCH --constraint=cpu
#SBATCH -q regular
#SBATCH -t 01:30:00
#SBATCH --array=851-1000
#SBATCH --reservation=_CAP_holi_dr3_altmtl
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128

source /global/common/software/desi/desi_environment.sh main
module load LSS/main
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main

scriptdir=/global/u2/d/desica/LSScode/LSS/scripts/mock_tools/
mocknum=$SLURM_ARRAY_TASK_ID

python $scriptdir/join_imaging_mask_z1.py --realization $mocknum --tracer ELG --mock generic --file_name /pscratch/sd/d/desica/DA3/mocks/holi_v4/ELG/webjax_v4.81_DR3/seed{seed}/imforFA0_Y5_noimagingmask_applied.fits --specver matterhorn-v2

python $scriptdir/join_imaging_mask_z1.py --realization $mocknum --tracer QSO --mock generic --file_name /pscratch/sd/d/desica/DA3/mocks/holi_v4/QSO/webjax_v4.81_DR3/seed{seed}/imforFA0_Y5_noimagingmask_applied.fits --specver matterhorn-v2

python $scriptdir/join_imaging_mask_z1.py --realization $mocknum --tracer LRG --mock generic --file_name /pscratch/sd/d/desica/DA3/mocks/holi_v4/LRG/webjax_v4.80_DR3/seed{seed}/imforFA0_Y5_noimagingmask_applied.fits --specver matterhorn-v2

python $scriptdir/add_contaminants_to_mock_z1.py --tracer ELG --realization $mocknum --mock generic --file_name /pscratch/sd/d/desica/DA3/mocks/holi_v4/ELG/webjax_v4.81_DR3/seed{seed}/forFA0.fits --contaminant_path /global/cfs/projectdirs/desi/mocks/cai/contaminants/DA3/matterhorn-v2/test/ --maxrea 9

python $scriptdir/add_contaminants_to_mock_z1.py --tracer QSO --realization $mocknum --mock generic --file_name /pscratch/sd/d/desica/DA3/mocks/holi_v4/QSO/webjax_v4.81_DR3/seed{seed}/forFA0.fits --contaminant_path /global/cfs/projectdirs/desi/mocks/cai/contaminants/DA3/matterhorn-v2/test/ --maxrea 9

python $scriptdir/concatenate_tracers_to_fba_z1.py --realization $mocknum --mock_version_forLRG webjax_v4.80_DR3 --mock_version_forELG webjax_v4.81_DR3 --mock_version_forQSO webjax_v4.81_DR3 --mock generic --output_path /pscratch/sd/d/desica/DA3/mocks/holi_v4/altmtl --input_path /pscratch/sd/d/desica/DA3/mocks/holi_v4
