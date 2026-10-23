#!/bin/bash
#SBATCH --time=00:15:00
#SBATCH --qos=regular
#SBATCH --cpus-per-task=256
#SBATCH --constraint=cpu
#SBATCH --array=400-999
#SBATCH --account=desi
#SBATCH --nodes=1
#SBATCH -J potaholi
#SBATCH --reservation=altmtl

source /global/common/software/desi/desi_environment.sh main
#PYTHONPATH=$PYTHONPATH:$HOME/LSS/py


#NEED TO BE TESTED IN CASE OF I/O ISSUES
#
#Execute it as: >> sbatch abamtl_getpota_sbatch_da2.sh DARK or sbatch abamtl_getpota_sbatch_da2.sh BRIGHT
#if LSSCODE has not been sourced previously, it should be the path to https://github.com/desihub/LSS/ in your system
#if [ "$1" == "DARK" ]; then
srun python getpotaDA2_mock.py --realization $SLURM_ARRAY_TASK_ID --mock holi --base_output /pscratch/sd/d/desica/DA2/mocks/holi_v3/ --base_input /pscratch/sd/d/desica/DA2/mocks/holi_v3/forFA$SLURM_ARRAY_TASK_ID.fits 

#elif [ "$1" == "BRIGHT" ]; then
#	python $LSSCODE/scripts/mock_tools/getpotaDA2_mock.py --realization $SLURM_ARRAY_TASK_ID --mock_version BGS_v2 --prog BRIGHT    
#else
#    echo "Need to define DARK or BRIGHT"
#fi





