#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=30:00
#SBATCH --qos=shared
#SBATCH --account=desi
#SBATCH --cpus-per-task=54
#SBATCH --constraint=cpu
#SBATCH --mem=108G
#SBATCH --array=15
#SBATCH -J getpota
#SBATCH -o ./stdout/%x_%a.o%j
#SBATCH -e ./stdout/%x_%a.e%j

source /global/common/software/desi/desi_environment.sh main
PYTHONPATH=$PYTHONPATH:$HOME/installed_packages/desihub/LSS/py

LSSCODE="$HOME/installed_packages/desihub/LSS/"

seed=$SLURM_ARRAY_TASK_ID
printf -v realization "%d" $seed

base_input="/pscratch/sd/j/jerryou/DA2/mocks/EZmock/forFA${seed}.fits"
base_output="/pscratch/sd/j/jerryou/DA2/mocks/EZmock/mock${realization}/"
tile_temp_dir="/pscratch/sd/j/jerryou/DA2/mocks/EZmock/mock${realization}/"
mkdir -p ${tile_temp_dir}

#NEED TO BE TESTED IN CASE OF I/O ISSUES
#
#Execute it as: >> sbatch abamtl_getpota_sbatch_da2.sh DARK or sbatch abamtl_getpota_sbatch_da2.sh BRIGHT
#if LSSCODE has not been sourced previously, it should be the path to https://github.com/desihub/LSS/ in your system
if [ "$1" == "DARK" ]; then
    time srun -n 1 --cpu-bind=none python $LSSCODE/scripts/mock_tools/getpotaDA2_mock.py --prog DARK --mock Generic --realization ${realization} --base_input ${base_input} --base_output ${base_output} --tile_temp_dir ${tile_temp_dir}

elif [ "$1" == "BRIGHT" ]; then
    time srun -n 1 --cpu-bind=none python $LSSCODE/scripts/mock_tools/getpotaDA2_mock.py --realization $SLURM_ARRAY_TASK_ID --mock_version BGS_v2 --prog BRIGHT    
else
    echo "Need to define DARK or BRIGHT"
fi





