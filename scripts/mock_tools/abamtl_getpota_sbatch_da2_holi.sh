#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --qos=regular
#SBATCH --cpus-per-task=256
#SBATCH --constraint=cpu
#SBATCH --array=500,601,602,603,604,605,606,607,608,609,610,611,612,613,614,615,616,618,619,620,621,622,623,624,625,626,627,628,629,630,631,632,633,634,635,636,637,638,639,640,641,642,643,644,645,646,648,649,650
#SBATCH --account=desi
#SBATCH --nodes=1

source /global/common/software/desi/desi_environment.sh main
#PYTHONPATH=$PYTHONPATH:$HOME/LSS/py


#NEED TO BE TESTED IN CASE OF I/O ISSUES
#
#Execute it as: >> sbatch abamtl_getpota_sbatch_da2.sh DARK or sbatch abamtl_getpota_sbatch_da2.sh BRIGHT
#if LSSCODE has not been sourced previously, it should be the path to https://github.com/desihub/LSS/ in your system
#if [ "$1" == "DARK" ]; then
srun python getpotaDA2_mock.py --realization $SLURM_ARRAY_TASK_ID --mock holi --base_output /pscratch/sd/d/desica/DA2/mocks/holi_v1/ --base_input /pscratch/sd/d/desica/DA2/mocks/holi_v1/forFA$SLURM_ARRAY_TASK_ID.fits 

#elif [ "$1" == "BRIGHT" ]; then
#	python $LSSCODE/scripts/mock_tools/getpotaDA2_mock.py --realization $SLURM_ARRAY_TASK_ID --mock_version BGS_v2 --prog BRIGHT    
#else
#    echo "Need to define DARK or BRIGHT"
#fi





