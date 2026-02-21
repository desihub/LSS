#!/bin/bash
#SBATCH -N 1
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -J sys_Y1ab2ndgen
#SBATCH -t 1:00:00
#SBATCH -L SCRATCH
#SBATCH --output=/pscratch/sd/a/arosado/slurm/Y1ab2ndgen_bgs_linweight_%A_%a.out
#SBATCH --array=0-24
start_time=$(date +%s)
set -e
echo $SLURM_ARRAY_TASK_ID

source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main 
export LSSDIR=$HOME/LSS
PYTHONPATH=$PYTHONPATH:$LSSDIR

DATA_VERSION='v0.6'
TRACER=BGS_BRIGHT-21.5

if [ $TRACER == BGS_BRIGHT-21.5 ]
then
    MOCKVERSION=/SecondGenMocks/AbacusSummitBGS/
else
    MOCKVERSION=/SecondGenMocks/AbacusSummit_v4_1/
fi
BASEDIR=/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/$MOCKVERSION
ALTMTL='y'
echo $BASEDIR

WEIGHT=WEIGHT_IMLIN
VALIDATION=$LSSDIR/scripts/validation/validation_improp_mock_altmtl.py

REAL=$SLURM_ARRAY_TASK_ID
echo mock $REAL

python $LSSDIR/scripts/mock_tools/addsys.py --realization $REAL --tracer $TRACER --imsys y --add_imsys_ran y --par y --base_dir $BASEDIR --use_altmtl $ALTMTL --data_version $DATA_VERSION

python $VALIDATION --tracer $TRACER --mockn $REAL --dataver $DATA_VERSION --mockversion $MOCKVERSION --weight_col $WEIGHT 
python $VALIDATION --tracer $TRACER --mockn $REAL --dataver $DATA_VERSION --mockversion $MOCKVERSION

# Record the end time
end_time=$(date +%s)

# Calculate the difference in seconds
time_diff=$((end_time - start_time))

# Convert seconds to hours, minutes, and seconds
hours=$((time_diff / 3600))
minutes=$(( (time_diff % 3600) / 60 ))
seconds=$((time_diff % 60))

echo "Script execution time: $hours hours $minutes minutes $seconds seconds"