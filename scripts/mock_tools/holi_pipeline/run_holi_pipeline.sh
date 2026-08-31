#!/usr/bin/env bash

# Use:
#  In python 3.11 env:
#  $ run_holi_pipeline.sh <absolute/path/file/par1> ...   
#
# Parameters :
#  1. absolute path and name of holi pipeline parameters in toml format, like holi_params.toml
#
# Processing:
#  1. Create log directory with date, hour of creation in directory log file
#  2. Copy file parameter and slumr array job sbatch_holi_pipeline.sh
#  3. Launch sbatch_holi_pipeline.sh


# set env
HOLI_PARS=$1
LSS_DIR=$(get_pars.py $HOLI_PARS LSS_dir)
LOGS_DIR=$(get_pars.py $HOLI_PARS logs_dir)

# create log dir
LOG_DIR="${LOGS_DIR}/$(date +holi_%y%m%d_%Hh%M)"
mkdir -p "$LOG_DIR"
mkdir -p "$LOG_DIR/logs"

# copy pipeline parameters 
cp $HOLI_PARS $LOG_DIR
cp $LSS_DIR/sbatch_holi_pipeline.sh $LOG_DIR

# launch SLURM array job 
cd $LOG_DIR
holi_pipe=$LSS_DIR/scripts/mock_tools/holi_pipeline/sbatch_holi_pipeline.sh 
echo "Launch Holi pipeline : $holi_pipe"
echo "       Log directory : $LOG_DIR"
sbatch $holi_pipe $HOLI_PARS
cd -
