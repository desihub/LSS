#!/usr/bin/env bash

# Description:
# ===========
# top level of Holi pipeline with slurm : step1-7 and submit step8
#
# job array to process chunk (size ntasks) of simulation 
# why _flat  ? 1 level of srun, no nested step to avoid resource problem

#SBATCH --ntasks=10
#SBATCH --cpus-per-task=24
#SBATCH --account=desi
#SBATCH --constraint=cpu
#SBATCH -q regular
#SBATCH -J Holi1-7
#SBATCH -t 0:59:00
#SBATCH --array=0-7
#SBATCH --output=holi1-7_%j.out
#SBATCH --error=holi1-7_%j.err
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=jcolley@lpnhe.in2p3.fr

#
# Script parameter
#
# first ID seed to process => SBATCH --array=0-x **must start with 0**
FIRST_ID=$1
LOG_DIR=$PWD
#
# Env
#
#TODO: use file parameters
LSS_DIR=/global/cfs/cdirs/desi/users/colley/LSS
#DS_DIR=/global/cfs/cdirs/desi/mocks/fa4acm/holi/webjax_v4.80
DS_DIR=/pscratch/sd/j/jcolley/holi/webjax_v4.80
HOLI_DIR=$LSS_DIR/scripts/mock_tools/holi_pipeline
cd $HOLI_DIR

# SLURM parameters with default value for interactive mode
NTASKS=${SLURM_NTASKS:-2}
NCPU_PT=${SLURM_CPUS_PER_TASK:-4}
ARRAY_RANK=${SLURM_ARRAY_TASK_ID:-0}
FIRST_ID_RANK=$((NTASKS*ARRAY_RANK + FIRST_ID))

#
# STEP 1 create catalogue
#
time srun -n $NTASKS -c $NCPU_PT \
--output="${LOG_DIR}/logs/step1-7_${FIRST_ID_RANK}_t%t.log" \
--error="${LOG_DIR}/logs/step1-7_${FIRST_ID_RANK}_t%t.log" \
./step1.sh $LSS_DIR $DS_DIR $FIRST_ID_RANK

#
# STEP 2 concatenate file input/output for BRICKMASK
#
## create file with first seed
input3=$DS_DIR/input_chunk_$ARRAY_RANK.txt
output3=$DS_DIR/output_chunk_$ARRAY_RANK.txt
cat $DS_DIR/input$FIRST_ID_RANK.txt > $input3
cat $DS_DIR/output$FIRST_ID_RANK.txt > $output3

## loop on id seed
for ((ids=1; ids<NTASKS; ids++)); do
    seed=$((ids + FIRST_ID_RANK))
    cat $DS_DIR/input$seed.txt >> $input3
    cat $DS_DIR/output$seed.txt >> $output3
    echo input$seed.txt
done


#
# STEP 3 brickmask
#

## env definition
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh dr1
module load cpu cray-fftw
export CFITSIO_DIR=/global/common/software/desi/users/naimgk/cfitsio
export LD_LIBRARY_PATH=$CFITSIO_DIR/lib:$LD_LIBRARY_PATH
#TODO: path hardcoding, use file parameter
export EXE_PATH=/global/cfs/cdirs/desi/users/colley/wd_brick/brickmask_mpi
export CONF_PATH=/global/cfs/cdirs/desi/users/colley/LSS/scripts/mock_tools/pipeline

## BRICKMASK NOTE that command line options have priority over this file.
ALL_CPU=$((NTASKS*NCPU_PT))
# output, error in log sbatch file
time srun --exclusive -n $ALL_CPU -c 1 --cpu-bind=cores $EXE_PATH/BRICKMASK -i $input3 -o $output3 -c $HOLI_DIR/brickmask.conf


#
# STEP 4 5 6 7
#
time srun -n $NTASKS -c $NCPU_PT \
--open-mode=append \
--output="${LOG_DIR}/logs/step1-7_${FIRST_ID_RANK}_t%t.log" \
--error="${LOG_DIR}/logs/step1-7_${FIRST_ID_RANK}_t%t.log" \
./step4567.sh $LSS_DIR $DS_DIR $FIRST_ID_RANK

#
# submit STEP 8 (can't used all CPUs and very long step => new job)
#
#sbatch --ntasks=$NTASKS ./slurm3_step8.sh  $LSS_DIR $DS_DIR $FIRST_ID_RANK