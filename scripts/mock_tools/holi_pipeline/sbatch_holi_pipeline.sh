#!/usr/bin/env bash

# Holi pipeline - SLURM array job driver for steps 1 to 8

# Goal
# ====
# Run the full Holi mock pipeline (catalogue creation, brickmask,
# imaging mask, contaminants, tracer concatenation, AMTL initialization
# and fiber assignment) for a batch of mock realizations (seeds).
#
# A SLURM array job spreads the seeds across array ranks; within each
# rank, "ntasks" SLURM tasks (srun) process "ntasks" seeds in parallel.


#SBATCH --array=0-0
# NOTE: the first seed ID to process is set by the "first_id" parameter
# (see get_pars.py $HOLI_PARS first_id below), not by this --array directive.
# The array range below only needs to start at 0: SBATCH --array=0-x

# ntasks is defied by memory necessary by process, here step1 arround 20-30 GB.
# actually 17 is the max
#SBATCH --ntasks=17     # number of seeds processed per array rank
#SBATCH --cpus-per-task=15 # 256/ntasks
#SBATCH --account=desi
#SBATCH --constraint=cpu
#SBATCH -q regular
## SBATCH -q shared
## SBATCH --mem-per-cpu=22G
#SBATCH -J HoliSimu
#SBATCH -t 35:00:00
#SBATCH --output=holi_%j.log
#SBATCH --error=holi_%j.log
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=jcolley@lpnhe.in2p3.fr
##SBATCH --mail-user=<user@mail.xx>

#
# Script parameter: path to the pipeline parameter file (read by get_pars.py)
#
export HOLI_PARS=$1
# NOTE: all export env variables are visible in script launched by srun.

#
# Environment
#
export LOG_DIR=$PWD
LSS_DIR=$(get_pars.py $HOLI_PARS LSS_dir)
DS_DIR=$(get_pars.py $HOLI_PARS mock_dir)

HOLI_DIR=$LSS_DIR/scripts/mock_tools/holi_pipeline
cd $HOLI_DIR

# use local package LSS
export PYTHONPATH=$LSS_DIR/py:$PYTHONPATH
export PATH=$LSS_DIR/bin:$HOLI_DIR:$PATH

FIRST_ID=$(get_pars.py $HOLI_PARS first_id)

# SLURM parameters, with fallback defaults for interactive (non-sbatch) runs
NTASKS=${SLURM_NTASKS:-2}
NCPU_PT=${SLURM_CPUS_PER_TASK:-1}
ARRAY_RANK=${SLURM_ARRAY_TASK_ID:-0}
FIRST_ID_RANK=$((NTASKS*ARRAY_RANK + FIRST_ID))

#
# STEP 1: create the mock catalogue for each seed in this rank
#
date; echo "Step 1: create catalog ELG,LRG, QSO"
time srun -n $NTASKS -c $NCPU_PT --kill-on-bad-exit=0 \
--output="${LOG_DIR}/logs/seed_${FIRST_ID_RANK}_t%t.log" \
--error="${LOG_DIR}/logs/seed_${FIRST_ID_RANK}_t%t.log" \
./step1.sh $LSS_DIR $DS_DIR $FIRST_ID_RANK

# debug step1.sh memory problem
#exit $?

#
# STEP 2: merge the per-seed input/output file lists (from step 1)
# into a single pair of files, so BRICKMASK can process the whole
# rank (all its seeds) in one call
#
date; echo "Step 2: concatenate files for BRICKMASK"
## initialize the chunk files with the first seed of this rank
input3=$DS_DIR/input_chunk_$ARRAY_RANK.txt
output3=$DS_DIR/output_chunk_$ARRAY_RANK.txt
cat $DS_DIR/input$FIRST_ID_RANK.txt > $input3
cat $DS_DIR/output$FIRST_ID_RANK.txt > $output3

## append the remaining seeds of this rank
for ((ids=1; ids<NTASKS; ids++)); do
    seed=$((ids + FIRST_ID_RANK))
    cat $DS_DIR/input$seed.txt >> $input3
    cat $DS_DIR/output$seed.txt >> $output3
    echo input$seed.txt
done


#
# STEP 3: run BRICKMASK on the merged input/output files for this rank
#
(
    date; echo "Step 3: BRICKMASK"
    ## environment setup
    source /global/common/software/desi/users/adematti/cosmodesi_environment.sh dr1
    module load cpu cray-fftw
    export CFITSIO_DIR=$(get_pars.py $HOLI_PARS brickmask.cfitsio)
    export LD_LIBRARY_PATH=$CFITSIO_DIR/lib:$LD_LIBRARY_PATH
    EXE_PATH=$(get_pars.py $HOLI_PARS brickmask.exe_dir)
    CONF_PATH=$(get_pars.py $HOLI_PARS brickmask.conf_dir)
    # use local package LSS, refresh after source env
    export PYTHONPATH=$LSS_DIR/py:$PYTHONPATH
    export PATH=$LSS_DIR/bin:$HOLI_DIR:$PATH


    ## NOTE: BRICKMASK command-line options take precedence over the values
    ## set in the configuration file (brickmask.conf).
    ALL_CPU=$((NTASKS*NCPU_PT))
    time srun --exclusive -n $ALL_CPU -c 1 --cpu-bind=cores \
    --output="${LOG_DIR}/logs/brickmask_${FIRST_ID_RANK}.log" \
    --error="${LOG_DIR}/logs/brickmask_${FIRST_ID_RANK}.log" \
    $EXE_PATH/BRICKMASK -i $input3 -o $output3 -c $CONF_PATH/brickmask.conf
)
# clean files
rm $input3 $output3

#
# STEPS 4-7: imaging mask join, contaminants, tracer concatenation and
# AMTL initialization, for each seed in this rank (see step4567.sh)
#
date; echo "Step 4-7: imaging mask join, contaminants, tracer concatenation and AMTL initialization"
time srun -n $NTASKS -c $NCPU_PT --kill-on-bad-exit=0 \
--open-mode=append \
--output="${LOG_DIR}/logs/seed_${FIRST_ID_RANK}_t%t.log" \
--error="${LOG_DIR}/logs/seed_${FIRST_ID_RANK}_t%t.log" \
./step4567.sh $LSS_DIR $DS_DIR $FIRST_ID_RANK

#
# STEP 8: Fiber Assignment
#
time srun -n $NTASKS -c $NCPU_PT --kill-on-bad-exit=0 \
--open-mode=append \
--output="${LOG_DIR}/logs/seed_${FIRST_ID_RANK}_t%t.log" \
--error="${LOG_DIR}/logs/seed_${FIRST_ID_RANK}_t%t.log" \
./step8.sh  $LSS_DIR $DS_DIR $FIRST_ID_RANK

#
# Comment: in fact I need to have all node for pipeline for memory requirement
#          so I have all CPU anyway 
#
# if (( NCPU_PT == 1 )); then
#     # "full" mode: only 1 CPU per task was requested, so step 8 can run
#     # here with srun without wasting any CPU
#     date; echo "Step 8: AltMTL"
#     time srun -n $NTASKS -c $NCPU_PT --kill-on-bad-exit=0 \
#     --open-mode=append \
#     --output="${LOG_DIR}/logs/seed_${FIRST_ID_RANK}_t%t.log" \
#     --error="${LOG_DIR}/logs/seed_${FIRST_ID_RANK}_t%t.log" \
#     ./step8.sh  $LSS_DIR $DS_DIR $FIRST_ID_RANK
# else
#     # "split" mode: more than 1 CPU per task was requested, but Fiber
#     # Assignment only uses 1 CPU per seed, so submit step 8 as a separate
#     # job with fewer CPUs per task instead of wasting the extra CPUs
#     sbatch --ntasks=$NTASKS \
#     ./sbatch8_AltMTL.sh  $HOLI_PARS $FIRST_ID_RANK $LOG_DIR
# fi
