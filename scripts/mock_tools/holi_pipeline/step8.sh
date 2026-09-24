#!/usr/bin/env bash

# from https://github.com/desihub/LSS/blob/main/scripts/mock_tools/DR2_altmtl_sbatch/altmtl_200249_holi_v3.sbatch

# Processing 1 seed with 1 CPU
# OpenMP disable ?

echo "================= Step 8: AltMTL"
date

#
# script parameters
#
LSS_DIR=$1
DS_DIR=$2   # root directory of mock with version
FIRST_ID=$3      # id seed to process

PROCID=${SLURM_PROCID:-0}
IDS=$((FIRST_ID+PROCID))

#
# Environment
#
HOLI_DIR=$LSS_DIR/scripts/mock_tools/holi_pipeline
export PATH=$HOLI_DIR:$PATH
#desi_env_vers=$(get_pars.py $HOLI_PARS amtl.desi_env_vers)
#source /global/common/software/desi/desi_environment.sh $desi_env_vers
source /global/common/software/desi/desi_environment.sh 26.3
# module load LSS/main
# use local package LSS, refresh after source env
export PYTHONPATH=$LSS_DIR/py:$PYTHONPATH
export PATH=$LSS_DIR/bin:$PATH
# desitarget ecsv to fits
target_dir=$(get_pars.py $HOLI_PARS TARGET_dir)
export PYTHONPATH=$target_dir/py:$PYTHONPATH
export PATH=$target_dir/bin:$PATH


export OMP_NUM_THREADS=1


ALTMTLHOME=$DS_DIR

simName="altmtl{mock_number:04d}"

printf -v outputMTLFinalDestination "$ALTMTLHOME/$simName/"

obscon='DARK'
survey='main'
ProcPerNode=1
numobs_from_ledger=''
redoFA=''
getosubp=''
debug=''
verbose=''
secondary=''
mock='--mock'
targfile="--targfile=${ALTMTLHOME}/forFA{mock_number:04d}.fits"
multiDate='--multiDate'
reproducing=''
mockid=$IDS
zfix="${ALTMTLHOME}/qsos/qso{mock_number:04d}.txt"

argstring="--altMTLBaseDir=$outputMTLFinalDestination --obscon=$obscon --survey=$survey --ProcPerNode=$ProcPerNode $numobs_from_ledger $redoFA $getosubp $debug $verbose $secondary $mock $targfile $multiDate $reproducing --mockid=$mockid --zfix=$zfix"

echo $argstring

#python $path2LSS/runAltMTLRealizations.py $argstring
runAltMTL.py $argstring

# python "$path2LSS/runAltMTLRealizations.py" $argstring &
# parent_pid=$!
# sleep 10
# worker_pid=$(pgrep -P "$parent_pid" | head -n 1)

# if [[ -z "$worker_pid" ]]; then
#   echo "No AltMTL worker process found for parent PID $parent_pid" >&2
#   wait "$parent_pid"
#   exit 1
# fi

# echo "Profiling worker PID $worker_pid"
# py-spy record --pid "$worker_pid" --duration 240 \
#   --rate 100 --format speedscope \
#   --output "profile_worker_${IDS}.json" &

# wait "$parent_pid"


# py-spy record --duration 120 --rate 140 --format speedscope \
#   --output "profile_worker_${IDS}.json" \
#   -- python runAltMTL.py $argstring