#!/usr/bin/env bash

# from https://github.com/desihub/LSS/blob/main/scripts/mock_tools/DR2_altmtl_sbatch/altmtl_200249_holi_v3.sbatch

# Goal:
# ====
# Process holi pipeline for step runAltMTLRealizations.py 
# **script slurm free** for test, debug, optimization


#
# script parameters
#
LSS_DIR=$1
DS_DIR=$2   # root directory of mock with version
FIRST_ID=$3      # id seed to process
SIZE_CHUNK=$4     # multiprocess
NCPUTASK=$5
PROCID=${SLURM_PROCID:-0}

IDS=$((FIRST_ID+PROCID))

#NCPU=$((SIZE_CHUNK*NCPUTASK))
NCPU=$NCPUTASK
#
# environment
#
source /global/common/software/desi/desi_environment.sh main
module load LSS/main

export OMP_NUM_THREADS=1
export OMP_NUM_THREADS=$NCPU
export OMP_NUM_THREADS=2

path2LSS=$LSS_DIR/bin

#ALTMTLHOME="/pscratch/sd/d/desica/DA2/mocks/holi_v3/"
ALTMTLHOME=$DS_DIR

simName="altmtl{mock_number:04d}"

printf -v outputMTLFinalDestination "$ALTMTLHOME/$simName/"

obscon='DARK'
survey='main'
#ProcPerNode=$NCPU
ProcPerNode=1
numobs_from_ledger=''
redoFA=''
getosubp=''
debug=''
verbose='--verbose'
secondary=''
mock='--mock'
targfile="--targfile=${ALTMTLHOME}/forFA{mock_number:04d}.fits"
multiDate='--multiDate'
reproducing=''
mockinit=$IDS
#mockend=$((IDS + SIZE_CHUNK))
mockend=$((IDS + 1))
mocklist=''
zfix="${ALTMTLHOME}/qsos/qso{mock_number:04d}.txt"

argstring="--altMTLBaseDir=$outputMTLFinalDestination --obscon=$obscon --survey=$survey --ProcPerNode=$ProcPerNode $numobs_from_ledger $redoFA $getosubp $debug $verbose $secondary $mock $targfile $multiDate $reproducing --mockmin=$mockinit --mockmax=$mockend --mocklist=$mocklist --zfix=$zfix"
echo "argstring for dateloop"
echo $argstring

#srun -n 1 -c $NCPU --cpu-bind=none python $path2LSS/runAltMTLRealizations.py $argstring
python $path2LSS/runAltMTLRealizations.py $argstring > $DS_DIR/fa_chunk_${IDS}.log  2>&1