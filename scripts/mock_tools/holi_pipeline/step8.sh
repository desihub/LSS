#!/usr/bin/env bash

# from https://github.com/desihub/LSS/blob/main/scripts/mock_tools/DR2_altmtl_sbatch/altmtl_200249_holi_v3.sbatch

# Processing 1 seed with 1 CPU
# OpenMP disable ?


#
# script parameters
#
LSS_DIR=$1
DS_DIR=$2   # root directory of mock with version
FIRST_ID=$3      # id seed to process

PROCID=${SLURM_PROCID:-0}
IDS=$((FIRST_ID+PROCID))

#
# environment
#
source /global/common/software/desi/desi_environment.sh main
module load LSS/main

export OMP_NUM_THREADS=1

path2LSS=$LSS_DIR/bin
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
mockinit=$IDS
mockend=$((IDS + 1))
mocklist=''
zfix="${ALTMTLHOME}/qsos/qso{mock_number:04d}.txt"

argstring="--altMTLBaseDir=$outputMTLFinalDestination --obscon=$obscon --survey=$survey --ProcPerNode=$ProcPerNode $numobs_from_ledger $redoFA $getosubp $debug $verbose $secondary $mock $targfile $multiDate $reproducing --mockmin=$mockinit --mockmax=$mockend --mocklist=$mocklist --zfix=$zfix"
echo "argstring for dateloop"
echo $argstring

python $path2LSS/runAltMTLRealizations.py $argstring > $DS_DIR/fa_chunk_${IDS}.log  2>&1