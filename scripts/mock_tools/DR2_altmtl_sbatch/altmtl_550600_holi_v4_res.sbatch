#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=desi
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=256
#SBATCH --constraint=cpu
##SBATCH -q debug
#SBATCH -q regular
#SBATCH -t 72:00:00
#SBATCH -o /pscratch/sd/d/desica/DA3/mocks/holi_v4/altmtl/%x.o%j
#SBATCH -e /pscratch/sd/d/desica/DA3/mocks/holi_v4/altmtl/%x.e%j
##SBATCH --dependency=afterany:46276251
#SBATCH --nodes=1
#SBATCH --reservation=altmtl

source /global/common/software/desi/desi_environment.sh main

module load LSS/main

export OMP_NUM_THREADS=1

path2LSS=$LSS/bin

ALTMTLHOME="/pscratch/sd/d/desica/DA3/mocks/holi_v4/altmtl/"
simName="altmtl{mock_number}"

printf -v outputMTLFinalDestination "$ALTMTLHOME/$simName/"

obscon='DARK'
survey='main'
ProcPerNode=128
numobs_from_ledger=''
redoFA=''
getosubp=''
debug=''
verbose='--verbose'
secondary=''
mock='--mock'
targfile="--targfile=${ALTMTLHOME}/forFA{mock_number}.fits"
multiDate='--multiDate'
reproducing=''
mockinit=550
mockend=600
mocklist=''
zfix="${ALTMTLHOME}/qsos/qso{mock_number}.txt"

argstring="--altMTLBaseDir=$outputMTLFinalDestination --obscon=$obscon --survey=$survey --ProcPerNode=$ProcPerNode $numobs_from_ledger $redoFA $getosubp $debug $verbose $secondary $mock $targfile $multiDate $reproducing --mockmin=$mockinit --mockmax=$mockend --mocklist=$mocklist --zfix=$zfix"
echo "argstring for dateloop"
echo $argstring

srun -n 1 --cpu-bind=none python $path2LSS/runAltMTLRealizations.py $argstring