#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=desi
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=256
#SBATCH --constraint=cpu
##SBATCH -q debug
#SBATCH -q regular
#SBATCH -t 48:00:00
##SBATCH --dependency=afterany:46276251
#SBATCH --nodes=1
#SBATCH --reservation=finish_dr2_mocks

source /global/common/software/desi/desi_environment.sh main

module load LSS/main

export OMP_NUM_THREADS=1

#path2LSS=$LSS/bin
path2LSS=/global/homes/d/desica/LSScode/LSS/bin

ALTMTLHOME="/pscratch/sd/d/desica/DA2/mocks/holi_bgs/"
simName="altmtl{mock_number}"

printf -v outputMTLFinalDestination "$ALTMTLHOME/$simName/"

obscon='BRIGHT'
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
mockinit=200
mockend=250
mocklist=''
#zfix="None"

argstring="--altMTLBaseDir=$outputMTLFinalDestination --obscon=$obscon --survey=$survey --ProcPerNode=$ProcPerNode $numobs_from_ledger $redoFA $getosubp $debug $verbose $secondary $mock $targfile $multiDate $reproducing --mockmin=$mockinit --mockmax=$mockend --mocklist=$mocklist"
echo "argstring for dateloop"
echo $argstring

srun -n 1 --cpu-bind=none python $path2LSS/runAltMTLRealizations.py $argstring