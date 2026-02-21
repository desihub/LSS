#!/bin/bash
#SBATCH --nodes=1
#SBATCH -q shared
#SBATCH --account=desi
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --constraint=cpu
#SBATCH -t 48:00:00
#SBATCH -J fba_onemock
#SBATCH -o ./stdout/%x.o%j
#SBATCH -e ./stdout/%x.e%j
##SBATCH --dependency=afterany:36265710

source /global/common/software/desi/desi_environment.sh main
export OMP_NUM_THREADS=1

PYTHONPATH=$PYTHONPATH:$HOME/installed_packages/desihub/LSS
path2LSS=/global/homes/j/jerryou/installed_packages/desihub/LSS/bin

ALTMTLHOME=/pscratch/sd/j/jerryou/DESI_Y3/SecondGenMocks/EZmock/
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
targfile="--targfile=/pscratch/sd/j/jerryou/DESI_Y3/SecondGenMocks/EZmock/forFA/forFA{mock_number}.fits"
multiDate='--multiDate'
reproducing=''
#mockinit=16    # 16, 80
#mockend=17
mockinit=80   
mockend=81
mocklist=''

argstring="--altMTLBaseDir=$outputMTLFinalDestination --obscon=$obscon --survey=$survey --ProcPerNode=$ProcPerNode $numobs_from_ledger $redoFA $getosubp $debug $verbose $secondary $mock $targfile $multiDate $reproducing --mockmin=$mockinit --mockmax=$mockend --mocklist=$mocklist"
echo "argstring for dateloop"
echo $argstring

srun -n 1 --cpu-bind=none python $path2LSS/runAltMTLRealizations.py $argstring
