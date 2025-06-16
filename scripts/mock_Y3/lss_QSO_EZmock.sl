#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=desi
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=256
#SBATCH --constraint=cpu
#SBATCH -q debug
##SBATCH -q regular
##SBATCH -t 2:00:00
#SBATCH -t 30:00
#SBATCH -J lss_qso
#SBATCH -o ./stdout/%x.o%j
#SBATCH -e ./stdout/%x.e%j
##SBATCH --dependency=afterany:35957648

#annotated mock altmtl pipeline
source /global/common/software/desi/desi_environment.sh main
export LSSCODE=$HOME/installed_packages/desihub/LSS
PYTHONPATH=$PYTHONPATH:$LSSCODE/py

specdata="loa-v1"

mocknum=4

codepath="$LSSCODE/scripts/mock_tools"
#codepath=./

base_altmtl_dir="/pscratch/sd/j/jerryou/"
simName="EZmock"

#--par added to all steps to make sure parallel processing is used

#the following steps are based on desihub/LSS/scripts/mock_tools/mock_altmtl_pipeline_dark_new.txt

#get QSO fulld (--fulld y), and masked data (--apply_veto y)
#redundant/unnecessary steps removed
#takes ~ 5mins, 4% cpu node memory
time srun -n 1 --cpu-bind=none python ${codepath}/mkCat_SecondGen_amtl.py --base_altmtl_dir ${base_altmtl_dir} --mockver ab_secondgen --mocknum $mocknum --survey DA2 --simName $simName --specdata ${specdata} --tracer QSO --fulld y --apply_veto y --par y


#switch to cosmodesi environment for QSO
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main

#get QSO clustering catalogs (--mkclusdat y --mkclusran y), split them NGC/SGC (--splitGC y), refactor/add FKP weights (--nz y)
#takes ~17 minutes, OOM required running only 6 randoms at a time
#takes ~100% cpu node memory
time srun -n 1 --cpu-bind=none python ${codepath}/mkCat_SecondGen_amtl.py --base_altmtl_dir ${base_altmtl_dir} --mockver ab_secondgen --mocknum $mocknum --survey DA2 --simName $simName --specdata ${specdata} --tracer QSO --mkclusdat y --mkclusran y --splitGC y --nz y --par y
