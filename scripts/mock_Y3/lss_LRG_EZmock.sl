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
#SBATCH -J lss_lrg
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

#the following steps are based on /global/homes/j/jerryou/installed_packages/desihub/LSS/scripts/mock_tools/mock_altmtl_pipeline_dark_new.txt

##combine information from the assignment files (--combd y) and real data spec files (--joindspec y)
##takes about ~13mins, 20% cpu node memory
time srun -n 1 --cpu-bind=none python ${codepath}/mkCat_SecondGen_amtl.py --base_altmtl_dir ${base_altmtl_dir} --mockver ab_secondgen --mocknum $mocknum --survey DA2 --simName $simName --add_gtl y --specdata ${specdata} --tracer dark --combd y --joindspec y --par y --usepota y


##get LRG fulld (--fulld y), and masked data (--apply_veto y) and randoms (--apply_veto_ran y -->n?), and add tileloc info to randoms (--add_tlcomp y -->n?)
##takes about ~11mins, 10% cpu node memory 
time srun -n 1 --cpu-bind=none python ${codepath}/mkCat_SecondGen_amtl.py --base_altmtl_dir ${base_altmtl_dir} --mockver ab_secondgen --mocknum $mocknum --survey DA2 --simName $simName --add_gtl y --specdata ${specdata} --tracer LRG --fulld y --apply_veto y --par y


##get LRG clustering catalogs (--mkclusdat y --mkclusran y), split them NGC/SGC (--splitGC y), refactor/add FKP weights (--nz y), (--add_gtl y-->n?)
##takes about ~12mins, 85% cpu node memory
time srun -n 1 --cpu-bind=none python ${codepath}/mkCat_SecondGen_amtl.py --base_altmtl_dir ${base_altmtl_dir} --mockver ab_secondgen --mocknum $mocknum --survey DA2 --simName $simName --specdata ${specdata} --tracer LRG --mkclusdat y --mkclusran y --splitGC y --nz y --par y
