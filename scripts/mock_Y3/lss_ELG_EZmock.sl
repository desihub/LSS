#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=desi
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=256
#SBATCH --constraint=cpu
#SBATCH -q debug
#SBATCH -t 30:00
#SBATCH -J lss_elg
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

#the following steps are based on /desihub/LSS/scripts/mock_tools/mock_altmtl_pipeline_dark_new.txt


##get ELG fulld (--fulld y), and masked data (--apply_veto y)
##takes ~25mins, 9% cpu node memory
time srun -n 1 --cpu-bind=none python ${codepath}/mkCat_SecondGen_amtl.py --base_altmtl_dir ${base_altmtl_dir} --mockver ab_secondgen --mocknum $mocknum --survey DA2 --simName $simName --specdata ${specdata} --tracer ELG_LOP --notqso y --fulld y --apply_veto y --par y

##get ELG clustering catalogs (--mkclusdat y --mkclusran y), split them NGC/SGC (--splitGC y), refactor/add FKP weights (--nz y)
##takes ~13 minutes (8 minutes just reading input random files; 4 minutes for 9 x2), ~100% cpu node memory
time srun -n 1 --cpu-bind=none python ${codepath}/mkCat_SecondGen_amtl.py --base_altmtl_dir ${base_altmtl_dir} --mockver ab_secondgen --mocknum $mocknum --survey DA2 --simName $simName --specdata ${specdata} --tracer ELG_LOP --notqso y --mkclusdat y --mkclusran y --splitGC y --nz y --par y
