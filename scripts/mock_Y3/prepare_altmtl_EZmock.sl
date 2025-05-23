#!/bin/bash
#SBATCH --time=30:00
#SBATCH --qos=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --constraint=cpu
#SBATCH --array=15
#SBATCH --account=desi
#SBATCH --mem=36G
#SBATCH -J prep_mock
#SBATCH -o ./stdout/%x_%a.o%j
#SBATCH -e ./stdout/%x_%a.e%j


#source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
source /global/common/software/desi/desi_environment.sh main

survey="DA2"

mockname="EZmock"
nzmask="y"
ELGsplit="n"

nproc=18
ztruecol="Z_COSMO"
zrsdcol="Z"
seed=$SLURM_ARRAY_TASK_ID

codepath="/global/homes/j/jerryou/installed_packages/desihub/LSS/scripts/mock_tools/"

for tracer in "LRG" "ELG" "QSO"; do

   input_mockpath="/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/EZmock/${tracer}/"
   
   input_mockfile=$(printf "EZmock_%s_complete_AbacusSummit_base_c000_ph000_NScomb_%04d.fits.gz" $tracer $seed)
   echo "input: ${input_mockpath}/${input_mockfile}"
   
   #let length=${#input_mockfile}-3
   
   odir="/pscratch/sd/j/jerryou/DA2/mocks/EZmock/${tracer}/"
   mkdir -p ${odir}
   
   ofilename=$(printf "EZmock_%s_%s_c000_ph000_NScomb_%04d.fits.gz" $tracer $survey $seed)
   output_fullpathfn="${odir}/${ofilename}"
   echo "output: ${output_fullpathfn}"
   
   
   time srun -n 1 python ${codepath}/prepare_mocks_Y3.py --survey $survey --input_mockpath ${input_mockpath} --input_mockfile ${input_mockfile} --output_fullpathfn ${output_fullpathfn} --nproc $nproc --tracer $tracer --ztruecol $ztruecol --zrsdcol $zrsdcol --mockname $mockname --nzmask $nzmask --ELGsplit $ELGsplit

done

# concatenate targets (LRG+ELG+QSO)
srun -n 1 --cpu-bind=none python /global/homes/j/jerryou/installed_packages/desihub/LSS/scripts/mock_Y3/concatenate_EZmock_targets.py $seed
