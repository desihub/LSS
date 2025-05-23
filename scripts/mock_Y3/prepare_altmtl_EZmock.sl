#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=desi
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --constraint=cpu
#SBATCH -q debug
#SBATCH -t 30:00
#SBATCH -J prepare_mock
#SBATCH -o ./stdout/%x.o%j
#SBATCH -e ./stdout/%x.e%j


source /global/common/software/desi/desi_environment.sh main

survey="DA2"
nproc=128
zcol='Z'

seed=2

for tracer in "LRG" "ELG" "QSO"; do
##for tracer in "QSO"; do

   input_mockpath="/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/EZmock/${tracer}/"
   
   input_mockfile=$(printf "EZmock_%s_complete_AbacusSummit_base_c000_ph000_NScomb_%04d.fits.gz" $tracer $seed)
   echo "input: ${input_mockpath}/${input_mockfile}"
   
   #let length=${#input_mockfile}-3
   
   odir="/pscratch/sd/j/jerryou/DESI_Y3/SecondGenMocks/EZmock/${tracer}/"
   mkdir -p ${odir}
   
   ofilename=$(printf "EZmock_%s_%s_c000_ph000_NScomb_%04d.fits.gz" $tracer $survey $seed)
   output_fullpathfn="${odir}/${ofilename}"
   echo "output: ${output_fullpathfn}"
   
   
   time srun -n 1 --cpu-bind=none python prepare_EZmock_Y3.py --survey $survey --input_mockpath ${input_mockpath} --input_mockfile ${input_mockfile} --output_fullpathfn ${output_fullpathfn} --nproc $nproc --tracer $tracer --zcol $zcol

done

# check files exist
srun -n 1 --cpu-bind=none python concatenate_EZmock_targets.py $seed
