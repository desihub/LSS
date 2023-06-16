#source /global/common/software/desi/desi_environment.sh main
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
PYTHONPATH=$PYTHONPATH:$LSSDIR/LSS/py #make sure to set $LSSDIR to wherever the LSS repo is installed, e.g., export LSSCODE=$HOME
module load cudatoolkit/11.4
module load pytorch/1.10.0


run_sysnet=$HOME/desicode/sysnetdev/scripts/app.py

export PYTHONPATH=$PYTHONPATH:$HOME/desicode/sysnetdev/

#!/bin/bash
# example to run as bash script 
# dnnp with pnll, linp with pnll or lin with mse
# bash run_sysnet.sh N ELG_LOPnotqso true false 1024 0.003 dnnp pnll true v0

# below is code for reading data and where output is
run=$1 # north or south
tracer=$2

do_LRfinder=$3 #false #for running the learning rate finder
do_nnrun=$4 # true #false

batchsize=$5
learnrate=$6
model=$7
loss=$8

version=${9} #catalog version
basedir=$10
output_nn=$basedir/$version/sysnet/${tracer}_${run}
#/global/cfs/projectdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/$version/sysnet/${tracer}_${run}
input_data=$basedir/$version/sysnet/prep_${tracer}_${run}.fits
#/global/cfs/projectdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/$version/sysnet/prep_${tracer}_${run}.fits
echo using $input_data and saving to $output_nn

# nn parameters
axes="all"
#(0 1 2 3 4 5 6 7 8 9 10 11 12 13 14) 
# ['EBV_CHIANG_SFDcorr','STARDENS','HALPHA','EBV_MPF_Mean_FW15','BETA_ML','HI', 'PSFSIZE_G','PSFSIZE_R','PSFSIZE_Z', 'PSFDEPTH_G','PSFDEPTH_R','PSFDEPTH_Z', 'GALDEPTH_G','GALDEPTH_R','GALDEPTH_Z']
#(0 1 8) EBV, STARDENS, EBV_DIFF
# (0 1 2 3 4 5 6 7 8)['EBV','STARDENS','GALDEPTH_R_extcorr','GALDEPTH_G_extcorr','GALDEPTH_Z_extcorr','PSFSIZE_R','PSFSIZE_G','PSFSIZE_Z'] + ['EBV_DIFF']
nchain=5
nepoch=200
nns=(4 20)
bsize=$batchsize #5000
lr=$learnrate #(0.003) #(0.007)
model=$model #dnnp
loss=$loss #pnll
etamin=0.00001
    
if [ "${do_LRfinder}" = true ]
then
   du -h $input_data
   python $run_sysnet -i ${input_data} -o ${output_nn} -ax ${axes[@]} -bs ${bsize} --model $model --loss $loss --nn_structure ${nns[@]} -fl
fi

if [ "${do_nnrun}" = true ]
then
   du -h $input_data
   echo using lr ${lr[${ph}]}
   python $run_sysnet -i ${input_data} -o ${output_nn} -ax ${axes[@]} -bs ${bsize} --model $model --loss $loss --nn_structure ${nns[@]} -lr ${lr[${ph}]} --eta_min $etamin -ne $nepoch -nc $nchain -k 
fi