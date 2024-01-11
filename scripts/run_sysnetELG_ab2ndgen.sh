#!/bin/bash

# below is code for reading data and where output is
run=$1 # north or south
tracer=$2

do_LRfinder=$3 #false #for running the learning rate finder
do_nnrun=$4 # true #false

batchsize=$5
learnrate=$6
model=$7
loss=$8

real=${9} #realization number
mockversion=${14}
if [ "${mockversion}" == "v0" ]
then
    output_nn=/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/mock$real/sysnet/${tracer}_${run}
    input_data=/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/mock$real/sysnet/prep_${tracer}_${run}.fits
else
    output_nn=/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/mock$real/$mockversion/sysnet/${tracer}_${run}
    input_data=/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/mock$real/$mockversion/sysnet/prep_${tracer}_${run}.fits
fi
echo using $input_data and saving to $output_nn


# nn parameters
axes="all"
nchain=${10}
nepoch=${11}
layers=${12}
units=${13}
nns=($layers $units) # structure of the neural net (# layers, # units)
bsize=$batchsize
lr=$learnrate
model=$model #dnnp
loss=$loss #pnll
etamin=0.00001

echo $run $tracer $do_LRfinder $do_nnrun $batchsize $learnrate
    
if [ "${do_LRfinder}" = true ]
then
   du -h $input_data
   sysnet-app -i ${input_data} -o ${output_nn} -ax ${axes[@]} -bs ${bsize} --model $model --loss $loss --nn_structure ${nns[@]} -fl
fi

if [ "${do_nnrun}" = true ]
then
   du -h $input_data
   echo using lr ${lr[${ph}]}
   sysnet-app -i ${input_data} -o ${output_nn} -ax ${axes[@]} -bs ${bsize} --model $model --loss $loss --nn_structure ${nns[@]} -lr ${lr[${ph}]} --eta_min $etamin -ne $nepoch -nc $nchain -k 
fi
