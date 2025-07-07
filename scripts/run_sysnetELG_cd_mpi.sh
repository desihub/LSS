#!/bin/bash

# example to run as bash script 
# dnnp with pnll, linp with pnll or lin with mse

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
basedir=${10}
output_nn=$basedir/$version/sysnet/${tracer}_${run}
input_data=$basedir/$version/sysnet/prep_${tracer}_${run}.fits
echo using $input_data and saving to $output_nn

# nn parameters
axes="all"
nchain=${11}
nepoch=${12}
layers=${13}
units=${14}
nns=($layers $units) # structure of the neural net (# layers, # units)
bsize=$batchsize
lr=$learnrate 
model=$model #dnnp
loss=$loss #pnll
etamin=0.00001

sysnet_app=$LSSDIR/LSS/py/LSS/imaging/sysnet_app_mpi.py

if [ "${do_LRfinder}" = true ]
then
   du -h $input_data
   srun -n 1 $sysnet_app -i ${input_data} -o ${output_nn} -ax ${axes[@]} -bs ${bsize} --model $model --loss $loss --nn_structure ${nns[@]} -fl
fi

if [ "${do_nnrun}" = true ]
then
   du -h $input_data
   echo using lr ${lr[${ph}]}
   time srun -n 5 $sysnet_app -i ${input_data} -o ${output_nn} -ax ${axes[@]} -bs ${bsize} --model $model --loss $loss --nn_structure ${nns[@]} -lr ${lr[${ph}]} --eta_min $etamin -ne $nepoch -nc $nchain -k 
fi
