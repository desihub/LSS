#!/bin/bash

for ((i=$2;i<=$3;i++ ))
do
 echo $SCRATCH/DA3/mocks/$1/altmtl$i/initled
 rm -r $SCRATCH/DA3/mocks/$1/altmtl$i/initled
 rm /pscratch/sd/d/desica/DA3/mocks/$1/altmtl$i/Univ000/fa/MAIN/*/*targ.fits
 rm /pscratch/sd/d/desica/DA3/mocks/$1/altmtl$i/Univ000/fa/MAIN/*/*.sh
 rm /pscratch/sd/d/desica/DA3/mocks/$1/altmtl$i/Univ000/fa/MAIN/*/*.pickle
done

