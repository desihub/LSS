#!/bin/bash

for ((i=$2;i<=$3;i++ ))
do
 echo $SCRATCH/DA2/mocks/$1/altmtl$i/initled
 rm -r $SCRATCH/DA2/mocks/$1/altmtl$i/initled
 rm /pscratch/sd/d/desica/DA2/mocks/holi_v2/altmtl$i/Univ000/fa/MAIN/*/*targ.fits
done

