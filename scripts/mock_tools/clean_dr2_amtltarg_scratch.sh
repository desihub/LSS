#!/bin/bash

for ((i=$2;i<=$3;i++ ))
do
 echo $SCRATCH/DA2/mocks/$1/altmtl$i/Univ000/fa/MAIN/
 rm $SCRATCH/DA2/mocks/$1/altmtl$i/Univ000/fa/MAIN/*/*targ.fits
done

