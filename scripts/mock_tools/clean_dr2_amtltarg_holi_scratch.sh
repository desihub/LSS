#!/bin/bash


for ((i=$1;i<=$2;i++ ))
do
 echo $SCRATCH/DA2/mocks/holi_v1/altmtl$1/Univ000/fa/MAIN/
 rm $SCRATCH/DA2/mocks/holi_v1/altmtl$1/Univ000/fa/MAIN/*/*targ.fits
done

