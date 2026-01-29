#!/bin/bash

for ((i=$3;i<=$4;i++ ))
do
 echo $SCRATCH/DA2/mocks/$1/altmtl$i/Univ000/main/$2
 rm -r $SCRATCH/DA2/mocks/$1/altmtl$i/Univ000/main/$2
done

