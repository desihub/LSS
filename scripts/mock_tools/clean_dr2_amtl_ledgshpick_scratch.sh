#!/bin/bash

for ((i=$2;i<=$3;i++ ))
do
 echo $SCRATCH/DA2/mocks/$1/altmtl$i/
 rm -r $SCRATCH/DA2/mocks/$1/altmtl$i/Univ000/main
 rm /pscratch/sd/d/desica/DA2/mocks/$1/altmtl$i/Univ000/fa/MAIN/*/*.sh
 rm /pscratch/sd/d/desica/DA2/mocks/$1/altmtl$i/Univ000/fa/MAIN/*/*.pickle
done

