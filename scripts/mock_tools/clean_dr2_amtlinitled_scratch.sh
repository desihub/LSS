#!/bin/bash

for ((i=$2;i<=$3;i++ ))
do
 echo $SCRATCH/DA2/mocks/$1/altmtl$i/initled
 rm -r $SCRATCH/DA2/mocks/$1/altmtl$i/initled
done

