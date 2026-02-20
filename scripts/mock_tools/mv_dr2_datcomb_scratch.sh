#!/bin/bash

for ((i=$2;i<=$3;i++ ))
do
 echo $SCRATCH/DA2/mocks/$1/altmtl$i/initled
 mv $SCRATCH/DA2/mocks/$1/altmtl$i/fba$i/* /global/cfs/cdirs/desi/mocks/cai/LSS/DA2/mocks/$1/altmtl$i/fba$i/
 mv $SCRATCH/DA2/mocks/$1/mock$i/* /global/cfs/cdirs/desi/mocks/cai/LSS/DA2/mocks/$1/altmtl$i/fba$i/
done

