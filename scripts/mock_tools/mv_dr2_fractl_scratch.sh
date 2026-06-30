#!/bin/bash

for ((i=$2;i<=$3;i++ ))
do
 echo $SCRATCH/DA2/mocks/$1/altmtl$i/
 mv $SCRATCH/DA2/mocks/$1/altmtl$i/loa-v1/mock$i/LSScats/*tlobs* /global/cfs/cdirs/desi/mocks/cai/LSS/DA2/mocks/$1/altmtl$i/loa-v1/mock$i/LSScats/
done

