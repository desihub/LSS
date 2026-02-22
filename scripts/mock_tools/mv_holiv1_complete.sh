#!/bin/bash
for ((i=$1;i<=$2;i++ ))
do
 mv /pscratch/sd/d/desica/holi_v1/mock$i/* /global/cfs/cdirs/desi/mocks/cai/LSS/DA2/mocks/holi_v1/altmtl$i/loa-v1/mock$i/LSScats/
 echo $i
done

