#!/bin/bash


for ((i=$1;i<=$2;i++ ))
do
 rm -r /pscratch/sd/d/desica/DA2/mocks/holi_v1/altmtl$i/loa-v1/mock$i
done

