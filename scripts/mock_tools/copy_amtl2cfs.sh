#!/bin/bash
for ((i=$1;i<=$2;i++ ))
do
 python scripts/mock_tools/mkCat_amtl.py --mocknum $i --tracer LRG --simName $3 --transfer_cfs
done

