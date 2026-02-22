#!/bin/bash
for ((i=$1;i<=$2;i++ ))
do
 python scripts/mock_tools/mkCat_amtl.py --mocknum $i --tracer LRG --simName GLAM-Uchuu/altmtl_cov --transfer_cfs
done

