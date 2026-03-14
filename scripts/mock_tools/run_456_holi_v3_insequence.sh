#!/bin/bash

#arguments are realization elg_version lrg_version qso_version output_directory mock_type

for ((i=$1;i<=$2;i++ ))
do
 bash scripts/mock_tools/run1_456.sh $i webjax_v4.81 webjax_v4.80 webjax_v4.80 /pscratch/sd/d/desica/DA2/mocks/holi_v3 holi

done