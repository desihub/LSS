#!/bin/bash

#arguments to run1_456.sh are realization elg_version lrg_version qso_version output_directory mock_type

max_jobs=$3
job_count=0

for ((i=$1;i<=$2;i++ ))
do
 if [ ! -f "/pscratch/sd/d/desica/DA2/mocks/holi_v3/forFA$i.fits" ]; then
     echo /pscratch/sd/d/desica/DA2/mocks/holi_v3/forFA$i.fits
     bash scripts/mock_tools/run1_456.sh $i webjax_v4.81 webjax_v4.80 webjax_v4.80 /pscratch/sd/d/desica/DA2/mocks/holi_v3 holi &

     ((job_count++))

     # Wait if we've reached max parallel jobs
     if [ "$job_count" -ge $max_jobs ]; then
         wait -n  # Wait for any one job to complete
         ((job_count--))
     fi
 fi
done

# Wait for all remaining jobs to complete
wait
