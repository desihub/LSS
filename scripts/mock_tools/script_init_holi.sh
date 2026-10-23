#!/bin/bash
source /global/common/software/desi/desi_environment.sh main
module load desitarget/3.0.0

for i in {0..49}
do
    printf -v padded_i "%04d" $i
    echo $i
srun -N 1 -C cpu --cpus-per-task=64 -t 00:40:00 --ntasks=1 --qos interactive --account desi python initialize_amtl_mocks_da2.py /pscratch/sd/d/desica/DA2/mocks/holi_v3/forFA$i.fits /pscratch/sd/d/desica/DA2/mocks/holi_v3/altmtl$i/ DARK
done
	 
