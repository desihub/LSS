#!/bin/bash
for i in {1..24} 
do
srun -N 1 -C cpu -t 03:50:00 -q interactive --mem=0 python initialize_amtl_mocks_da2.py $DESI_ROOT/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummit_v4_1/forFA$i.fits $DESI_ROOT/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummit_v4_1/altmtl$i DARK
done
