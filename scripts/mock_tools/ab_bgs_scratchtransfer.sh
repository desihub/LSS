#!/bin/bash
for (( i=$1;i<=$2;i++ ))
do
 cp /pscratch/sd/a/ajross/AbacusSummitBGS/SecondGenMocks/AbacusSummit/mock$i/* /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummitBGS/SecondGenMocks/AbacusSummit/mock$i/

done

