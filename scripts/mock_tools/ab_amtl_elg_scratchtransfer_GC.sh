#!/bin/bash
for (( i=$1;i<=$2;i++ ))
do
 cp /pscratch/sd/a/ajross/Y1/mocks/SecondGenMocks/AbacusSummit/altmtl$i/mock$i/LSScats/*GC* /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/altmtl$i/mock$i/LSScats/

done

