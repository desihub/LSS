#!/bin/bash
for (( i=$1;i<=$2;i++ ))
do
chmod 775 /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/altmtl$i/mock$i/LSScats/*
done

