#!/bin/bash
for (( i=$1;i<=$2;i++ ))
do
 chgrp desi /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/EZmock/FFA/mock$i*/*
done