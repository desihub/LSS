#!/bin/bash

num=$(($1 + 1)) 

python /pscratch/sd/a/acarnero/codes/LSS/scripts/mock_tools/prepare_mocks_Y3_dark.py --mockver ab_secondgen --mockpath /global/cfs/cdirs/desi/cosmosim/SecondGenMocks/AbacusSummit/CutSky_v4_1 --realmin $1 --realmax $num --isProduction y --split_snapshot y --new_version AbacusSummit_v4_1

