#!/bin/bash

set -e

#annotated mock altmtl pipeline
source /global/common/software/desi/desi_environment.sh main
export LSSCODE=$HOME
PYTHONPATH=$PYTHONPATH:$LSSCODE/LSS

mocknum=$1

specrel=kibo-v1

#--par added to all steps to make sure parallel processing is used

#combine information from the assignment files (--combd y) and real data spec files (--joindspec y)
#srun -N 1 -C cpu -t 02:00:00 --qos interactive --account desi python $LSSCODE/LSS/scripts/mock_tools/mkCat_SecondGen_amtl.py --base_output /global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummit_v4_1/altmtl$mocknum --mockver ab_secondgen --mocknum $mocknum --survey DA2 --add_gtl y --specdata $specrel --tracer dark --targDir /dvs_ro/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummit_v4_1 --combd y --joindspec y --par y --outmd notscratch

#get the dark_*full_noveto.ran.fits files
srun -N 1 -C cpu -t 02:00:00 --qos interactive --account desi python $LSSCODE/LSS/scripts/mock_tools/mkCat_SecondGen_amtl.py --base_output /global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummit_v4_1/altmtl$mocknum --mockver ab_secondgen --mocknum $mocknum --survey DA2  --specdata $specrel --tracer dark --targDir /dvs_ro/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummit_v4_1 --fullr y --par y --outmd notscratch

#get LRG fulld (--fulld y), and masked data (--apply_veto y) 
srun -N 1 -C cpu -t 02:00:00 --qos interactive --account desi python $LSSCODE/LSS/scripts/mock_tools/mkCat_SecondGen_amtl.py --base_output /global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummit_v4_1/altmtl$mocknum --mockver ab_secondgen --mocknum $mocknum --survey DA2  --specdata $specrel --tracer LRG --targDir /dvs_ro/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummit_v4_1 --fulld y --apply_veto y --par y --outmd notscratch

#mask randoms (--apply_veto_ran y), and add tileloc info to randoms (--add_tlcomp y)
srun -N 1 -C cpu -t 02:00:00 --qos interactive --account desi python $LSSCODE/LSS/scripts/mock_tools/mkCat_SecondGen_amtl.py --base_output /global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummit_v4_1/altmtl$mocknum --mockver ab_secondgen --mocknum $mocknum --survey DA2  --specdata $specrel --tracer LRG --targDir /dvs_ro/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummit_v4_1 --fulld n --apply_veto n --apply_veto_ran y --add_tlcomp y --par y --outmd notscratch

#get LRG clustering catalogs (--mkclusdat y --mkclusran y), split them NGC/SGC (--splitGC y), refactor/add FKP weights (--nz y)
srun -N 1 -C cpu -t 02:00:00 --qos interactive --account desi python $LSSCODE/LSS/scripts/mock_tools/mkCat_SecondGen_amtl.py --base_output /global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummit_v4_1/altmtl$mocknum --mockver ab_secondgen --mocknum $mocknum --survey DA2  --specdata $specrel --tracer LRG --targDir /dvs_ro/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummit_v4_1 --mkclusdat y --mkclusran y --splitGC y --nz y --par y --outmd notscratch


#get QSO fulld (--fulld y), and masked data (--apply_veto y) and randoms (--apply_veto_ran y), and add tileloc info to randoms (--add_tlcomp y)
srun -N 1 -C cpu -t 02:00:00 --qos interactive --account desi python $LSSCODE/LSS/scripts/mock_tools/mkCat_SecondGen_amtl.py --base_output /global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummit_v4_1/altmtl$mocknum --mockver ab_secondgen --mocknum $mocknum --survey DA2 --specdata $specrel --tracer QSO --targDir /dvs_ro/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummit_v4_1 --fulld y --apply_veto y --apply_veto_ran y --add_tlcomp y --par y --outmd notscratch


#switch to cosmodesi environment for the rest of the codes
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
export LSSCODE=$HOME #change to wherever LSS code got cloned
PYTHONPATH=$PYTHONPATH:$LSSCODE/LSS


#get QSO clustering catalogs (--mkclusdat y --mkclusran y), split them NGC/SGC (--splitGC y), refactor/add FKP weights (--nz y)
srun -N 1 -C cpu -t 02:00:00 --qos interactive --account desi python $LSSCODE/LSS/scripts/mock_tools/mkCat_SecondGen_amtl.py --base_output /global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummit_v4_1/altmtl$mocknum --mockver ab_secondgen --mocknum $mocknum --survey DA2  --specdata $specrel --tracer QSO --targDir /dvs_ro/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummit_v4_1 --mkclusdat y --mkclusran y --splitGC y --nz y --par y --outmd notscratch

#get ELG fulld (--fulld y), and masked data (--apply_veto y) and randoms (--apply_veto_ran y), and add tileloc info to randoms (--add_tlcomp y)
srun -N 1 -C cpu -t 02:00:00 --qos interactive --account desi python $LSSCODE/LSS/scripts/mock_tools/mkCat_SecondGen_amtl.py --base_output /global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummit_v4_1/altmtl$mocknum --mockver ab_secondgen --mocknum $mocknum --survey DA2  --specdata $specrel --tracer ELG_LOP --notqso y --targDir /dvs_ro/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummit_v4_1 --fulld y --apply_veto y --apply_veto_ran y --add_tlcomp y --par y --outmd notscratch

#get ELG clustering catalogs (--mkclusdat y --mkclusran y), split them NGC/SGC (--splitGC y), refactor/add FKP weights (--nz y)
srun -N 1 -C cpu -t 02:00:00 --qos interactive --account desi python $LSSCODE/LSS/scripts/mock_tools/mkCat_SecondGen_amtl.py --base_output /global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummit_v4_1/altmtl$mocknum --mockver ab_secondgen --mocknum $mocknum --survey DA2 --specdata $specrel --tracer ELG_LOP --notqso y --targDir /dvs_ro/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummit_v4_1 --mkclusdat y --mkclusran y --splitGC y --nz y --par y --outmd notscratch
