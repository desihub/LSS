#!/bin/bash
set -e
source /global/common/software/desi/desi_environment.sh main
module load LSS/main
export LSSCODE=$LSS
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
mocknum=$1
scriptdir=/global/homes/d/desica/LSScode/LSS/scripts
sim=holi_v4/altmtl
survey=DA3
surveycat=DA2 #this will make it process the DR2 footprint
#PYTHONPATH=/global/homes/d/desica/LSScode/LSS/py:$PYTHONPATH


$scriptdir/sysnetELG_LOPnotqso_zbins.sh '' ELG_LOPnotqso $SCRATCH/DA2/mocks/$sim/altmtl$mocknum/loa-v1/mock$mocknum/LSScats/
$scriptdir/sysnetELG_LOPnotqso_zbins.sh '' ELGnotqso $SCRATCH/DA2/mocks/$sim/altmtl$mocknum/loa-v1/mock$mocknum/LSScats/

