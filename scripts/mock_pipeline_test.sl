#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=desi
#SBATCH --constraint=cpu
#SBATCH -q debug
#SBATCH -t 00:30:00

source /global/common/software/desi/users/adematti/cosmodesi_environment.sh test
export LSScode=$HOME

scripts/sysnetELG_LOPnotqso_zbins.sh '' ELG_LOPnotqso /pscratch/sd/a/ajross/DA2/mocks/SecondGenMocks/AbacusSummit_v4_1/altmtl5/loa-v1/mock5/LSScats/
