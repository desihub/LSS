#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=desi
#SBATCH --constraint=cpu
#SBATCH -q debug
#SBATCH -t 00:30:00

source /global/common/software/desi/users/adematti/cosmodesi_environment.sh test
export LSSCODE=$HOME
srun python $LSSCODE/LSS/scripts/mock_tools/mkCat_amtl.py --mockver ab_secondgen --mocknum 5 --survey DA2 --specdata loa-v1 --tracer LRG  --doimlin y --replace_syscol --par y --imsys_zbin fine
srun python $LSSCODE/LSS/scripts/mock_tools/mkCat_amtl.py --mockver ab_secondgen --mocknum 5 --survey DA2 --specdata loa-v1 --tracer ELG_LOP --notqso y --par y --prep4sysnet y
#scripts/sysnetELG_LOPnotqso_zbins.sh '' ELG_LOPnotqso /pscratch/sd/a/ajross/DA2/mocks/SecondGenMocks/AbacusSummit_v4_1/altmtl5/loa-v1/mock5/LSScats/
