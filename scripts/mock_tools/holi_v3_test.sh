#!/bin/bash
#SBATCH --nodes=1

source /global/common/software/desi/desi_environment.sh main
module load LSS/main

declare -A conf_version

conf_version[QSO]="webjax_v4.80"
conf_version[ELG]="webjax_v4.81"
conf_version[LRG]="webjax_v4.80"

declare -A sversion
sversion[QSO]="v4.80"
sversion[ELG]="v4.81"
sversion[LRG]="v4.80"


#echo $version


declare -A nzfile

nzfile[QSO]="/global/cfs/cdirs/desi/mocks/cai/holi/webjax_v4.80/nzref_da2_qso.txt"
nzfile[LRG]="/global/cfs/cdirs/desi/mocks/cai/holi/webjax_v4.80/nzref_da2_lrg.txt"
nzfile[ELG]="/global/cfs/cdirs/desi/mocks/cai/holi/webjax_v4.81/nzref_da2_elg_N.txt,/global/cfs/cdirs/desi/mocks/cai/holi/webjax_v4.81/nzref_da2_elg_S.txt"

tracers=('LRG','ELG','QSO')
SLURM_ARRAY_TASK_ID=500
seed=$(printf "seed%04d" "$SLURM_ARRAY_TASK_ID")
for tracer in "${tracers[@]}"; do

  version="${conf_version[$tracer]}"
  short_version="${sversion[$tracer]}"
  nzname="${nzfile[$tracer]}"
  echo $nzname
  echo $short_version
  echo $version
  echo $seed
  python $LSS/scripts/mock_tools/prepare_mocks_Y3_test1.py --survey DA2 --specdata loa-v1 --mockname holi --input_mockpath /global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/ --input_mockfile holi_"$tracer"_"$short_version"_GCcomb_clustering.dat.h5 --tracer $tracer --zrsdcol Z --output_fullpathfn /global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/"$tracer"/forFA0_Y3_noimagingmask_applied.fits --save_mock_nz n --nzfilename $nzname --need_nz_calib y
done
