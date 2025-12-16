#!/bin/bash

declare -A conf_version

conf_version[QSO]="v4.00"
conf_version[ELG]="v5.0"
conf_version[LRG]="v4.00"

tracer="QSO"

version="${conf_version[$tracer]}"
echo $version


declare -A nzfile

nzfile[QSO]="/global/cfs/cdirs/desi/mocks/cai/holi/v4.00/nzref_da2_qso.txt"
nzfile[LRG]="/global/cfs/cdirs/desi/mocks/cai/holi/v4.00/nzref_da2_lrg.txt"
nzfile[ELG]="/global/cfs/cdirs/desi/mocks/cai/holi/v5.0/nzref_da2_elg_N.txt,/global/cfs/cdirs/desi/mocks/cai/holi/v5.0/nzref_da2_elg_S.txt"

nzname="${nzfile[$tracer]}"


for i in {0..1000}; do
    seed=$(printf "seed%04d" "$i")

    if [ -f "/global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/holi_"$tracer"_"$version"_GCcomb_clustering.dat.h5" ]; then
	    echo "Processing $seed"
    		#if [ -f "/global/cfs/cdirs/desi/mocks/cai/holi/v4.00/$seed/QSO/forFA0_Y3_noimagingmask_applied.fits" ]; then
#    	echo "File exists"
#else
	srun -N 1 -C cpu -t 01:00:00 --qos interactive --account desi python prepare_mocks_Y3_test1.py --survey DA2 --specdata loa-v1 --mockname holi --input_mockpath /global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/ --input_mockfile holi_"$tracer"_"$version"_GCcomb_clustering.dat.h5 --tracer $tracer --zrsdcol Z --output_fullpathfn /global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/"$tracer"_v2/forFA0_Y3_noimagingmask_applied.fits --save_mock_nz n --nzfilename $nzname --need_nz_calib y
    fi
done
