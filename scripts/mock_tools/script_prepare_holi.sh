#!/bin/bash

declare -A conf_version

conf_version[QSO]="v4.00"
conf_version[ELG]="v5.0_Y5"
conf_version[LRG]="v4.00"

tracer="ELG"

version="${conf_version[$tracer]}"
echo $version


declare -A nzfile

nzfile[QSO]="/global/cfs/cdirs/desi/mocks/cai/holi/v4.00/nzref_da2_qso.txt"
nzfile[LRG]="/global/cfs/cdirs/desi/mocks/cai/holi/v4.00/nzref_da2_lrg.txt"
nzfile[ELG]="/global/cfs/cdirs/desi/mocks/cai/holi/v5.0_Y5/nzref_da2_elg_N.txt,/global/cfs/cdirs/desi/mocks/cai/holi/v5.0_Y5/nzref_da2_elg_S.txt"

echo $version




for i in {0..1000}; do
    seed=$(printf "seed%04d" "$i")

    if [ -f "/global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/holi_"$tracer"_"$version"_GCcomb_clustering.dat.h5" ]; then
	   	echo "Processing $seed"
		if [ -f "/global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/$tracer/forFA0_Y3_noimagingmask_applied.fits" ]; then
    	echo "File exists"
else
	python prepare_mocks_Y3_test1.py --survey DA2 --specdata loa-v1 --mockname holi --input_mockpath /global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/ --input_mockfile holi_"$tracer"_v5.0_GCcomb_clustering.dat.h5 --tracer $tracer --zrsdcol Z --output_fullpathfn /global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/"$tracer"/forFA0_Y3_noimagingmask_applied.fits --save_mock_nz n --nzfilename $nzname --need_nz_calib y
		fi
	    	fi
done

tracer="LRG"
nzname="${nzfile[$tracer]}"
version="${conf_version[$tracer]}"

for i in {0..1000}; do
    seed=$(printf "seed%04d" "$i")

    if [ -f "/global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/holi_"$tracer"_"$version"_GCcomb_clustering.dat.h5" ]; then
                echo "Processing $seed"
                if [ -f "/global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/$tracer/forFA0_Y3_noimagingmask_applied.fits" ]; then
        echo "File exists"
else
        python prepare_mocks_Y3_test1.py --survey DA2 --specdata loa-v1 --mockname holi --input_mockpath /global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/ --input_mockfile holi_"$tracer"_"$version"_GCcomb_clustering.dat.h5 --tracer $tracer --zrsdcol Z --output_fullpathfn /global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/"$tracer"/forFA0_Y3_noimagingmask_applied.fits --save_mock_nz n --nzfilename $nzname --need_nz_calib y
                fi
                fi
done

tracer="QSO"
nzname="${nzfile[$tracer]}"
version="${conf_version[$tracer]}"

for i in {0..1000}; do
    seed=$(printf "seed%04d" "$i")

    if [ -f "/global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/holi_"$tracer"_"$version"_GCcomb_clustering.dat.h5" ]; then
                echo "Processing $seed"
                if [ -f "/global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/$tracer/forFA0_Y3_noimagingmask_applied.fits" ]; then
        echo "File exists"
else
	python prepare_mocks_Y3_test1.py --survey DA2 --specdata loa-v1 --mockname holi --input_mockpath /global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/ --input_mockfile holi_"$tracer"_"$version"_GCcomb_clustering.dat.h5 --tracer $tracer --zrsdcol Z --output_fullpathfn /global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/"$tracer"/forFA0_Y3_noimagingmask_applied.fits --save_mock_nz n --nzfilename $nzname --need_nz_calib y

		fi
                fi
done


        ##srun -N 1 -C cpu -t 01:00:00 --qos interactive --account desi python prepare_mocks_Y3_test1.py --survey DA2 --specdata loa-v1 --mockname holi --input_mockpath /global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/ --input_mockfile holi_"$tracer"_"$version"_GCcomb_clustering.dat.h5 --tracer $tracer --zrsdcol Z --output_fullpathfn /global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/"$tracer"/forFA0_Y3_noimagingmask_applied.fits --save_mock_nz n --nzfilename $nzname --need_nz_calib y
