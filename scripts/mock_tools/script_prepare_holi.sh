#!/bin/bash


declare -A conf_version

conf_version[QSO]="webjax_v4.80"
conf_version[ELG]="webjax_v4.80"
conf_version[LRG]="webjax_v4.80"


#echo $version


declare -A nzfile

nzfile[QSO]="/global/cfs/cdirs/desi/mocks/cai/holi/webjax_v4.80/nzref_da2_qso.txt"
nzfile[LRG]="/global/cfs/cdirs/desi/mocks/cai/holi/webjax_v4.80/nzref_da2_lrg.txt"
nzfile[ELG]="/global/cfs/cdirs/desi/mocks/cai/holi/webjax_v4.80/nzref_da2_elg_N.txt,/global/cfs/cdirs/desi/mocks/cai/holi/webjax_v4.80/nzref_da2_elg_S.txt"

#echo $version

tracer="ELG"

version="${conf_version[$tracer]}"
nzname="${nzfile[$tracer]}"

seed="seed0000"
srun -n 1 -C cpu -t 01:00:00 --cpus-per-task=64 --qos interactive --account desi python prepare_mocks_Y3_test1.py --survey DA2 --specdata loa-v1 --mockname holi --input_mockpath /global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/ --input_mockfile holi_ELG_v4.80_GCcomb_clustering.dat.h5 --tracer ELG --zrsdcol Z --output_fullpathfn /global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/"$tracer"/forFA0_Y3_noimagingmask_applied.fits --save_mock_nz y --nzfilename $nzname --need_nz_calib y

for i in {1..49}; do
    seed=$(printf "seed%04d" "$i")
if [ -f "/global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/$tracer/forFA0_Y3_noimagingmask_applied.fits" ]; then
    	echo "/global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/$tracer/forFA0_Y3_noimagingmask_applied.fits exists"
else
srun -n 1 -C cpu -t 01:00:00 --cpus-per-task=64 --qos interactive --account desi python prepare_mocks_Y3_test1.py --survey DA2 --specdata loa-v1 --mockname holi --input_mockpath /global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/ --input_mockfile holi_ELG_v4.80_GCcomb_clustering.dat.h5 --tracer ELG --zrsdcol Z --output_fullpathfn /global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/"$tracer"/forFA0_Y3_noimagingmask_applied.fits --save_mock_nz n --nzfilename $nzname --need_nz_calib y
		fi
done

tracer="LRG"
nzname="${nzfile[$tracer]}"
version="${conf_version[$tracer]}"

seed="seed0000"
srun -n 1 -C cpu -t 01:00:00 --cpus-per-task=64 --qos interactive --account desi python prepare_mocks_Y3_test1.py --survey DA2 --specdata loa-v1 --mockname holi --input_mockpath /global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/ --input_mockfile holi_"$tracer"_"$version"_GCcomb_clustering.dat.h5 --tracer $tracer --zrsdcol Z --output_fullpathfn /global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/"$tracer"/forFA0_Y3_noimagingmask_applied.fits --save_mock_nz y --nzfilename $nzname --need_nz_calib y

for i in {1..49}; do
    seed=$(printf "seed%04d" "$i")

#    if [ -f "/global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/holi_"$tracer"_"$version"_GCcomb_clustering.dat.h5" ]; then
##                echo "Processing $seed"
                if [ -f "/global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/$tracer/forFA0_Y3_noimagingmask_applied.fits" ]; then
        echo "File exists"
else
srun -n 1 -C cpu -t 01:00:00 --cpus-per-task=64 --qos interactive --account desi python prepare_mocks_Y3_test1.py --survey DA2 --specdata loa-v1 --mockname holi --input_mockpath /global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/ --input_mockfile holi_"$tracer"_"$version"_GCcomb_clustering.dat.h5 --tracer $tracer --zrsdcol Z --output_fullpathfn /global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/"$tracer"/forFA0_Y3_noimagingmask_applied.fits --save_mock_nz n --nzfilename $nzname --need_nz_calib y
                fi
#                fi
done

tracer="QSO"
nzname="${nzfile[$tracer]}"
version="${conf_version[$tracer]}"
seed="seed0000"

srun -n 1 -C cpu -t 01:00:00 --cpus-per-task=64 --qos interactive --account desi python prepare_mocks_Y3_test1.py --survey DA2 --specdata loa-v1 --mockname holi --input_mockpath /global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/ --input_mockfile holi_"$tracer"_"$version"_GCcomb_clustering.dat.h5 --tracer $tracer --zrsdcol Z --output_fullpathfn /global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/"$tracer"/forFA0_Y3_noimagingmask_applied.fits --save_mock_nz y --nzfilename $nzname --need_nz_calib y

for i in {1..49}; do
    seed=$(printf "seed%04d" "$i")

    #if [ -f "/global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/holi_"$tracer"_"$version"_GCcomb_clustering.dat.h5" ]; then
    #            echo "Processing $seed"
            if [ -f "/global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/$tracer/forFA0_Y3_noimagingmask_applied.fits" ]; then
    echo "File exists"
else
srun -n 1 -C cpu -t 01:00:00 --cpus-per-task=64 --qos interactive --account desi python prepare_mocks_Y3_test1.py --survey DA2 --specdata loa-v1 --mockname holi --input_mockpath /global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/ --input_mockfile holi_"$tracer"_"$version"_GCcomb_clustering.dat.h5 --tracer $tracer --zrsdcol Z --output_fullpathfn /global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/"$tracer"/forFA0_Y3_noimagingmask_applied.fits --save_mock_nz n --nzfilename $nzname --need_nz_calib y

		fi
#               fi
done

