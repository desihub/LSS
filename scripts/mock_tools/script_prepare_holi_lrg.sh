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

tracer="LRG"

version="${conf_version[$tracer]}"
nzname="${nzfile[$tracer]}"

seed="seed0001"
python prepare_mocks_Y3_test1.py --survey DA2 --specdata loa-v1 --mockname holi --input_mockpath /global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/ --input_mockfile holi_LRG_v4.80_GCcomb_clustering.dat.h5 --tracer LRG --zrsdcol Z --output_fullpathfn /global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/"$tracer"/forFA0_Y3_noimagingmask_applied.fits --save_mock_nz y --nzfilename $nzname --need_nz_calib y

for i in 2 3 4 5 6 7 8 9 10 11 25 26 54 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 81 82 83 120 121 126 127 128 129 130 133 134 135 136 139 140 141 183 184 185 186 187 201 202 240 241 242 248 249 271 280 281 286 287 288 289 290 293 296 297 323 324 334 335 336 337 338 339 340 341 342 343 344 372 373 380 381 382 383 384 387 389 390 391 400 434 435 436 437 438 439 440 441 442 476 477 484 485 486 487 496 532 533 535 566 567 568 569 570 571 572 592 624 625 626 648 649 650 651 652 665 666 673 718 719 720 721 722 744 745 758 759 760 761 765 766 767 784 787 803 804 812 813 814 817 830 842 843 859 860 861 876 877 878 889 890 891 892 893 894 906 907 908 910 927 928 929 930 931 932 933 934 935 936 939 953 954 965 966 967 968 969 970 971 972 973 974 975 976 977 978 979; do
    seed=$(printf "seed%04d" "$i")
#    if [ -f "/global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/holi_"$tracer"_"$version"_GCcomb_clustering.dat.h5" ]; then
#	   	echo "Processing $seed"
		if [ -f "/global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/$tracer/forFA0_Y3_noimagingmask_applied.fits" ]; then
    	echo "/global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/$tracer/forFA0_Y3_noimagingmask_applied.fits exists"
else
	##nzname="/global/cfs/cdirs/desi/mocks/cai/holi/webjax_v4.80/nzref_da2_lrg_$seed.txt"
	python prepare_mocks_Y3_test1.py --survey DA2 --specdata loa-v1 --mockname holi --input_mockpath /global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/ --input_mockfile holi_LRG_v4.80_GCcomb_clustering.dat.h5 --tracer LRG --zrsdcol Z --output_fullpathfn /global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/"$tracer"/forFA0_Y3_noimagingmask_applied.fits --save_mock_nz n --nzfilename $nzname --need_nz_calib y
		fi
#	    	fi
done

#tracer="LRG"
#nzname="${nzfile[$tracer]}"
#version="${conf_version[$tracer]}"

#for i in {0..1000}; do
#    seed=$(printf "seed%04d" "$i")

#    if [ -f "/global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/holi_"$tracer"_"$version"_GCcomb_clustering.dat.h5" ]; then
#                echo "Processing $seed"
#                if [ -f "/global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/$tracer/forFA0_Y3_noimagingmask_applied.fits" ]; then
#        echo "File exists"
#else
#        python prepare_mocks_Y3_test1.py --survey DA2 --specdata loa-v1 --mockname holi --input_mockpath /global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/ --input_mockfile holi_"$tracer"_"$version"_GCcomb_clustering.dat.h5 --tracer $tracer --zrsdcol Z --output_fullpathfn /global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/"$tracer"/forFA0_Y3_noimagingmask_applied.fits --save_mock_nz n --nzfilename $nzname --need_nz_calib y
#                fi
#                fi
#done

#tracer="QSO"
#nzname="${nzfile[$tracer]}"
#version="${conf_version[$tracer]}"

#for i in {0..15}; do
#    seed=$(printf "seed%04d" "$i")

    #if [ -f "/global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/holi_"$tracer"_"$version"_GCcomb_clustering.dat.h5" ]; then
    #            echo "Processing $seed"
    #            if [ -f "/global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/$tracer/forFA0_Y3_noimagingmask_applied.fits" ]; then
    #    echo "File exists"
#else
#	python prepare_mocks_Y3_test1.py --survey DA2 --specdata loa-v1 --mockname holi --input_mockpath /global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/ --input_mockfile holi_"$tracer"_"$version"_GCcomb_clustering.dat.h5 --tracer $tracer --zrsdcol Z --output_fullpathfn /global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/"$tracer"/forFA0_Y3_noimagingmask_applied.fits --save_mock_nz n --nzfilename $nzname --need_nz_calib y

#		fi
#               fi
#done
#------------------------------------------
#for i in {21..29}; do
#    seed=$(printf "seed%04d" "$i")

    #if [ -f "/global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/holi_"$tracer"_"$version"_GCcomb_clustering.dat.h5" ]; then
    #            echo "Processing $seed"
    #            if [ -f "/global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/$tracer/forFA0_Y3_noimagingmask_applied.fits" ]; then
    #    echo "File exists"
#else
#	python prepare_mocks_Y3_test1.py --survey DA2 --specdata loa-v1 --mockname holi --input_mockpath /global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/ --input_mockfile holi_"$tracer"_"$version"_GCcomb_clustering.dat.h5 --tracer $tracer --zrsdcol Z --output_fullpathfn /global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/"$tracer"/forFA0_Y3_noimagingmask_applied.fits --save_mock_nz n --nzfilename $nzname --need_nz_calib y

#		fi
#               fi
#done

        ##srun -N 1 -C cpu -t 01:00:00 --qos interactive --account desi python prepare_mocks_Y3_test1.py --survey DA2 --specdata loa-v1 --mockname holi --input_mockpath /global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/ --input_mockfile holi_"$tracer"_"$version"_GCcomb_clustering.dat.h5 --tracer $tracer --zrsdcol Z --output_fullpathfn /global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/"$tracer"/forFA0_Y3_noimagingmask_applied.fits --save_mock_nz n --nzfilename $nzname --need_nz_calib y
