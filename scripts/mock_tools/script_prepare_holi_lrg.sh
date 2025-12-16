#!/bin/bash

for i in {451..500} {601..650}; do
#for i in {202..300} {451..499} {601..650}; do
    seed=$(printf "seed%04d" "$i")
    echo "Processing $seed"
    if [ -f "/global/cfs/cdirs/desi/mocks/cai/holi/v5.0/$seed/ELG/forFA0_Y3_noimagingmask_applied.fits" ]; then
    #if [ -f "/global/cfs/cdirs/desi/mocks/cai/holi/v4.00/$seed/QSO/forFA0_Y3_noimagingmask_applied.fits" ]; then
    	echo "File exists"
else
	srun -N 1 -C cpu -t 02:00:00 --qos interactive --account desi python prepare_mocks_Y3_test1.py --survey DA2 --specdata loa-v1 --mockname holi --input_mockpath /global/cfs/cdirs/desi/mocks/cai/holi/v5.0/$seed --input_mockfile holi_ELG_v5.0_GCcomb_clustering.dat.h5 --tracer ELG --zrsdcol Z --output_fullpathfn /global/cfs/cdirs/desi/mocks/cai/holi/v5.0/$seed/ELG/forFA0_Y3_noimagingmask_applied.fits --save_mock_nz n --nzfilename /global/cfs/cdirs/desi/mocks/cai/holi/v5.0/nzref_da2_elg_N.txt,/global/cfs/cdirs/desi/mocks/cai/holi/v5.0/nzref_da2_elg_S.txt --need_nz_calib y

    #python prepare_mocks_Y3_test1.py --survey DA2 --specdata loa-v1 --mockname holi --input_mockpath /global/cfs/cdirs/desi/mocks/cai/holi/v4.00/$seed/ --input_mockfile holi_QSO_v4.00_GCcomb_clustering.dat.h5 --tracer QSO --zrsdcol Z --output_fullpathfn /global/cfs/cdirs/desi/mocks/cai/holi/v4.00/$seed/QSO/forFA0_Y3_noimagingmask_applied.fits --save_mock_nz n --nzfilename /global/cfs/cdirs/desi/mocks/cai/holi/v4.00/nzref_da2_qso.txt --need_nz_calib y
    fi
done
