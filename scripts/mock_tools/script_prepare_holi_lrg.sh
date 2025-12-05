#!/bin/bash

for i in {202..300} {451..499} {601..650}; do
    seed=$(printf "seed%04d" "$i")
    echo "Processing $seed"
    if [ -f "/global/cfs/cdirs/desi/mocks/cai/holi/v4.00/$seed/QSO/forFA0_Y3_noimagingmask_applied.fits" ]; then
    	echo "File exists"
else
    python prepare_mocks_Y3_test1.py --survey DA2 --specdata loa-v1 --mockname holi --input_mockpath /global/cfs/cdirs/desi/mocks/cai/holi/v4.00/$seed/ --input_mockfile holi_QSO_v4.00_GCcomb_clustering.dat.h5 --tracer QSO --zrsdcol Z --output_fullpathfn /global/cfs/cdirs/desi/mocks/cai/holi/v4.00/$seed/QSO/forFA0_Y3_noimagingmask_applied.fits --save_mock_nz n --nzfilename /global/cfs/cdirs/desi/mocks/cai/holi/v4.00/nzref_da2_qso.txt --need_nz_calib y
    fi
done
