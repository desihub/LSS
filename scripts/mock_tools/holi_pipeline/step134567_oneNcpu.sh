#!/usr/bin/env bash

# Goal:
# ====
# Process holi pipeline for step 1 to 7 for one seed and all tracer ELG, LRG, QSO
# The possibility of using a RAM disk for file-based I/O data flow should be investigated.
#   - certainly useful for BRICKMASK which has intensive I/O, but
#     huge file in / out (same file with only 3 colmuns added) 
#      => can need all memory of the node
#      => may be possible with multi-node job ...
#      => first implementation, some processing ok for only one seed but almost kill process ... 

# Loadbalancing:
# =============
# 1 simulation/seed/mock with N CPUs, where N >= 4 
# use parallel submission with slurm
date
#
# script parameters
#
LSS_DIR=$1
DS_DIR=$2   # root directory of mock with version
IDS=$3      # id seed to process
NCPU=${SLURM_CPUS_PER_TASK:-4}
NCPU_M2=$((NCPU-2))

# Node-local temporary space; fallback is /tmp if /dev/shm is unavailable.
#RAMDISK=${RAMDISK:-/dev/shm/$USER/${IDS:-$$}}
#mkdir -p "$RAMDISK"

#
# internal variables
#
PROC_DIR=$LSS_DIR/scripts/mock_tools/holi_pipeline
cd $PROC_DIR
ls -trl
seed=$(printf "seed%04d" "$IDS")
# dataflow name input
in_4=imforFA0_Y3_noimagingmask_applied.fits
in_5=forFA0.fits
in_6=forFA0_withcontaminants.fits
in_7=forFA0_concat.fits

###############################################
# Holi Pipeline step 1 to 7
###############################################

#
# step 1 prepare mocks
#
echo "================= step 1"
# TODO: manage version ?
# TODO: use file parameters 
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main

export PYTHONPATH=$LSS_DIR/py:$PYTHONPATH
export PATH=$LSS_DIR/bin:$PATH

declare -A conf_version
conf_version[QSO]="webjax_v4.80"
conf_version[ELG]="webjax_v4.80"
conf_version[LRG]="webjax_v4.80"

declare -A nzfile
nzfile[QSO]="$DS_DIR/nzref_da2_qso.txt"
nzfile[LRG]="$DS_DIR/nzref_da2_lrg.txt"
nzfile[ELG]="$DS_DIR/nzref_da2_elg_N.txt,$DS_DIR/nzref_da2_elg_S.txt"

tracer="ELG"
version="${conf_version[$tracer]}"
nzname="${nzfile[$tracer]}"
input_mockpath=/global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/
input_mockfile=holi_"$tracer"_v4.80_GCcomb_clustering.dat.h5
out1_ELG=$DS_DIR/$seed/"$tracer"/forFA0_Y3_noimagingmask_applied.fits
#out1_ELG=$RAMDISK/$seed/"$tracer"/forFA0_Y3_noimagingmask_applied.fits

max_gal=-1
        
#time srun -n 1 -c $NCPU_M2 --exclusive python ./prepare_mocks_Y3_test1.py --limit_for_test $max_gal --survey DA2 --specdata loa-v1 --mockname holi --input_mockpath $input_mockpath --input_mockfile $input_mockfile --tracer ELG --zrsdcol Z --output_fullpathfn $out1_ELG --save_mock_nz n --nzfilename $nzname --need_nz_calib y  &
time python ./prepare_mocks_Y3_test1.py --limit_for_test $max_gal --survey DA2 --specdata loa-v1 --mockname holi --input_mockpath $input_mockpath --input_mockfile $input_mockfile --tracer ELG --zrsdcol Z --output_fullpathfn $out1_ELG --save_mock_nz n --nzfilename $nzname --need_nz_calib y  &


tracer="LRG"
nzname="${nzfile[$tracer]}"
version="${conf_version[$tracer]}"
input_mockpath=/global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/
input_mockfile=holi_"$tracer"_v4.80_GCcomb_clustering.dat.h5
out1_LRG=$DS_DIR/$seed/"$tracer"/forFA0_Y3_noimagingmask_applied.fits
#out1_LRG=$RAMDISK/$seed/"$tracer"/forFA0_Y3_noimagingmask_applied.fits
#time python ./prepare_mocks_Y3_test1.py --survey DA2 --specdata loa-v1 --mockname holi --input_mockpath $input_mockpath --input_mockfile $input_mockfile --tracer QSO --zrsdcol Z --output_fullpathfn $out1_QSO --save_mock_nz n --nzfilename $nzname --need_nz_calib y &
#time srun -n 1 -c 1 --exclusive time python ./prepare_mocks_Y3_test1.py --limit_for_test $max_gal --survey DA2 --specdata loa-v1 --mockname holi --input_mockpath $input_mockpath --input_mockfile $input_mockfile --tracer LRG --zrsdcol Z --output_fullpathfn $out1_LRG --save_mock_nz n --nzfilename $nzname --need_nz_calib y &
time python ./prepare_mocks_Y3_test1.py --limit_for_test $max_gal --survey DA2 --specdata loa-v1 --mockname holi --input_mockpath $input_mockpath --input_mockfile $input_mockfile --tracer LRG --zrsdcol Z --output_fullpathfn $out1_LRG --save_mock_nz n --nzfilename $nzname --need_nz_calib y &



tracer="QSO"
nzname="${nzfile[$tracer]}"
version="${conf_version[$tracer]}"
input_mockpath=/global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/
input_mockfile=holi_"$tracer"_v4.80_GCcomb_clustering.dat.h5
out1_QSO=$DS_DIR/$seed/"$tracer"/forFA0_Y3_noimagingmask_applied.fits
#out1_QSO=$RAMDISK/$seed/"$tracer"/forFA0_Y3_noimagingmask_applied.fits

#time srun -n 1 -c 1 --exclusive python ./prepare_mocks_Y3_test1.py --limit_for_test $max_gal --survey DA2 --specdata loa-v1 --mockname holi --input_mockpath $input_mockpath --input_mockfile $input_mockfile --tracer QSO --zrsdcol Z --output_fullpathfn $out1_QSO --save_mock_nz n --nzfilename $nzname --need_nz_calib y &
time python ./prepare_mocks_Y3_test1.py --limit_for_test $max_gal --survey DA2 --specdata loa-v1 --mockname holi --input_mockpath $input_mockpath --input_mockfile $input_mockfile --tracer QSO --zrsdcol Z --output_fullpathfn $out1_QSO --save_mock_nz n --nzfilename $nzname --need_nz_calib y &


wait


#
# step 3 brickmask
#
echo "================= step 3"
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh dr1
export PYTHONPATH=$LSS_DIR/py:$PYTHONPATH
export PATH=$LSS_DIR/bin:$PATH

#module load cpu cray-fftw
export CFITSIO_DIR=/global/common/software/desi/users/naimgk/cfitsio
export LD_LIBRARY_PATH=$CFITSIO_DIR/lib:$LD_LIBRARY_PATH
#TODO: path hardcoding, use file parameter
export EXE_PATH=/global/cfs/cdirs/desi/users/colley/wd_brick/brickmask_mpi
export CONF_PATH=/global/cfs/cdirs/desi/users/colley/LSS/scripts/mock_tools/pipeline

out3_ELG=$DS_DIR/$seed/ELG/$in_4
#out3_ELG=$RAMDISK/$seed/ELG/$in_4
out3_LRG=$DS_DIR/$seed/LRG/$in_4
#out3_LRG=$RAMDISK/$seed/LRG/$in_4
out3_QSO=$DS_DIR/$seed/QSO/$in_4
#out3_QSO=$RAMDISK/$seed/QSO/$in_4
# BRICKMASK NOTE that command line options have priority over this file.
# Unnecessary entries can be left unset.
#srun -n 1 -c 3 --cpu-bind=cores $EXE_PATH/BRICKMASK -i $out1_ELG -o $out3_ELG -c $PROC_DIR/brickmask.conf


# merge input output for one call to BRICKMASK 
echo $out1_ELG > $DS_DIR/input$IDS.txt
echo $out3_ELG > $DS_DIR/output$IDS.txt
echo $out1_LRG >> $DS_DIR/input$IDS.txt
echo $out3_LRG >> $DS_DIR/output$IDS.txt
echo $out1_QSO >> $DS_DIR/input$IDS.txt
echo $out3_QSO >> $DS_DIR/output$IDS.txt

time srun --overlap -n $NCPU -c 1 --cpu-bind=cores $EXE_PATH/BRICKMASK -i $DS_DIR/input$IDS.txt -o $DS_DIR/output$IDS.txt -c $PROC_DIR/brickmask.conf

rm $DS_DIR/input$IDS.txt $DS_DIR/output$IDS.txt
#rm -rf $RAMDISK


#
# step 4 mask
#
echo "================= step 4"
out4_ELG=$DS_DIR/$seed/ELG/$in_5
out4_LRG=$DS_DIR/$seed/LRG/$in_5
out4_QSO=$DS_DIR/$seed/QSO/$in_5
# time srun -n 1 -c $NCPU_M2 --exclusive ./join_imaging_mask_stdpars.py --inputs $out3_ELG --outputs $out4_ELG &
# time srun -n 1 -c 1 --exclusive ./join_imaging_mask_stdpars.py --inputs $out3_LRG --outputs $out4_LRG &
# time srun -n 1 -c 1 --exclusive ./join_imaging_mask_stdpars.py --inputs $out3_QSO --outputs $out4_QSO &
time ./join_imaging_mask_stdpars.py --inputs $out3_ELG --outputs $out4_ELG &
time ./join_imaging_mask_stdpars.py --inputs $out3_LRG --outputs $out4_LRG &
time ./join_imaging_mask_stdpars.py --inputs $out3_QSO --outputs $out4_QSO &
wait

#
# step 5 contaminant, **not for LRG**
#
echo "================= step 5"
date
out5_ELG=$DS_DIR/$seed/ELG/$in_6
out5_QSO=$DS_DIR/$seed/QSO/$in_6
# time srun -n 1 -c $NCPU_M2 --exclusive ./add_contaminants_to_mock_stdpars.py --inputs $out4_ELG --outputs $out5_ELG &
# time srun -n 1 -c 2 --exclusive ./add_contaminants_to_mock_stdpars.py --inputs $out4_QSO --outputs $out5_QSO &
time ./add_contaminants_to_mock_stdpars.py --inputs $out4_ELG --outputs $out5_ELG &
time ./add_contaminants_to_mock_stdpars.py --inputs $out4_QSO --outputs $out5_QSO &
wait

#
# step 6 concatenate tracers
#
echo "================= step 6"
date
in6_ELG=$out5_ELG
in6_QSO=$out5_QSO
in6_LRG=$out4_LRG
out6=$DS_DIR/$(printf "forFA%04d.fits" "$IDS")
#time srun -n 1 -c $NCPU --exclusive ./concatenate_tracers_to_fba_stdpars.py --inputs $in6_ELG $in6_LRG $in6_QSO --outputs $out6 --id_seed $IDS
time ./concatenate_tracers_to_fba_stdpars.py --inputs $in6_ELG $in6_LRG $in6_QSO --outputs $out6 --id_seed $IDS



#
# step 7 initialize amtl mocks
#
echo "================= step 7"
source /global/common/software/desi/desi_environment.sh main
module load desitarget/3.0.0
export PYTHONPATH=$LSS_DIR/py:$PYTHONPATH
export PATH=$LSS_DIR/bin:$PATH

altmtlxxxx=$(printf "altmtl%04d" "$IDS")
out7=$DS_DIR/$altmtlxxxx
echo $out7

#time srun -n 1 -c $NCPU --exclusive ./initialize_amtl_mocks_da2_stdpars.py --inputs $out6 --outputs $out7 --obscon DARK
time ./initialize_amtl_mocks_da2_stdpars.py --inputs $out6 --outputs $out7 --obscon DARK
date
exit 0
