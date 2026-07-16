#!/usr/bin/env bash

# Goal:
# ====
# Process holi pipeline for step 1 to 7 for one seed and all tracer ELG, LRG, QSO
# The possibility of using a RAM disk for file-based I/O data flow should be investigated.
# **script slurm free** for test, debug, optimization

# Loadbalancing:
# =============
# 1 simulation/seed/mock with 1 CPU
# use parallel submission with slurm

#
# script parameters
#
LSS_DIR=$1
DS_DIR=$2   # root directory of mock with version
IDS=$3      # id seed to process

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
        
time python ./prepare_mocks_Y3_test1.py --survey DA2 --specdata loa-v1 --mockname holi --input_mockpath $input_mockpath --input_mockfile $input_mockfile --tracer ELG --zrsdcol Z --output_fullpathfn $out1_ELG --save_mock_nz n --nzfilename $nzname --need_nz_calib y

tracer="LRG"
nzname="${nzfile[$tracer]}"
version="${conf_version[$tracer]}"
input_mockpath=/global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/
input_mockfile=holi_"$tracer"_v4.80_GCcomb_clustering.dat.h5
out1_LRG=$DS_DIR/$seed/"$tracer"/forFA0_Y3_noimagingmask_applied.fits

time python ./prepare_mocks_Y3_test1.py --survey DA2 --specdata loa-v1 --mockname holi --input_mockpath $input_mockpath --input_mockfile $input_mockfile --tracer LRG --zrsdcol Z --output_fullpathfn $out1_LRG --save_mock_nz n --nzfilename $nzname --need_nz_calib y

tracer="QSO"
nzname="${nzfile[$tracer]}"
version="${conf_version[$tracer]}"
input_mockpath=/global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/
input_mockfile=holi_"$tracer"_v4.80_GCcomb_clustering.dat.h5
out1_QSO=$DS_DIR/$seed/"$tracer"/forFA0_Y3_noimagingmask_applied.fits

time python ./prepare_mocks_Y3_test1.py --survey DA2 --specdata loa-v1 --mockname holi --input_mockpath $input_mockpath --input_mockfile $input_mockfile --tracer QSO --zrsdcol Z --output_fullpathfn $out1_QSO --save_mock_nz n --nzfilename $nzname --need_nz_calib y

#
# step 3 brickmask
#
echo "================= step 3"
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh dr1
#module load cpu cray-fftw
export CFITSIO_DIR=/global/common/software/desi/users/naimgk/cfitsio
export LD_LIBRARY_PATH=$CFITSIO_DIR/lib:$LD_LIBRARY_PATH
#TODO: path hardcoding, use file parameter
export EXE_PATH=/global/cfs/cdirs/desi/users/colley/wd_brick/brickmask
export CONF_PATH=/global/cfs/cdirs/desi/users/colley/LSS/scripts/mock_tools/pipeline

out3_ELG=$DS_DIR/$seed/ELG/$in_4
out3_LRG=$DS_DIR/$seed/LRG/$in_4
out3_QSO=$DS_DIR/$seed/QSO/$in_4
# BRICKMASK NOTE that command line options have priority over this file.
# Unnecessary entries can be left unset.
#srun -n 1 -c 3 --cpu-bind=cores $EXE_PATH/BRICKMASK -i $out1_ELG -o $out3_ELG -c $PROC_DIR/brickmask.conf

echo $out1_ELG > input$IDS.txt
echo $out3_ELG > output$IDS.txt
time $EXE_PATH/BRICKMASK -i input$IDS.txt -o output$IDS.txt -c $PROC_DIR/brickmask.conf

echo $out1_LRG > input$IDS.txt
echo $out3_LRG > output$IDS.txt
time $EXE_PATH/BRICKMASK -i input$IDS.txt -o output$IDS.txt -c $PROC_DIR/brickmask.conf

echo $out1_QSO > input$IDS.txt
echo $out3_QSO > output$IDS.txt
time $EXE_PATH/BRICKMASK -i input$IDS.txt -o output$IDS.txt -c $PROC_DIR/brickmask.conf

rm input$IDS.txt output$IDS.txt
#
# step 4 mask
#
echo "================= step 4"
out4_ELG=$DS_DIR/$seed/ELG/$in_5
out4_LRG=$DS_DIR/$seed/LRG/$in_5
out4_QSO=$DS_DIR/$seed/QSO/$in_5
time ./join_imaging_mask_stdpars.py --inputs $out3_ELG --outputs $out4_ELG 
time ./join_imaging_mask_stdpars.py --inputs $out3_LRG --outputs $out4_LRG 
time ./join_imaging_mask_stdpars.py --inputs $out3_QSO --outputs $out4_QSO 


#
# step 5 contaminant, **not for LRG**
#
echo "================= step 5"
date
out5_ELG=$DS_DIR/$seed/ELG/$in_6
out5_QSO=$DS_DIR/$seed/QSO/$in_6
time ./add_contaminants_to_mock_stdpars.py --inputs $out4_ELG --outputs $out5_ELG 
time ./add_contaminants_to_mock_stdpars.py --inputs $out4_QSO --outputs $out5_QSO 


#
# step 6 concatenate tracers
#
echo "================= step 6"
date
in6_ELG=$out5_ELG
in6_QSO=$out5_QSO
in6_LRG=$out4_LRG
out6=$DS_DIR/$(printf "forFA%04d.fits" "$IDS")
time ./concatenate_tracers_to_fba_stdpars.py --inputs $in6_ELG $in6_LRG $in6_QSO --outputs $out6
#time ./concatenate_tracers_to_fba_stdpars.py --inputs $in6_LRG $in6_LRG $in6_QSO --outputs $out6


#
# step 7 initialize amtl mocks
#
echo "================= step 7"
source /global/common/software/desi/desi_environment.sh main
module load desitarget/3.0.0

altmtlxxxx=$(printf "altmtl%04d" "$IDS")
out7=$DS_DIR/$altmtlxxxx
echo $out7
# use 1 CPU, may be memory duplication in multiproc mode
time ./initialize_amtl_mocks_da2_stdpars.py --inputs $out6 --outputs $out7 --obscon DARK

exit 0
