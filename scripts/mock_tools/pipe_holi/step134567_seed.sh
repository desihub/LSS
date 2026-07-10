#!/bin/bash

# 
# process holi pipeline for step 1 to 7 for one seed and all tracer ELG, LRG, QSO
# use 3 CPUs, one for earch tracer
# The possibility of using a RAM disk for file-based I/O data flow should be investigated.

#
# script parameters
#
LSS_DIR=$1
IDS=$2     # id seed to process
DS_DIR=$3 

#
# Parameters
#
PROC_DIR=$LSS_DIR/scripts/mock_tools/pipe_holi
seed=$(printf "seed%04d" "$IDS")
# dataflow name input
in_4=imforFA0_Y3_noimagingmask_applied.fits
in_5=forFA0.fits
in_6=forFA0_withcontaminants.fits
in_7=forFA0_concat.fits

###############################################
# Pipeline
###############################################

# -n or --ntasks 1 else slurm/perlmutter launch 255 task by default
# -c or --cpus-per-task 1 because python script are mono-process

#
# step 1 prepare mocks
#

# TODO: many option for step 1 create a file parameters
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
        
srun -n 1 -c 1 python $PROC_DIR/prepare_mocks_Y3_test1.py --survey DA2 --specdata loa-v1 --mockname holi --input_mockpath $input_mockpath --input_mockfile $input_mockfile --tracer ELG --zrsdcol Z --output_fullpathfn $out1_ELG --save_mock_nz n --nzfilename $nzname --need_nz_calib y &

tracer="LRG"
nzname="${nzfile[$tracer]}"
version="${conf_version[$tracer]}"
input_mockpath=/global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/
input_mockfile=holi_"$tracer"_v4.80_GCcomb_clustering.dat.h5
out1_LRG=$DS_DIR/$seed/"$tracer"/forFA0_Y3_noimagingmask_applied.fits

srun -n 1 -c 1 python $PROC_DIR/prepare_mocks_Y3_test1.py --survey DA2 --specdata loa-v1 --mockname holi --input_mockpath $input_mockpath --input_mockfile $input_mockfile --tracer LRG --zrsdcol Z --output_fullpathfn $out1_LRG --save_mock_nz n --nzfilename $nzname --need_nz_calib y &

tracer="QSO"
nzname="${nzfile[$tracer]}"
version="${conf_version[$tracer]}"
input_mockpath=/global/cfs/cdirs/desi/mocks/cai/holi/$version/$seed/
input_mockfile=holi_"$tracer"_v4.80_GCcomb_clustering.dat.h5
out1_QSO=$DS_DIR/$seed/"$tracer"/forFA0_Y3_noimagingmask_applied.fits

srun -n 1 -c 1 python $PROC_DIR/prepare_mocks_Y3_test1.py --survey DA2 --specdata loa-v1 --mockname holi --input_mockpath $input_mockpath --input_mockfile $input_mockfile --tracer QSO --zrsdcol Z --output_fullpathfn $out1_QSO --save_mock_nz n --nzfilename $nzname --need_nz_calib y &

wait

#
# step 3 brickmask
#
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh dr1
module load cpu cray-fftw
export CFITSIO_DIR=/global/common/software/desi/users/naimgk/cfitsio
export LD_LIBRARY_PATH=$CFITSIO_DIR/lib:$LD_LIBRARY_PATH
export EXE_PATH=/global/cfs/cdirs/desi/users/colley/brickmask
export CONF_PATH=/global/cfs/cdirs/desi/users/colley/LSS/scripts/mock_tools/pipeline

out3_ELG=$DS_DIR/$seed/ELG/$in_4
out3_LRG=$DS_DIR/$seed/LRG/$in_4
out3_QSO=$DS_DIR/$seed/QSO/$in_4
# BRICKMASK NOTE that command line options have priority over this file.
# Unnecessary entries can be left unset.
srun -n 1 -c 3 --cpu-bind=cores $EXE_PATH/BRICKMASK -i $out1_ELG -o $out3_ELG -c $PROC_DIR/brickmask.conf
srun -n 1 -c 3 --cpu-bind=cores $EXE_PATH/BRICKMASK -i $out1_LRG -o $out3_LRG -c $PROC_DIR/brickmask.conf
srun -n 1 -c 3 --cpu-bind=cores $EXE_PATH/BRICKMASK -i $out1_QSO -o $out3_QSO -c $PROC_DIR/brickmask.conf

#
# step 4 mask
#
out4_ELG=$DS_DIR/$seed/ELG/$in_5
out4_LRG=$DS_DIR/$seed/LRG/$in_5
out4_QSO=$DS_DIR/$seed/QSO/$in_5
srun -n 1 -c 1 $PROC_DIR/join_imaging_mask_stdpars.py --inputs $out3_ELG --outputs $out4_ELG &
srun -n 1 -c 1 $PROC_DIR/join_imaging_mask_stdpars.py --inputs $out3_LRG --outputs $out4_LRG &
srun -n 1 -c 1 $PROC_DIR/join_imaging_mask_stdpars.py --inputs $out3_QSO --outputs $out4_QSO &
wait 

#
# step 5 contaminant, **not for LRG**
#
echo "step 5"
date
out5_ELG=$DS_DIR/$seed/ELG/$in_6
out5_QSO=$DS_DIR/$seed/QSO/$in_6
srun -n 1 -c 1 $PROC_DIR/add_contaminants_to_mock_stdpars.py --inputs $out4_ELG --outputs $out5_ELG &
srun -n 1 -c 1 $PROC_DIR/add_contaminants_to_mock_stdpars.py --inputs $out4_QSO --outputs $out5_QSO &
wait

#
# step 6 concatenate tracers
#
echo "step 6"
date
in6_ELG=$out5_ELG
in6_QSO=$out5_QSO
in6_LRG=$out4_LRG
out6=$DS_DIR/$(printf "forFA%04d.fits" "$IDS")
srun -n 1 -c 1 $PROC_DIR/concatenate_tracers_to_fba_stdpars.py --inputs $in6_ELG $in6_LRG $in6_QSO --outputs $out6
date

#
# step 7 initialize amtl mocks
#
source /global/common/software/desi/desi_environment.sh main
module load desitarget/3.0.0

altmtlxxxx=$(printf "altmtl%04d" "$IDS")
out7=$DS_DIR/$altmtlxxxx
echo $out7
# use 1 CPU, may be memory duplication in multiproc mode
srun -n 1 -c 1 $PROC_DIR/initialize_amtl_mocks_da2_stdpars.py --inputs $out6 --outputs $out7 --obscon DARK

exit 0
