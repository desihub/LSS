#!/usr/bin/env bash

# Goal:
# ====

# Loadbalancing:
# =============

echo "$(date), Step 1"

#
# script parameters
#
LSS_DIR=$1
DS_DIR=$2           # root directory of mock with version
FIRST_ID_RANK=$3    # first id seed in rank job array

#
# SLURM parameters
#
# task rank
PROCID=${SLURM_PROCID:-0}
# ID seed mock to process in task rank
IDS=$((FIRST_ID_RANK + PROCID))

# NCPU=${SLURM_CPUS_PER_TASK:-1}
# NCPU_M2=$((NCPU-2))

#
# Internal variables
#
seed=$(printf "seed%04d" "$IDS")
# dataflow name input
in_4=imforFA0_Y3_noimagingmask_applied.fits
in_5=forFA0.fits
in_6=forFA0_withcontaminants.fits
in_7=forFA0_concat.fits

#
# step 1 prepare mocks
#
echo "================= step 1"
## Env
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
export PYTHONPATH=$LSS_DIR/py:$LSS_DIR/scripts/mock_tools:$PYTHONPATH
export PATH=$LSS_DIR/bin:$PATH
PROC_DIR=$LSS_DIR/scripts/mock_tools/holi_pipeline
cd $PROC_DIR

## TEST/DEBUG: negative value for no debug
max_gal=$(get_pars.py $HOLI_PARS max_gal)

## simulation
input_ref=$(get_pars.py $HOLI_PARS input_ref)
SURV=$(get_pars.py $HOLI_PARS prepare_mocks.survey)

declare -A nzfile
nzfile[QSO]="$DS_DIR/nzref_da2_qso.txt"
nzfile[LRG]="$DS_DIR/nzref_da2_lrg.txt"
nzfile[ELG]="$DS_DIR/nzref_da2_elg_N.txt,$DS_DIR/nzref_da2_elg_S.txt"

#
tracer="ELG"
#
version="${conf_version[$tracer]}"
nzname="${nzfile[$tracer]}"
input_mockpath=$input_ref/$seed/
input_mockfile=holi_"$tracer"_v4.80_GCcomb_clustering.dat.h5
out1_ELG=$DS_DIR/$seed/"$tracer"/forFA0_Y3_noimagingmask_applied.fits
time python ./prepare_mocks_Y3.py --limit_for_test $max_gal --survey $SURV --specdata loa-v1 --mockname holi --input_mockpath $input_mockpath --input_mockfile $input_mockfile --tracer ELG --zrsdcol Z --output_fullpathfn $out1_ELG --save_mock_nz n --nzfilename $nzname --need_nz_calib y  
status_elg=$?

#
tracer="LRG"
#
nzname="${nzfile[$tracer]}"
version="${conf_version[$tracer]}"
input_mockpath=$input_ref/$seed/
input_mockfile=holi_"$tracer"_v4.80_GCcomb_clustering.dat.h5
out1_LRG=$DS_DIR/$seed/"$tracer"/forFA0_Y3_noimagingmask_applied.fits
time python ./prepare_mocks_Y3.py --limit_for_test $max_gal --survey $SURV --specdata loa-v1 --mockname holi --input_mockpath $input_mockpath --input_mockfile $input_mockfile --tracer LRG --zrsdcol Z --output_fullpathfn $out1_LRG --save_mock_nz n --nzfilename $nzname --need_nz_calib y 
status_lrg=$?

#
tracer="QSO"
#
nzname="${nzfile[$tracer]}"
version="${conf_version[$tracer]}"
input_mockpath=$input_ref/$seed/
input_mockfile=holi_"$tracer"_v4.80_GCcomb_clustering.dat.h5
out1_QSO=$DS_DIR/$seed/"$tracer"/forFA0_Y3_noimagingmask_applied.fits
time python ./prepare_mocks_Y3.py --limit_for_test $max_gal --survey $SURV --specdata loa-v1 --mockname holi --input_mockpath $input_mockpath --input_mockfile $input_mockfile --tracer QSO --zrsdcol Z --output_fullpathfn $out1_QSO --save_mock_nz n --nzfilename $nzname --need_nz_calib y 
status_qso=$?

#
# concatenate input / output of current seed
#
out3_ELG=$DS_DIR/$seed/ELG/$in_4
out3_LRG=$DS_DIR/$seed/LRG/$in_4
out3_QSO=$DS_DIR/$seed/QSO/$in_4

if [[ $status_elg -eq 0 && $status_lrg -eq 0 && $status_qso -eq 0 ]]; then
	printf '%s\n' "$out1_ELG" "$out1_LRG" "$out1_QSO" > "$DS_DIR/input$IDS.txt"
	printf '%s\n' "$out3_ELG" "$out3_LRG" "$out3_QSO" > "$DS_DIR/output$IDS.txt"
    exit 0
else
	: > "$DS_DIR/input$IDS.txt"
	: > "$DS_DIR/output$IDS.txt"
    exit 1
fi
