#!/usr/bin/env bash

# Goal:
# ====

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
IDS=$((FIRST_ID_RANK + PROCID ))

# NCPU=${SLURM_CPUS_PER_TASK:-4}
# NCPU_M2=$((NCPU-2))

#
# internal variables
#
HOLI_DIR=$LSS_DIR/scripts/mock_tools/holi_pipeline
cd $HOLI_DIR
seed=$(printf "seed%04d" "$IDS")
# dataflow name input
in_4=imforFA0_Y3_noimagingmask_applied.fits
in_5=forFA0.fits
in_6=forFA0_withcontaminants.fits
in_7=forFA0_concat.fits

#
# step 4 mask
#
echo "================= step 4"
## Env
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh dr1
export PYTHONPATH=$LSS_DIR/py:$PYTHONPATH
export PATH=$LSS_DIR/bin:$PATH
## in4=out3, out4
out3_ELG=$DS_DIR/$seed/ELG/$in_4
out3_LRG=$DS_DIR/$seed/LRG/$in_4
out3_QSO=$DS_DIR/$seed/QSO/$in_4
out4_ELG=$DS_DIR/$seed/ELG/$in_5
out4_LRG=$DS_DIR/$seed/LRG/$in_5
out4_QSO=$DS_DIR/$seed/QSO/$in_5
## join tracers
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
time ./concatenate_tracers_to_fba_stdpars.py --inputs $in6_ELG $in6_LRG $in6_QSO --outputs $out6 --id_seed $IDS

#
# step 7 initialize amtl mocks
#
echo "================= step 7"
date
## Env
source /global/common/software/desi/desi_environment.sh main
module load desitarget/3.0.0
export PYTHONPATH=$LSS_DIR/py:$PYTHONPATH
export PATH=$LSS_DIR/bin:$PATH

altmtlxxxx=$(printf "altmtl%04d" "$IDS")
out7=$DS_DIR/$altmtlxxxx
echo $out7

time ./initialize_amtl_mocks_da2_stdpars.py --inputs $out6 --outputs $out7 --obscon DARK

date
exit 0
