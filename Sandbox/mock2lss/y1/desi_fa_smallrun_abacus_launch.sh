#!/bin/bash

OUTDIR=$1    # /global/cfs/cdirs/desi/users/raichoor/fiberassign-abacus1year-20211116/
PROGRAM=$2   # bright or dark
FOOTPRINT=$3 # ngc or sgc or all 
STEPS=$4     # tiles,targ or fa
NUMPROC=$5   # 1 if tiles,targ or stats; 32 if fa


source /global/cfs/cdirs/desi/software/desi_environment.sh master
module swap desimodel/master
module swap fiberassign/5.3.0

#
#INFN= /global/project/projectdirs/desi/users/mvargas/for_EZmocks/mtlz_ABACUS-000_alltracers_1year.fits
TILESFN=/project/projectdirs/desi/users/mvargas/for_EZmocks/Tiles_year1_formated.fits

RUNDATE=2021-12-04T00:00:00+00:00

SURVEY=main
DTVER=1.1.1

if [[ "$FOOTPRINT" == "ngc" ]]
then
    RADEC=85,300,-90,90
fi
if [[ "$FOOTPRINT" == "sgc" ]]
then
    RADEC=300,85,-90,90
fi
if [[ "$FOOTPRINT" == "all" ]]
then
    RADEC=0,360,-90,90
fi


if [[ "$PROGRAM" == "bright" ]]
then
    NPASS=4
fi
if [[ "$PROGRAM" == "dark" ]]
then
    NPASS=7
fi


#MYOUTDIR=$OUTDIR/$SURVEY-$DTVER-$PROGRAM-$FOOTPRINT
MYOUTDIR=$OUTDIR/

CMD="./desi_fa_smallrun_abacus --infn $6 --outdir $MYOUTDIR --radec $RADEC --rundate $RUNDATE --survey $SURVEY --program $PROGRAM --npass $NPASS --dtver $DTVER --tilesfn $TILESFN --numproc $NUMPROC --steps $STEPS"

#CMD="./desi_fa_smallrun_abacus --infn $INFN --outdir $MYOUTDIR --radec $RADEC --rundate $RUNDATE --survey $SURVEY --program $PROGRAM --npass $NPASS --dtver $DTVER --tilesfn $TILESFN --numproc $NUMPROC --steps $STEPS"
echo $CMD
eval $CMD
