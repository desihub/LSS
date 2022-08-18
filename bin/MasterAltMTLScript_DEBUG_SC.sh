#!/bin/bash
start=`date +%s.%N`
#All Boolean True/False parameters are 0 for False or 1 for True
#So python interprets them correctly

#Options for InitializeAltMTLs

#Random seed. Change to any integer you want (or leave the same)
seed=31415
#Number of realizations to generate. Ideally a multiple of 64 for bitweights
#However, you can choose smaller numbers for debugging
ndir=2
#Set to true(1) if you want to clobber already existing files for Alt MTL generation
overwrite=0
#Observing conditions to generate MTLs for (should be all caps "DARK" or "BRIGHT")
obscon='DARK'
#Survey to generate MTLs for (should be lowercase "sv3" or "main", sv2, sv1, and cmx are untested and will likely fail)
survey='main'
#Where to generate MTLs. Automatically formats number of MTLs into directory name but you can change this
date='061522'
printf -v outputMTLDirBase "$CSCRATCH/alt_mtls_archivedatesortTest%s_%03ddirs_%sRepro/" $date $ndir $survey
hpListFile='MainSurveyHPList.txt'
#These two options only are considered if the obscon is bright
#First option indicates whether to shuffle the top level priorities
#of BGS_FAINT/BGS_FAINT_HIP. Second option indicates what fraction/percent
#of BGS_FAINT to promote. Default is 20%
shuffleBrightPriorities=0
PromoteFracBGSFaint=0.2
#location of original MTLs to shuffle
exampleledgerbase=/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/
#Options for DateLoopAltMTL and runAltMTLParallel

#Quick Restart (i.e. reset the MTLs by copying the saved original shuffled files). 
#Default = 0/False. Set equal to 1 if you want to restart from the first observations
qR=0
#Number of observation dates to loop through
#Defaults to 40 dates for SV3
NObsDates=3
#Number of nodes to run on. This will launch up to 64*N jobs 
#if that number of alternate universes have already been generated
#Defaults to 1 for 64 directories
NNodes=1

#getosubp: grab subpriorities from the original (exampleledgerbase) MTLs
#This should only be turned on for testing/debugging purposes

getosubp=0

#Include secondary targets?
secondary=0

numobs_from_ledger=1
#Force redo fiber assignment if it has already been done. 
redoFA=0


#Options for MakeBitweightsParallel
#True/False(1/0) as to whether to split bitweight calculation
#among nodes by MPI between realizations
splitByReal=0
#Split the calculation of bitweights into splitByChunk
#chunks of healpixels. 
splitByChunk=100

#Set to true if you want to clobber already existing bitweight files
overwrite2=1
#Actual running of scripts

#Copy this script to output directory for reproducbility
thisFileName=$outputMTLDirBase/$0

echo $thisFileName

if [ -f "$thisFileName" ]
then
    echo "File is found. Checking to see it is identical to the original."
    cmp  $0 $thisFileName
    comp=$?
    if  [[ $comp -eq 1 ]]
    then 
        echo "Files are not identical."
        echo "If this is intended, please delete the original copied script at $thisFileName"
        echo "If this is unintended, you can reuse the original copied script at that same location"
        echo "goodbye"
        exit 3141
    elif [[ $comp -eq 0 ]] 
    then
        echo "files are same, continuing"
    else 
        echo "Something has gone very wrong. Exit code for cmp was $a"
        exit $a
    fi
else
   echo "Copied script is not found. Copying now, making directories as needed."
   mkdir -p $outputMTLDirBase
   cp $0 $outputMTLDirBase
fi

echo 'moving on to python scripts (REMOVE BEFORE PUSHING)'
printf -v OFIM "%s/InitializeAltMTLsParallelOutput_%sRepro%s.out" $outputMTLDirBase $survey $date

srun --cpu-bind=none --nodes=$NNodes -C haswell -A desi --qos=interactive -t 04:00:00 --mem=120000 InitializeAltMTLsParallel_DEBUG.py $seed $ndir $overwrite $obscon $survey $outputMTLDirBase $hpListFile $shuffleBrightPriorities $PromoteFracBGSFaint $exampleledgerbase $NNodes >& $OFIM
endInit=`date +%s.%N`

exit 123141
if [ $? -ne 0 ]; then
    exit 1234
fi
printf -v OFDL "%s/dateLoopAltMTLOutput_%sRepro%s.out" $outputMTLDirBase $survey $date
runtimeInit=$( echo "$endInit - $start" | bc -l )
bash dateLoopAltMTL.sh $qR $NObsDates $NNodes $outputMTLDirBase $secondary $obscon $survey $numobs_from_ledger $redoFA $getosubp  >& $OFDL
endDL=`date +%s.%N`

exit 31415926

if [ $? -ne 0 ]; then
    exit 12345
fi
if [ $splitByReal -ne 0 ]; then
    printf -v OFBW "%s/MakeBitweightsOutputCase1%sRepro%s.out" $outputMTLDirBase $survey $date
    srun --cpu-bind=none --nodes=$NNodes -C haswell -A desi --qos=interactive -t 04:00:00 --mem=120000 MakeBitweights.py $survey $obscon $ndir $splitByReal $splitByChunk $hpListFile $outputMTLDirBase $overwrite2 >& $OFBW
else
    printf -v OFBW "%s/MakeBitweightsOutputCase2%sRepro%s.out" $outputMTLDirBase $survey $date
    srun --cpu-bind=none --nodes=1 -C haswell -A desi --qos=interactive -t 04:00:00 --mem=120000 MakeBitweights.py $survey $obscon $ndir $splitByReal $splitByChunk $hpListFile $outputMTLDirBase $overwrite2 >& $OFBW
fi

endBW=`date +%s.%N`



runtimeInit=$( echo "$endInit - $start" | bc -l )
runtimeDateLoop=$( echo "$endDL - $endInit" | bc -l )
runtimeBitweights=$( echo "$endBW - $endDL" | bc -l )

echo "runtime for initialization"
echo $runtimeInit
echo "runtime for Dateloop of $NObsDates days"
echo $runtimeDateLoop
echo "runtime for making bitweights"
echo $runtimeBitweights