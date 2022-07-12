#!/bin/bash
start=`date +%s.%N`
#All Boolean True/False parameters are 0 for False or 1 for True
#So python interprets them correctly

#Location where you have cloned the LSS Repo
path2LSS=~/.local/desicode/LSS/bin/

#Flags for debug/verbose mode/profiling code time usage
debug=0
verbose=0
profile=0

#ALTMTLHOME is a home directory for all of your alternate MTLs. Default is your scratch directory
#There will be an environment variable $ALTMTLHOME for the "survey alt MTLs"
#However, you should specify your own directory to a. not overwrite the survey alt MTLs 
# and b. keep your alt MTLs somewhere that you have control/access
if [[ "${NERSC_HOST}" == "cori" ]]; then
    ALTMTLHOME=$CSCRATCH
elif [[ "${NERSC_HOST}" == "perlmutter" ]]; then
    ALTMTLHOME=$PSCRATCH
else
    echo "Something went wrong. Goodbye"
    exit 1234
fi

#simName is the subdirectory within ALTMTLHOME where this specific set of alt MTLs will be written
simName="alt_mtls_OrigTimingTest%s_%03ddirs_%sRepro" 

#Options for InitializeAltMTLs

#Random seed. Change to any integer you want (or leave the same)
#If seed is different between two otherwise identical runs, the initial MTLs will also be different
#seed is also saved in output directory
seed=31415

#Number of realizations to generate. Ideally a multiple of 64 for bitweights
#However, you can choose smaller numbers for debugging
ndir=2

#Set to true(1) if you want to clobber already existing files for Alt MTL generation
overwrite=0

#Observing conditions for generating MTLs (should be all caps "DARK" or "BRIGHT")
obscon='DARK'

#Survey to generate MTLs for (should be lowercase "sv3" or "main", sv2, sv1, and cmx are untested and will likely fail)
survey='sv3'

#For rundate formatting in simName, either manually modify the string below 
#to be the desired date or comment that line out and uncomment the 
#following line to autogenerate date strings.
#To NOT use any date string specification, use the third line,  an empty string
datestring='071222'
#datestring=`date +%y%m%d`
#datestring=''

#Can save time in MTL generation by first writing files to local tmp directory and then copying over later
#usetmp=True will use the local tmp directory and usetmp=False will directly write to your output directory
usetmp=True

if [ usetmp ]
then
    outputMTLDirBaseBase=`mktemp -d /dev/shm/JLtempdirXXXX`
else 
    outputMTLDirBaseBase=$ALTMTLHOME
fi
printf -v outputMTLDirBase "$outputMTLDirBaseBase/$simName/" $datestring $ndir $survey
printf -v outputMTLFinalDestination "$ALTMTLHOME/$simName/" $datestring $ndir $survey

hpListFile='SV3HPList.txt'

#These two options only are considered if the obscon is bright
#First option indicates whether to shuffle the top level priorities
#of BGS_FAINT/BGS_FAINT_HIP. Second option indicates what fraction/percent
#of BGS_FAINT to promote. Default is 20%
shuffleBrightPriorities=0
PromoteFracBGSFaint=0.2

# location of original MTLs to shuffle.
# Default directory is a read only mount of the CFS filesystem
# You can only access that directory from compute nodes. 
# Do NOT use the commented out directory (the normal mount of CFS)
# unless the read only mount is broken
exampleledgerbase=/dvs_ro/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/
#exampleledgerbase=/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/

#Options for DateLoopAltMTL and runAltMTLParallel

#Quick Restart (i.e. reset the MTLs by copying the saved original shuffled files). 
#Default = 0/False. Set equal to 1 if you want to restart from the first observations
qR=0
#Number of observation dates to loop through
#Defaults to 33 dates for SV3
NObsDates=33
#Number of nodes to run on. This will launch up to 64*N jobs 
#if that number of alternate universes have already been generated
#Defaults to 1 for 64 directories
NNodes=1

#getosubp: grab subpriorities from the original (exampleledgerbase) MTLs
#This should only be turned on for testing/debugging purposes
getosubp=1

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
        echo "If this is intended, please delete or edit the original copied script at $thisFileName"
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
fi

if [ -d "$outputMTLFinalDestination" ]
then
    echo "output final directory exists"
else
   echo "output final directory does not exist. Creating and copying script there"
   mkdir -p $outputMTLFinalDestination
   cp $0 $outputMTLFinalDestination
fi

echo 'moving on to python scripts (REMOVE BEFORE PUSHING)'
printf -v OFIM "%s/InitializeAltMTLsParallelOutput_%sRepro%s.out" $outputMTLFinalDestination $survey $date

srun --nodes=$NNodes -C haswell -A desi --qos=interactive -t 04:00:00 --mem=120000 $path2LSS/InitializeAltMTLsParallel.py $seed $ndir $overwrite $obscon $survey $outputMTLDirBase $hpListFile $shuffleBrightPriorities $PromoteFracBGSFaint $exampleledgerbase $NNodes $usetmp "$outputMTLFinalDestination/Univ{0:03d}" >& $OFIM

if [ $? -ne 0 ]; then
    exit 1234
    endInit=`date +%s.%N`
    runtimeInit=$( echo "$endInit - $start" | bc -l )
    echo "runtime for initialization"
    echo $runtimeInit
fi

endInit=`date +%s.%N`
runtimeInit=$( echo "$endInit - $start" | bc -l )
echo "runtime for initialization"
echo $runtimeInit

printf -v OFDL "%s/dateLoopAltMTLOutput_%sRepro%s.out" $outputMTLFinalDestination $survey $datestring
runtimeInit=$( echo "$endInit - $start" | bc -l )
nohup bash $path2LSS/dateLoopAltMTL.sh $qR $NObsDates $NNodes $outputMTLFinalDestination $secondary $obscon $survey $numobs_from_ledger $redoFA $getosubp  >& $OFDL
endDL=`date +%s.%N`

if [ $? -ne 0 ]; then
    runtimeDateLoop=$( echo "$endDL - $endInit" | bc -l )
    echo "runtime for Dateloop of $NObsDates days"
    echo $runtimeDateLoop
    exit 12345
fi

runtimeDateLoop=$( echo "$endDL - $endInit" | bc -l )
echo "runtime for Dateloop of $NObsDates days"
echo $runtimeDateLoop

if [ $splitByReal -ne 0 ]; then
    printf -v OFBW "%s/MakeBitweightsOutputCase1%sRepro%s.out" $outputMTLFinalDestination $survey $datestring
    srun --nodes=$NNodes -C haswell -A desi --qos=interactive -t 04:00:00 --mem=120000 MakeBitweights.py $survey $obscon $ndir $splitByReal $splitByChunk $hpListFile $outputMTLDirBase $overwrite2 >& $OFBW
else
    printf -v OFBW "%s/MakeBitweightsOutputCase2%sRepro%s.out" $outputMTLFinalDestination $survey $datestring
    srun --nodes=1 -C haswell -A desi --qos=interactive -t 04:00:00 --mem=120000 MakeBitweights.py $survey $obscon $ndir $splitByReal $splitByChunk $hpListFile $outputMTLDirBase $overwrite2 >& $OFBW
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