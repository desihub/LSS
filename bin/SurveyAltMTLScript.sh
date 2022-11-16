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

#Uncomment the following line to set your own/nonscratch directory
#ALTMTLHOME=/path/to/your/directory/


if [[ "${NERSC_HOST}" == "cori" ]]; then
    CVal='haswell'
    QVal='interactive'
    ProcPerNode=32
    if [[ -z "${ALTMTLHOME}" ]]; then
        ALTMTLHOME=$CSCRATCH
    else
        echo "ALTMTLHOME Already set. ALTMTLHOME=$ALTMTLHOME"
    fi
elif [[ "${NERSC_HOST}" == "perlmutter" ]]; then
    srunConfig='-C cpu -q regular'
    CVal='cpu'
    QVal='interactive'
    ProcPerNode=128
    if [[ -z "${ALTMTLHOME}" ]]; then
        ALTMTLHOME=$PSCRATCH
    else
        echo "ALTMTLHOME Already set. ALTMTLHOME=$ALTMTLHOME"
    fi

else
    echo "This code is only supported on NERSC Cori and NERSC Perlmutter. Goodbye"
    exit 1234
fi


#simName is the subdirectory within ALTMTLHOME where this specific set of alt MTLs will be written
simName="$USER"_SV3TestStartEndDateAllDates

#Options for InitializeAltMTLs

#Random seed. Change to any integer you want (or leave the same)
#If seed is different between two otherwise identical runs, the initial MTLs will also be different
#seed is also saved in output directory
seed=314159

#Number of realizations to generate. Ideally a multiple of 64 for bitweights
#However, you can choose smaller numbers for debugging
ndir=2

#Set to true(1) if you want to clobber already existing files for Alt MTL generation
overwrite=0

#Observing conditions for generating MTLs (should be all caps "DARK" or "BRIGHT")
obscon='DARK'

#Survey to generate MTLs for (should be lowercase "sv3" or "main", sv2, sv1, and cmx are untested and will likely fail)
survey='sv3'
startDate=20210406
endDate=20210625

#For rundate formatting in simName, either manually modify the string below 
#to be the desired date or comment that line out and uncomment the 
#following line to autogenerate date strings.
#To NOT use any date string specification, use the third line, an empty string
#datestring='071322'
#datestring=`date +%y%m%d`
datestring=''

#Can save time in MTL generation by first writing files to local tmp directory and then copying over later
#usetmp=True will use the local tmp directory and usetmp=False will directly write to your output directory
usetmp=True

if [ $usetmp ]
then
    outputMTLDirBaseBase=`mktemp -d /dev/shm/"$USER"_tempdirXXXX`
else 
    outputMTLDirBaseBase=$ALTMTLHOME
fi
printf -v outputMTLDirBase "$outputMTLDirBaseBase/$simName/" $datestring $ndir $survey
printf -v outputMTLFinalDestination "$ALTMTLHOME/$simName/" $datestring $ndir $survey

#List of healpixels to create Alt MTLs for
#hpListFile="$path2LSS/MainSurveyHPList.txt"
hpListFile="$path2LSS/SV3HPList.txt"

#These two options only are considered if the obscon is BRIGHT
#First option indicates whether to shuffle the top level priorities
#of BGS_FAINT/BGS_FAINT_HIP. Second option indicates what fraction/percent
#of BGS_FAINT to promote to BGS_FAINT_HIP. Default is 20%, same as SV3
shuffleBrightPriorities=0
PromoteFracBGSFaint=0.2

# location of original MTLs to shuffle.
# Default directory is a read only mount of the CFS filesystem
# You can only access that directory from compute nodes. 
# Do NOT use the commented out directory (the normal mount of CFS)
# unless the read only mount is broken
#exampleledgerbase=/dvs_ro/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/
exampleledgerbase=/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/

#Options for DateLoopAltMTL and runAltMTLParallel

#Quick Restart (i.e. reset the MTLs by copying the saved original shuffled files). 
#Default = 0/False. Set equal to 1 if you want to restart from the first observations
qR=0

#Number of observation dates to loop through
#Defaults to 40 dates for SV3
NObsDates=40
#Number of nodes to run on. This will launch up to 64*N jobs 
#if that number of alternate universes have already been generated
#Defaults to 4 for 128 directories
NNodes=1

#getosubp: grab subpriorities from the original (exampleledgerbase) MTLs
#This should only be turned on for SV testing/debugging purposes
#This should not be required for main survey debugging. 
getosubp=0

#shuffleSubpriorities(reproducing) must be set to 1(0) to ensure 
#subpriorities are shuffled. debug mode for main survey
#will only require these flags to be set to 0(1) and not the getosubp flag
shuffleSubpriorities=1
reproducing=0

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

if [ $getosubp -gt 0 ]
then
    touch $outputMTLFinalDestination/GetOSubpTrue
fi

echo 'moving on to python scripts (REMOVE BEFORE PUSHING)'
printf -v OFIM "%s/Initialize%sAltMTLsParallelOutput_%sRepro%s.out" $outputMTLFinalDestination $obscon $survey $date

srun --nodes=$NNodes -C $CVal -q $QVal -A desi -t 04:00:00 --mem=120000 $path2LSS/InitializeAltMTLsParallel.py $seed $ndir $overwrite $obscon $survey $outputMTLDirBase $hpListFile $shuffleBrightPriorities $PromoteFracBGSFaint $exampleledgerbase $NNodes $usetmp "$outputMTLFinalDestination/Univ{0:03d}" $shuffleSubpriorities $reproducing $debug $verbose $ProcPerNode $startDate $endDate >& $OFIM
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

printf -v OFDL "%s/dateLoop%sAltMTLOutput_%sRepro%s.out" $outputMTLFinalDestination $obscon $survey $datestring

runtimeInit=$( echo "$endInit - $start" | bc -l )

nohup bash $path2LSS/dateLoopAltMTL.sh $qR $NObsDates $NNodes $outputMTLFinalDestination $secondary $obscon $survey $numobs_from_ledger $redoFA $getosubp $path2LSS $CVal $QVal $debug $verbose $ProcPerNode >& $OFDL

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
    printf -v OFBW "%s/MakeBitweights%sOutputCase1%sRepro%s.out" $outputMTLFinalDestination $obscon $survey $datestring
    srun --nodes=1 -C $CVal -q $QVal -A desi -t 04:00:00 --mem=120000 $path2LSS/MakeBitweights.py $survey $obscon $ndir $splitByReal $splitByChunk $hpListFile $outputMTLFinalDestination $overwrite2 >& $OFBW
else
    printf -v OFBW "%s/MakeBitweights%sOutputCase2%sRepro%s.out" $outputMTLFinalDestination $obscon $survey $datestring
    srun --nodes=1 -C $CVal -q $QVal -A desi -t 04:00:00 --mem=120000 $path2LSS/MakeBitweights.py $survey $obscon $ndir $splitByReal $splitByChunk $hpListFile $outputMTLFinalDestination $overwrite2 >& $OFBW
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
