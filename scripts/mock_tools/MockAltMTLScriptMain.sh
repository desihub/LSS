#!/bin/bash
start=`date +%s.%N`

#ALTMTLHOME is a home directory for all of your alternate MTLs. Default is your scratch directory
#There will be an environment variable $ALTMTLHOME for the "survey alt MTLs"
#However, you should specify your own directory to a. not overwrite the survey alt MTLs 
# and b. keep your alt MTLs somewhere that you have control/access

#Uncomment the following line to set your own/nonscratch directory
ALTMTLHOME=/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit

#Mock realization
mockNumber=0

#simName is the subdirectory within ALTMTLHOME where this specific set of alt MTLs will be written
simName="altmtl$mockNumber"

#Location where you have cloned the LSS Repo
path2LSS=/pscratch/sd/a/acarnero/codes/LSS/bin/

#Number of realizations to generate. Ideally a multiple of 64 for bitweights
#However, you can choose smaller numbers for debugging
ndir=1

#Observing conditions for generating MTLs (should be all caps "DARK" or "BRIGHT")
obscon='DARK'

#Survey to generate MTLs for (should be lowercase "sv3" or "main", sv2, sv1, and cmx are untested and will likely fail)
survey='main'
#survey='sv3'
#StartDate options are default None (empty strings). Uncommenting the second options will set them to the SV3 start and end dates. 
startDate=''
endDate='--endDate=20220613' #'' june 13 2022 20220613

#2022-09-13 this is in data
#startDate=20210406
#endDate=20210625

# location of original MTLs to shuffle.
# Default directory is a read only mount of the CFS filesystem
# You can only access that directory from compute nodes. 
# Do NOT use the commented out directory (the normal mount of CFS)
# unless the read only mount is broken
#exampleLedgerBase=/dvs_ro/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/
#exampleLedgerBase=/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/
#exampleLedgerBase=/pscratch/sd/j/jlasker/MockAMTLY1/FirstGenMocks/AbacusSummit/mtls/
exampleLedgerBase=$ALTMTLHOME/$simName/initled
#Options for DateLoopAltMTL and runAltMTLParallel

#List of healpixels to create Alt MTLs for
hpListFile="$exampleLedgerBase/hpxlist_dark.txt"
#hpListFile="$path2LSS/MainSurveyHPList.txt"
#hpListFile="$path2LSS/SV3HPList.txt"


#If running from mocks, must set target directory. 
#Otherwise this is optional
#targfile='--targfile=/global/cfs/cdirs/desi/target/catalogs/dr9/1.1.1/targets/main/resolve/' #Main survey target directory
#targfile='--targfile=/cscratch/sd/j/jlasker/MockAMTLY1/FirstGenMocks/AbacusSummit/forFA1.fits' #WITHOUT PHOTSYS
#targfile='--targfile=/pscratch/sd/j/jlasker/MockAMTLY1/FirstGenMocks/AbacusSummit/TargetsWithNumobs_012322.fits' #WITHOUT PHOTSYS
targfile="--targfile=$ALTMTLHOME/forFA$mockNumber.fits"

# Flags for debug/verbose mode/profiling code time usage. 
# Uncomment second set of options to turn on the modes
#debug=''
#verbose=''
#profile=''
debug='--debug'
verbose='--verbose'
profile='--profile'

#Uncomment second option if running on mocks
#mock=''
mock='--mock'


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




#Options for InitializeAltMTLs

#Random seed. Change to any integer you want (or leave the same)
#If seed is different between two otherwise identical runs, the initial MTLs will also be different
#seed is also saved in output directory
seed=151

#Uncomment second option if you want to clobber already existing files for Alt MTL generation
overwrite=''
#overwrite='--overwrite'

#For rundate formatting in simName, either manually modify the string below 
#to be the desired date or comment that line out and uncomment the 
#following line to autogenerate date strings.
#To NOT use any date string specification, use the third line, an empty string
#datestring='071322'
#datestring=`date +%y%m%d`
datestring=''

#Can save time in MTL generation by first writing files to local tmp directory and then copying over later
#uncommenting the second option will directly write to your output directory
usetmp=''
#usetmp='--dontUseTemp'

if [ -z $usetmp ]
then
    outputMTLDirBaseBase=`mktemp -d /dev/shm/"$USER"_tempdirXXXX`
else 
    outputMTLDirBaseBase=$ALTMTLHOME
fi
printf -v outputMTLDirBase "$outputMTLDirBaseBase/$simName/" $datestring $ndir $survey
printf -v outputMTLFinalDestination "$ALTMTLHOME/$simName/" $datestring $ndir $survey

#These two options only are considered if the obscon is BRIGHT
#First option indicates whether to shuffle the top level priorities
#of BGS_FAINT/BGS_FAINT_HIP. Uncomment section option to turn off shuffling of bright time priorities
#Second option indicates what fraction/percent
#of BGS_FAINT to promote to BGS_FAINT_HIP. Default is 20%, same as SV3

shuffleBrightPriorities='--shuffleBrightPriorities'
#shuffleBrightPriorities=''

PromoteFracBGSFaint=0.2

#Quick Restart (i.e. reset the MTLs by copying the saved original shuffled files). 
#Default = Empty String/False. Uncomment second option if you want to restart from the first observations
echo "Fix QR resetting for new argparse usage"
qR=''
#qR='-qr'

#Number of observation dates to loop through
#Defaults to 40 dates for SV3
NObsDates=400

#Number of nodes to run on. This will launch up to 64*N jobs 
#if that number of alternate universes have already been generated
#Calculated automatically from number of sims requested and number of processes per node. Be careful if setting manually
NNodes=$(( ($ndir + $ProcPerNode - 1 )/$ProcPerNode ))
#echo $NNodes
#getosubp: grab subpriorities from the original (exampleledgerbase) MTLs
#This should only be turned on for SV testing/debugging purposes
#This should not be required for main survey debugging. 
getosubp=""
#getosubp='--getosubp'

#shuffleSubpriorities(reproducing) must be left as empty strings to ensure 
#subpriorities are shuffled. debug mode for main survey
#will only require these flags to be set by uncommenting second options

dontShuffleSubpriorities=''
reproducing=''
#dontShuffleSubpriorities='--dontShuffleSubpriorities'
#reproducing='--reproducing'
#Include secondary targets?
secondary=''
#secondary='--secondary'

#Default is use numobs from ledger. Uncomment second option to set numobs NOT from ledger
numobs_from_ledger=''
#numobs_from_ledger='--NumObsNotFromLedger'

#Uncomment second line to force redo fiber assignment if it has already been done. 
redoFA=''
#redoFA='--redoFA'


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
thisFileName=$outputMTLFinalDestination/$0

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
   mkdir -p $outputMTLFinalDestination
   cp $SLURM_SUBMIT_DIR $0 $outputMTLFinalDestination/$0
fi

if [ -d "$outputMTLFinalDestination" ]
then
    echo "output final directory exists"
    echo $outputMTLFinalDestination
else
   echo "output final directory does not exist. Creating and copying script there"
   mkdir -p $outputMTLFinalDestination
   cp $0 $outputMTLFinalDestination
fi
echo 'string county thingies'
echo ${getosubp}

if [ -z ${getosubp} ]
then
    echo 'getosubp not set'
else
    echo 'getosubp set'
    touch $outputMTLFinalDestination/GetOSubpTrue
fi

printf -v OFIM "%s/Initialize%sAltMTLsParallelOutput_%sRepro%s.out" $outputMTLFinalDestination $obscon $survey $date

echo "srun --nodes=$NNodes -C $CVal -q $QVal -A desi -t 02:00:00 --mem=120000 $path2LSS/InitializeAltMTLsParallel.py --seed=$seed --ndir=$ndir  --obscon=$obscon --survey=$survey --outputMTLDirBase=$outputMTLDirBase --promoteFracBGSFaint=$PromoteFracBGSFaint --HPListFile=$hpListFile --exampleLedgerBase=$exampleLedgerBase --ProcPerNode=$ProcPerNode     --finalDir="$outputMTLFinalDestination/Univ{0:03d}" $overwrite $shuffleBrightPriorities $usetmp $dontShuffleSubpriorities $reproducing $debug $verbose  $startDate $endDate >& $OFIM"
#srun --nodes=$NNodes -C $CVal -q $QVal -A desi -t 01:00:00 --mem=120000 $path2LSS/InitializeAltMTLsParallel.py --seed=$seed --ndir=$ndir  --obscon=$obscon --survey=$survey --outputMTLDirBase=$outputMTLDirBase --PromoteFracBGSFaint=$PromoteFracBGSFaint --HPListFile=$hpListFile --exampleLedgerBase=$exampleLedgerBase --ProcPerNode=$ProcPerNode --finalDir="$outputMTLFinalDestination/Univ{0:03d}" $overwrite $shuffleBrightPriorities $usetmp $dontShuffleSubpriorities $reproducing $debug $verbose  $startDate $endDate >& $OFIM
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
argstring="--altMTLBaseDir=$outputMTLFinalDestination --obscon=$obscon --survey=$survey --ProcPerNode=$ProcPerNode $numobs_from_ledger $redoFA $getosubp $debug $verbose $secondary $mock $targfile"
nohup bash $path2LSS/dateLoopAltMTL.sh $NObsDates $NNodes $path2LSS $CVal $QVal $qR $argstring  >& $OFDL

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

if [ "$ndir" -gt 1 ]; then 
	if [ $splitByReal -ne 0 ]; then
    		printf -v OFBW "%s/MakeBitweights%sOutputCase1%sRepro%s.out" $outputMTLFinalDestination $obscon $survey $datestring
    		srun --nodes=1 -C $CVal -q $QVal -A desi -t 04:00:00 --mem=120000 $path2LSS/MakeBitweights.py $survey $obscon $ndir $splitByReal $splitByChunk $hpListFile $outputMTLFinalDestination $overwrite2 >& $OFBW
	else
    		printf -v OFBW "%s/MakeBitweights%sOutputCase2%sRepro%s.out" $outputMTLFinalDestination $obscon $survey $datestring
    		srun --nodes=1 -C $CVal -q $QVal -A desi -t 04:00:00 --mem=120000 $path2LSS/MakeBitweights.py $survey $obscon $ndir $splitByReal $splitByChunk $hpListFile $outputMTLFinalDestination $overwrite2 >& $OFBW
	fi
	endBW=`date +%s.%N`
fi




runtimeInit=$( echo "$endInit - $start" | bc -l )
runtimeDateLoop=$( echo "$endDL - $endInit" | bc -l )

if [ "$ndir" -gt 1 ]; then
	runtimeBitweights=$( echo "$endBW - $endDL" | bc -l )
fi

echo "runtime for initialization"
echo $runtimeInit
echo "runtime for Dateloop of $NObsDates days"
echo $runtimeDateLoop
if [ "$ndir" -gt 1 ]; then  
	echo "runtime for making bitweights"
	echo $runtimeBitweights
fi


