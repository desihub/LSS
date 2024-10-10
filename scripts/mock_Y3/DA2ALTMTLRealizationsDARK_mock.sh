#!/bin/bash
start=`date +%s.%N`

##TEMPrealization=0

#simName is the subdirectory within ALTMTLHOME where this specific set of alt MTLs will be written
#simName=JL_DebugReprocReprod2
simName="altmtl{mock_number}"
#Location where you have cloned the LSS Repo
path2LSS=/pscratch/sd/a/acarnero/codes/LSS/bin/

# Flags for debug/verbose mode/profiling code time usage. 
# Uncomment second set of options to turn on the modes
debug=''
#verbose=''
profile=''
#debug='--debug'
verbose='--verbose'
#profile='--profile'


#Uncomment second option if running on mocks
#mock=''
mock='--mock'


#Uncomment the following line to set your own/nonscratch directory
#ALTMTLHOME=/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/altmtl/
#ALTMTLHOME=/pscratch/sd/a/acarnero/test_main/
ALTMTLHOME=/global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummit_v4_1/

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
seed=14126579
#seed=3593589
#Number of realizations to generate. Ideally a multiple of 64 for bitweights
#However, you can choose smaller numbers for debugging
#Mock realization
mockinit=5
mockend=6
let ndir=$mockend-$mockinit 


#Uncomment second option if you want to clobber already existing files for Alt MTL generation
overwrite=''
#overwrite='--overwrite'

#Observing conditions for generating MTLs (should be all caps "DARK" or "BRIGHT")
obscon='DARK'
#obscon='BRIGHT'

#Survey to generate MTLs for (should be lowercase "sv3" or "main", sv2, sv1, and cmx are untested and will likely fail)
#survey='sv3'
survey='main'
# options are default None (empty strings). Uncommenting the second options will set them to the Y1 start and end dates. 
startDate=''
#endDate=''
#startDate='2021-05-13T08:15:37+00:00'
#endDate='2022-06-24T00:00:00+00:00'

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

#List of healpixels to create Alt MTLs for
#hpListFile="$path2LSS/MainSurveyHPList_mock.txt"
##TEMPhpListFile="/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_v3/altmtl$realization/initled/hpxlist_dark.txt"
#hpListFile="$path2LSS/MainSurveyHPList.txt"
#hpListFile="$path2LSS/DebugMainHPList.txt"
#hpListFile="$path2LSS/SV3HPList.txt"

#These two options only are considered if the obscon is BRIGHT
#First option indicates whether to shuffle the top level priorities
#of BGS_FAINT/BGS_FAINT_HIP. Uncomment section option to turn off shuffling of bright time priorities
#Second option indicates what fraction/percent
#of BGS_FAINT to promote to BGS_FAINT_HIP. Default is 20%, same as SV3

#shuffleBrightPriorities='--shuffleBrightPriorities'
shuffleBrightPriorities=''


shuffleELGPriorities=''
#shuffleELGPriorities='--shuffleELGPriorities'

#PromoteFracBGSFaint=0.2
PromoteFracBGSFaint=0.0
#PromoteFracELG=0.1
PromoteFracELG=0.

# location of original MTLs to shuffle.
# Default directory is a read only mount of the CFS filesystem
# You can only access that directory from compute nodes. 
# Do NOT use the commented out directory (the normal mount of CFS)
# unless the read only mount is broken
##TEMPexampleLedgerBase=/dvs_ro/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_v3/altmtl$realization/initled
#exampleLedgerBase=/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/
#exampleLedgerBase=/pscratch/sd/j/jlasker/MockAMTLY1/FirstGenMocks/AbacusSummit/mtls/
#exampleLedgerBase=$SCRATCH/MockAMTLY1/FirstGenMocks/AbacusSummit/mtls/
#Options for DateLoopAltMTL and runAltMTLParallel

#Quick Restart (i.e. reset the MTLs by copying the saved original shuffled files). 
#Default = Empty String/False. Uncomment second option if you want to restart from the first observations
#PLEASE DO NOT CHANGEME
echo "Fix QR resetting for new argparse usage"
qR=''
#qR='-qr'

#Number of observation dates to loop through
#Defaults to 40 dates for SV3
NObsDates=99999

# Whether to submit a new job with dateLoopAltMTL for each date
# or to submit a single job
# multiDate=0
multiDate='--multiDate'
echo 'setting QVal here for debug. Fix later.'
#QVal='debug'
QVal='regular'
#QVal='interactive'
#



#Number of nodes to run on. This will launch up to 64*N jobs 
#if that number of alternate universes have already been generated
#Calculated automatically from number of sims requested and number of processes per node. Be careful if setting manually
NNodes=$(( ($ndir + $ProcPerNode - 1 )/$ProcPerNode ))
#echo $NNodes
#getosubp: grab subpriorities from the original (exampleledgerbase) MTLs
#This should only be turned on for SV testing/debugging purposes
#This should not be required for main survey debugging. 
getosubp=''
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


#If running from mocks, must set target directory. 
#Otherwise this is optional
#targfile='' #CHANGEME IF RUNNING ON MOCKS
#targfile='--targfile=/global/cfs/cdirs/desi/target/catalogs/dr9/1.1.1/targets/main/resolve/' #Main survey target directory
#targfile="--targfile=/pscratch/sd/a/acarnero/test_main/forFA{mock_number}.fits"
targfile="--targfile=$DESI_ROOT/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummit_v4_1/forFA{mock_number}.fits"
#targfile='--targfile=CHANGEME IF RUNNING ON MOCKS' #/pscratch/sd/j/jlasker/MockAMTLY1/FirstGenMocks/AbacusSummit/forFA2.fits' 


#Default is use numobs from ledger. Uncomment second option to set numobs NOT from ledger
numobs_from_ledger=''
#numobs_from_ledger='--NumObsNotFromLedger'

#Uncomment second line to force redo fiber assignment if it has already been done. 
redoFA=''
#redoFA='--redoFA'


#Options for MakeBitweightsParallel
#True/False(1/0) as to whether to split bitweight calculation
#among nodes by MPI between realizations
#splitByReal=1

#Split the calculation of bitweights into splitByChunk
#chunks of healpixels. 
#splitByChunk=1

#Set to true (1) if you want to clobber already existing bitweight files
overwrite2=''
#overwrite2='--overwrite'
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

if [ -z $getosubp ]
then
    touch $outputMTLFinalDestination/GetOSubpTrue
fi

printf -v OFDL "%s/dateLoop%sAltMTLOutput_%sRepro%s.out" $outputMTLFinalDestination $obscon $survey $datestring

runtimeInit=$( echo "$endInit - $start" | bc -l )
argstring="--altMTLBaseDir=$outputMTLFinalDestination --obscon=$obscon --survey=$survey --ProcPerNode=$ProcPerNode $numobs_from_ledger $redoFA $getosubp $debug $verbose $secondary $mock $targfile $multiDate $reproducing --mockmin=$mockinit --mockmax=$mockend"
echo 'argstring for dateloop'
echo $argstring
nohup bash $path2LSS/dateLoopAltMTLBugFix_mock_batch.sh $NObsDates $NNodes $path2LSS $CVal $QVal $qR $argstring  >& $OFDL

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
exit 54321



runtimeDateLoop=$( echo "$endDL - $endInit" | bc -l )

echo "runtime for Dateloop of $NObsDates days"
echo $runtimeDateLoop
