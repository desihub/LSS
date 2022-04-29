#!/bin/bash
#All Boolean True/False parameters are 0 for False or 1 for True
#So python interprets them correctly

#Options for InitializeAltMTLs

#Random seed. Change to any integer you want (or leave the same)
seed=12345
#Number of realizations to generate. Ideally a multiple of 64 for bitweights
#However, you can choose smaller numbers for debugging
ndir=256
#Set to true(1) if you want to clobber already existing files for Alt MTL generation
overwrite=0
#Observing conditions to generate MTLs for (should be all caps "DARK" or "BRIGHT")
obscon='DARK'
#Survey to generate MTLs for (should be lowercase "sv3" or "main", sv2, sv1, and cmx are untested and will likely fail)
survey='sv3'
#mockrea=1 THIS IS THE MOCK REALIZATION

#mockrea=0 #THIS IS FOR EACH MOCK REALIZATION

for mockrea in {16..24} 
###for mockrea in {0..24} 
do
#Where to generate MTLs. Automatically formats number of MTLs into directory name but you can change this
printf -v outputMTLDirBase "$CSCRATCH/alt_mtls_masterScriptTest_%03ddirs_rea%03d/" $ndir $mockrea
hpListFile='SV3HPList.txt'
#hpListFile='SV3HPList_mock.txt'
#These two options only are considered if the obscon is bright
#First option indicates whether to shuffle the top level priorities
#of BGS_FAINT/BGS_FAINT_HIP. Second option indicates what fraction/percent
#of BGS_FAINT to promote. Default is 20%
shuffleBrightPriorities=0
PromoteFracBGSFaint=0.2
#location of original MTLs to shuffle
printf -v exampleledgerbase "$CSCRATCH/mtl_test/init_mock%03d/" $mockrea
#Options for DateLoopAltMTL and runAltMTLParallel

#Quick Restart (i.e. reset the MTLs by copying the saved original shuffled files). 
#Default = 0/False. Set equal to 1 if you want to restart from the first observations
qR=0
#Number of observation dates to loop through
#Defaults to 40 dates for SV3
NObsDates=40
#Number of nodes to run on. This will launch up to 64*N jobs 
#if that number of alternate universes have already been generated
#Defaults to 1 for 64 directories
NNodes=4


#Include secondary targets?
secondary=0

numobs_from_ledger=1
#Force redo fiber assignment if it has already been done. 
redoFA=1


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



bash dateLoopAltMTL_mock.sh $qR $NObsDates $NNodes $outputMTLDirBase $secondary $obscon $survey $numobs_from_ledger $redoFA >& dateLoopAltMTLOutput.out
if [ $? -ne 0 ]; then
    exit 12345
fi


if [ $splitByReal -ne 0 ]; then
    srun --nodes=$NNodes -C haswell -A desi --qos=interactive -t 04:00:00 --mem=120000 MakeBitweights_mock.py $survey $obscon $ndir $splitByReal $splitByChunk $hpListFile $outputMTLDirBase $overwrite2 $exampleledgerbase >& MakeBitweightsOutput.out
else
    srun --nodes=1 -C haswell -A desi --qos=interactive -t 04:00:00 --mem=120000 MakeBitweights_mock.py $survey $obscon $ndir $splitByReal $splitByChunk $hpListFile $outputMTLDirBase $overwrite2 $exampleledgerbase >& MakeBitweightsOutput.out
fi

done
