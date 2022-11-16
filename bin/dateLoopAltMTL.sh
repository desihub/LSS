#!/bin/bash

echo "All Arguments"
echo $@


#Quick Restart (i.e. reset the MTLs by copying the saved original shuffled files). 
#Default = 0/False. Set equal to 1 if you want to restart from the first observations
qR=$1

#Number of observation dates to loop through
NObsDates=$2
#Number of nodes to run on. This will launch up to 64*N jobs 
#if that number of alternate universes have already been generated
NNodes=$3

#Base directory for the alternate MTLs created in the InitializeAltMTLs script
altmtlbasedir=$4

#Include secondary targets?
secondary=$5

#Observing conditions to process the observations
obscon=$6

#Survey whose observations you are processing
survey=$7

numobs_from_ledger=$8

#Force redo fiber assignment if it has already been done. 
redoFA=$9

#getosubp: grab subpriorities from the original (exampleledgerbase) MTLs
#This should only be turned on for testing/debugging purposes

getosubp=${10}

path2LSS=${11}

CVal=${12}

QVal=${13}

debug=${14}

verbose=${15}

ProcPerNode=${16}

echo $NObsDates

for i in $(seq 0 1 $NObsDates)
do
    echo " NextDate"
    echo ""
    echo ""
    echo ""
    echo $i
    echo ""
    echo ""
    echo ""
    #srun  runAltMTLParallel.py $i
    echo "calling the function"
    echo "srun --nodes=$NNodes -C $CVal --qos=$QVal -A desi -t 02:00:00 $path2LSS/runAltMTLParallel.py $NNodes $qR $altmtlbasedir $secondary $obscon $survey $numobs_from_ledger $redoFA $getosubp $debug $verbose $ProcPerNode" 
    #srun --cpu-bind=none --nodes=$NNodes -C $CVal -q $QVal -A desi -t 02:00:00 $path2LSS/runAltMTLParallel.py $NNodes $qR $altmtlbasedir $secondary $obscon $survey $numobs_from_ledger $redoFA $getosubp $debug $verbose $ProcPerNode
    echo "CVal"
    echo $CVal
    echo "QVal"
    echo $QVal
    srun --nodes=$NNodes -C $CVal --qos=$QVal -A desi -t 02:00:00 $path2LSS/runAltMTLParallel.py $NNodes $qR $altmtlbasedir $secondary $obscon $survey $numobs_from_ledger $redoFA $getosubp $debug $verbose $ProcPerNode
    retcode=$?
    qR=0 #DO NOT CHANGE. This prevents further restarts after the first if qR is set to 1 at top.
    echo 'retcode'
    echo $retcode
    if [ $retcode -ne 0 ]; then
        echo 'something went wrong'
        echo $retcode
        exit 1234
    fi

done