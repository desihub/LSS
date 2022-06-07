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
    srun --nodes=$NNodes -C haswell -A desi --qos=interactive -t 02:00:00 runAltMTLParallel_mock.py $NNodes $qR $altmtlbasedir $secondary $obscon $survey $numobs_from_ledger $redoFA 
    qR=0 #DO NOT CHANGE. This prevents further restarts after the first if qR is set to 1 at top.
    if [ $? -ne 0 ]; then
        exit 1234
    fi
done
