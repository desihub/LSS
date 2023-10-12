#!/bin/bash

echo "All Arguments"
echo $@

NObsDates=$1

NNodes=$2

path2LSS=$3

CVal=$4

QVal=$5

argstring=${@:6}

echo 'argstring'
echo "$argstring"


#for i in $(seq 0 1 $NObsDates)
#do
#    echo " NextDate"
#    echo ""
#    echo ""
#    echo ""
#    echo $i
#    echo ""
#    echo ""
#    echo ""
if [ $QVal = 'interactive' ]; 
then

    srun --nodes=$NNodes -C $CVal --qos=$QVal -A desi -t 03:00:00 --dependency=afterany:16533190 $path2LSS/runAltMTLParallel.py $argstring
fi
if [ $QVal = 'regular' ]; 
then

    srun --nodes=$NNodes -C $CVal --qos=$QVal -A desi -t 12:00:00 --dependency=afterany:16533190 $path2LSS/runAltMTLParallel.py $argstring
fi

if [ $QVal = 'debug' ]; 
then

    srun --nodes=$NNodes -C $CVal --qos=$QVal -A desi -t 00:15:00 --dependency=afterany:16533190 $path2LSS/runAltMTLParallel.py $argstring
fi
#retcode=$?
#qR=0 #DO NOT CHANGE. This prevents further restarts after the first if qR is set to 1 at top.
#if [ $retcode -ne 0 ]; then
#    echo 'something went wrong'
#    echo $retcode
#    exit 1234
#fi

#done