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
    
    $path2LSS/runAltMTLParallaure.py $argstring
#    srun --immediate=14400 --nodes=$NNodes -C $CVal --qos=$QVal -A desi -t 02:00:00 $path2LSS/runAltMTLParallel.py $argstring
    retcode=$?
    qR=0 #DO NOT CHANGE. This prevents further restarts after the first if qR is set to 1 at top.
    if [ $retcode -ne 0 ]; then
        echo 'something went wrong'
        echo $retcode
        exit 1234
    fi

done
