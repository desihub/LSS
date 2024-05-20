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

    srun --nodes=$NNodes -C $CVal --qos=$QVal -A desi -t 02:00:00 -e mylog_error.txt --dependency=afterany:17881308 python $path2LSS/runAltMTLParallel.py $argstring
fi
if [ $QVal = 'regular' ]; 
then

    srun --nodes=$NNodes -C $CVal --qos=$QVal -A desi -t 24:00:00 -e mylog_error.txt -K --dependency=afterany:24949299 python  $path2LSS/runAltMTLParallel.py $argstring
fi

if [ $QVal = 'debug' ]; 
then

    srun --nodes=$NNodes -C $CVal --qos=$QVal -A desi -t 00:15:00 -e mylog_error.txt --dependency=afterany:17881308 python  $path2LSS/runAltMTLParallel.py $argstring
    #srun --nodes=$NNodes -C $CVal --qos=$QVal -A desi -t 00:15:00 --mem-per-cpu=1500 -K -e mylog_error.txt --dependency=afterany:17881308  $path2LSS/runAltMTLParallel.py $argstring
fi

#if [ $QVal = 'shared' ];
#then

#    srun --ntasks=128 -C $CVal --qos=$QVal -A desi -t 24:00:00 --mem=64G $path2LSS/runAltMTLParallel.py $argstring
#fi

#retcode=$?
#qR=0 #DO NOT CHANGE. This prevents further restarts after the first if qR is set to 1 at top.
#if [ $retcode -ne 0 ]; then
#    echo 'something went wrong'
#    echo $retcode
#    exit 1234
#fi

#done
