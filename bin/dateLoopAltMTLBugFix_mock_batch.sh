#!/bin/bash


export OMP_NUM_THREADS=1
### THIS IS WHAT COMES FROM dataLoopAltMTLBugFix_mock -------------------------------
#argstring="--altMTLBaseDir=$outputMTLFinalDestination --obscon=$obscon --survey=$survey --ProcPerNode=$ProcPerNode $numobs_from_ledger $redoFA $getosubp $debug $verbose $secondary $mock $targfile $multiDate $reproducing --mockmin=$mockinit --mockmax=$mockend --initpath=$initpath"
#echo 'argstring for dateloop'
#echo $argstring
#nohup bash $path2LSS/dateLoopAltMTLBugFix_mock.sh $NObsDates $NNodes $path2LSS $CVal $QVal $qR $argstring  >& $OFDL
#--------------------

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

    srun --nodes=1 -n 1 -C $CVal --qos=interactive -A desi -t 04:00:00 --dependency=afterany:46220452  python $path2LSS/runAltMTLRealizations.py $argstring
    #srun --nodes=$NNodes -C $CVal --qos=$QVal -A desi -t 03:00:00 --dependency=afterany:17881308 $path2LSS/runAltMTLParallel.py $argstring
fi
if [ $QVal = 'regular' ]; 
then
    #echo "srun -c 1 -n 1 -C $CVal --qos=shared -A desi -t 48:00:00 --mem=128G python $path2LSS/runAltMTLRealizations.py $argstring"
	echo "srun --nodes=1 -n 1 -C $CVal --qos=regular -A desi -t 48:00:00 python $path2LSS/runAltMTLRealizations.py $argstring"
	srun --nodes=1 -n 1 -C $CVal --qos=regular -A desi -t 48:00:00 python $path2LSS/runAltMTLRealizations.py $argstring
	#python $path2LSS/runAltMTLRealizations.py $argstring
    ##srun --nodes=$NNodes -C $CVal --qos=$QVal -A desi -t 12:00:00 --dependency=afterany:22562205 $path2LSS/runAltMTLRealizations.py $argstring
    #srun --nodes=$NNodes -C $CVal --qos=$QVal -A desi -t 12:00:00 --dependency=afterany:17881308 $path2LSS/runAltMTLParallel.py $argstring
fi

if [ $QVal = 'debug' ]; 
then

    srun --nodes=$NNodes -C $CVal --qos=$QVal -A desi -t 00:15:00 python $path2LSS/runAltMTLRealizations.py $argstring
    #srun --nodes=$NNodes -C $CVal --qos=$QVal -A desi -t 00:15:00 --dependency=afterany:17881308 $path2LSS/runAltMTLParallel.py $argstring
fi
#retcode=$?
#qR=0 #DO NOT CHANGE. This prevents further restarts after the first if qR is set to 1 at top.
#if [ $retcode -ne 0 ]; then
#    echo 'something went wrong'
#    echo $retcode
#    exit 1234
#fi

#done
