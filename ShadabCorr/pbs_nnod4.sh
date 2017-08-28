#PBS -l nodes=1:ppn=16 #Dont change anything here independent of number of nodes
#PBS -q physics
#PBS -j eo
#PBS -N CFmultinode
#PBS -l walltime=20:00:00

echo 'This job started on: ' `date`
cd $PBS_O_WORKDIR

module load python27
module load python27-extras


nodeid=$((4-$PBS_ARRAYID)) #The four here should be replace by the number of nodes you want to use

#The nnode variable should be same as number of nodes you want to use
python Runme_Correlation.py -data datafile -rand randomfile -out outroot  -nproc 16 -nnode 4 -nodeid $nodeid 

echo 'This job ended on: ' `date`
