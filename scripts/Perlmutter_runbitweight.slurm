#!/bin/bash
#SBATCH -N 1
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -J bitweights
#SBATCH --mail-user=your@email.com
#SBATCH --mail-type=ALL
#SBATCH -t 1:00:00

#OpenMP settings:
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

#run the application:
source /global/common/software/desi/desi_environment.sh 
PYTHONPATH=$PYTHONPATH:$HOME/LSS

srun -n 32 -c 8 --cpu_bind=cores $HOME/LSS/bin/mpi_bitweights --mtl inputfolder/targ.fits --tiles inputfolder/tiles.fits --format fits --outdir outputfolder --realizations 128
