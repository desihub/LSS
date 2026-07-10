#!/bin/bash
#SBATCH --ntasks=1 # force srun with option -n 1 , ie 1 task
#SBATCH --cpus-per-task=3
#SBATCH --account=desi
#SBATCH --constraint=cpu
#SBATCH -q regular
#SBATCH -J FA456
#SBATCH -t 0:20:00
#SBATCH --array=12-80
#SBATCH --output=holi_%j.out
#SBATCH --error=holi_%j.err
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=jcolley@lpnhe.in2p3.fr

#
# Environment
#
#source /global/common/software/desi/desi_environment.sh main
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main

LSS_DIR=/global/cfs/cdirs/desi/users/colley/LSS
export PYTHONPATH=$LSS_DIR/py:$PYTHONPATH
export PATH=$LSS_DIR/bin:$PATH

PROC_DIR=$LSS_DIR/scripts/mock_tools/pipe_holi

#
# Parameters
#
seed=$(printf "seed%04d" "$SLURM_ARRAY_TASK_ID")
#seed=seed0000
DS_DIR=/global/cfs/cdirs/desi/mocks/fa4acm/holi/webjax_v4.80
# dataflow name input
in_4=imforFA0_Y3_noimagingmask_applied.fits
in_5=forFA0.fits
in_6=forFA0_withcontaminants.fits
in_7=forFA0_concat.fits

#
# Pipeline
#

echo "step 4 mask"
date
echo $DS_DIR/$seed/ELG/$in_4
# step 4 mask
# -n or --ntasks 1 else slurm/perlmutter launch 255 task by default
# -c or --cpus-per-task 1 because python script are mono-process
srun -n 1 -c 1 $PROC_DIR/join_imaging_mask_stdpars.py --inputs $DS_DIR/$seed/ELG/$in_4 --outputs $DS_DIR/$seed/ELG/$in_5 &
srun -n 1 -c 1 $PROC_DIR/join_imaging_mask_stdpars.py --inputs $DS_DIR/$seed/LRG/$in_4 --outputs $DS_DIR/$seed/LRG/$in_5 &
srun -n 1 -c 1 $PROC_DIR/join_imaging_mask_stdpars.py --inputs $DS_DIR/$seed/QSO/$in_4 --outputs $DS_DIR/$seed/QSO/$in_5 &
wait

# step 5 contaminant
echo "step 5"
date
srun -n 1 -c 1 $PROC_DIR/add_contaminants_to_mock_stdpars.py --inputs $DS_DIR/$seed/ELG/$in_5 --outputs $DS_DIR/$seed/ELG/$in_6 &
srun -n 1 -c 1 $PROC_DIR/add_contaminants_to_mock_stdpars.py --inputs $DS_DIR/$seed/QSO/$in_5 --outputs $DS_DIR/$seed/QSO/$in_6 &
wait

# step 6 concatenate tracers
echo "step 6"
date
in_ELG=$DS_DIR/$seed/ELG/$in_6
in_LRG=$DS_DIR/$seed/LRG/$in_5
in_QSO=$DS_DIR/$seed/QSO/$in_6
out_6=$(printf "forFA%04d.fits" "$SLURM_ARRAY_TASK_ID")
srun -n 1 -c 1 $PROC_DIR/concatenate_tracers_to_fba_stdpars.py --inputs $in_ELG $in_LRG $in_QSO --outputs $DS_DIR/$out_6
date
