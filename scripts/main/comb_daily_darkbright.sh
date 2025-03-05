#!/bin/bash
#script to run both dark and bright combination steps (dark might have to be resubmitted to finish it)
#set environment externally like described here: https://desi.lbl.gov/trac/wiki/SurveyOps/LSSdaily

srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSSCODE/LSS/scripts/main/combdata_main.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec daily --prog dark 
srun -N 1 -C cpu -t 04:00:00 --qos interactive --account desi python $LSSCODE/LSS/scripts/main/combdata_main.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --verspec daily --prog bright 
