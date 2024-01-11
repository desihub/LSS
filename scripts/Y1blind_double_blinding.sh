#!/bin/bash

# Script for applying double blinding to DESI data

# -------------- Environment Setup --------------
# Load DESI environment
source /global/common/software/desi/desi_environment.sh main
# Load custom DESI environment for cosmodesi
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
# Swap to MPI-enabled pyrecon version
module swap pyrecon/main pyrecon/mpi
# Add custom python scripts to PYTHONPATH
PYTHONPATH=$PYTHONPATH:$HOME/LSS/py 

# -------------- Configuration Variables --------------
# Define tracer, survey, version, and directories
TRACER='LRG'
SURVEY='Y1'
VERSPEC='iron'
VER='test/blinded'
INDIR='/global/cfs/cdirs/desi/survey/catalogs/'
OUTDIR='/global/cfs/cdirs/desi/users/uendert/desi_blinding/double_blinded/test_mkclusdat_y_mkclusran_y/'

# Define cosmology parameters for double blinding
w0=-0.90
wa=0.26
fnl=20.0

# Define the name of the log file based on timestamp to ensure uniqueness
LOGFILE="blinding_$(date +%Y%m%d_%H%M%S).log"
# Create the log file
touch "$LOGFILE"

# -------------- Main Processing Function --------------
run_blinding() {
    srun -n 128 -N 1 -C cpu -t 04:00:00 -q interactive python scripts/main/apply_blinding_main_fromfile_fcomp_double_blinding.py \
        --type $TRACER \
        --survey $SURVEY \
        --basedir_in $INDIR \
        --basedir_out $OUTDIR \
        --verspec $VERSPEC \
        --version $VER \
        --baoblind y \
        --mkclusdat y \
        --mkclusran y \
        --resamp y \
        --maxr 18 \
        --dorecon y \
        --rsdblind y \
        --fnlblind y \
        --getFKP y \
        --specified_w0 $w0 \
        --specified_wa $wa \
        --specified_fnl $fnl \
        --get_par_mode specified >> "$LOGFILE"
}

# Save the settings to the log file before running the main function
echo "Running blinding process with the following settings:" | tee -a "$LOGFILE"
echo "Tracer: $TRACER" | tee -a "$LOGFILE"
echo "Survey: $SURVEY" | tee -a "$LOGFILE"
echo "Verversion spec: $VERSPEC" | tee -a "$LOGFILE"
echo "Version: $VER" | tee -a "$LOGFILE"
echo "Input Directory: $INDIR" | tee -a "$LOGFILE"
echo "Output Directory: $OUTDIR" | tee -a "$LOGFILE"
echo "Cosmology parameters - w0: $w0, wa: $wa, fnl: $fnl" | tee -a "$LOGFILE"
echo "Starting at: $(date)" | tee -a "$LOGFILE"

# Run the main process
run_blinding

# Output completion to log file
echo "Blinding process completed at: $(date)" | tee -a "$LOGFILE"
