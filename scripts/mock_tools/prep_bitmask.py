import os

def mk_inputandoutput_fn(file_root,qsodir='/global/cfs/cdirs/desi/mocks/cai/holi/webjax_v4.80/',lrgdir='/global/cfs/cdirs/desi/mocks/cai/holi/webjax_v4.80/',elgdir='/global/cfs/cdirs/desi/mocks/cai/holi/webjax_v4.81/',qsomin=0,qsomax=0,lrgmin=0,lrgmax=0,elgmin=0,elgmax=0):
    fn_inputs = open(file_root+'_in.txt','w')
    fn_outputs = open(file_root+'_out.txt','w')
    for i in range(qsomin,qsomax):
        in_fn = qsodir+'seed'+str(i).zfill(4)+'/QSO/forFA0_Y3_noimagingmask_applied.fits'
        if os.path.isfile(in_fn):
            fn_inputs.write(in_fn+'\n')
            out_fn = qsodir+'seed'+str(i).zfill(4)+'/QSO/imforFA0_Y3_noimagingmask_applied.fits'
            fn_outputs.write(fn_outputs+'\n')
        else:
            print(in_fn+' not found')
    for i in range(lrgmin,lrgmax):
        in_fn = lrgdir+'seed'+str(i).zfill(4)+'/LRG/forFA0_Y3_noimagingmask_applied.fits'
        if os.path.isfile(in_fn):
            fn_inputs.write(in_fn+'\n')
            out_fn = lrgdir+'seed'+str(i).zfill(4)+'/LRG/imforFA0_Y3_noimagingmask_applied.fits'
            fn_outputs.write(fn_outputs+'\n')
        else:
            print(in_fn+' not found')
    for i in range(elgmin,elgmax):
        in_fn = elgdir+'seed'+str(i).zfill(4)+'/ELG/forFA0_Y3_noimagingmask_applied.fits'
        if os.path.isfile(in_fn):
            fn_inputs.write(in_fn+'\n')
            out_fn = lrgdir+'seed'+str(i).zfill(4)+'/ELG/imforFA0_Y3_noimagingmask_applied.fits'
            fn_outputs.write(fn_outputs+'\n')
        else:
            print(in_fn+' not found')

    fn_inputs.close()
    fn_outputs.close()
    
    

def mk_confile(conf_fn,file_root):
	input_fn = file_root+'_in.txt'
	output_fn = file_root+'_out.txt'
	fo = open(file_root+conf_fn,'w')
	output = '''# Configuration file for BRICKMASK (default: `brickmask.conf').
# Format: keyword = value # comment
#     or: keyword = [element1, element2]
#    see: https://github.com/cheng-zhao/libcfg for details.
# Some of the entries allow expressions, see
#         https://github.com/cheng-zhao/libast for details.
# NOTE that command line options have priority over this file.
# Unnecessary entries can be left unset.

BRICK_LIST      = /dvs_ro/cfs/cdirs/cosmo/data/legacysurvey/dr9/randoms/survey-bricks-dr9-randoms-0.48.0.fits
	# Filename for the FITS table with the list of all bricks
	# and their regions (N or S), see
	# www.legacysurvey.org/dr9/files/#survey-bricks-dr9-randoms-0-48-0-fits
BRICK_PATH      = /dvs_ro/cfs/cdirs/cosmo/data/legacysurvey/dr9
	# String, top level directory of maskbit and nexp files, see
	# https://www.legacysurvey.org/dr9/files/#for-web-access
MASKBIT_NULL    = 1
	# Integer, bit code for objects outside all maskbit bricks (unset: 1).
NEXP_NULL       = 0
	# Integer, bit code for objects without exposure counts (unset: 0).

INPUT_FILES     = ''' + input_fn + '''
	# Filename of an ASCII file storing paths of input catalogs.
	# Formats and columns of all input files must be identical.
	# Each row of the ASCII file specifies the path of an input catalog.
	# Each space in the paths must be escaped by a leading '\' character.
	# Lines starting with '#' are omitted.
FILE_TYPE       = 1
	# Integer, format of the input and output catalogs (default: 0).
	# The allowed values are:
	# * 0: ASCII text file;
	# * 1: FITS table.
ASCII_COMMENT   = 
	# Character indicating comment lines for ASCII-format catalog (unset: '').
COORD_COLUMN    = ["RA", "DEC"]
	# 2-element integer or string array, columns of (RA,Dec) for `INPUT`.
	# They must be integers indicating the column numbers (starting from 1) for
	# an ASCII file, or strings indicating the column names for a FITS table.
OUTPUT_FILES    = ''' + output_fn + '''
	# Filename of an ASCII file storing paths of output catalogs.
	# Each row of the ASCII file specifies the path of an output catalog that
	# corresponds to the input catalog in `INPUT_FILES` at the same row.
	# Each space in the paths must be escaped by a leading '\' character.
	# Lines starting with '#' are omitted.
OUTPUT_COLUMN   = 
	# Integer or String arrays, columns to be saved to `OUTPUT`.
	# If not set, all columns of `INPUT` are saved in the original order.
	# Note that maskbits (and optionally subsample IDs) are always saved
	# as the last column (or last two columns).
MASKBIT_COLUMN  = "MASKBITS"
	# String, name of the maskbit column in the FITS-format `OUTPUT`.
NEXP_COLUMNS    = ["NOBS_G", "NOBS_R", "NOBS_Z"]
	# String array, names of the nexp columns in the FITS-format `OUTPUT`.
OVERWRITE       = 1
	# Flag indicating whether to overwrite existing files, integer (unset: 0).
	# Allowed values are:
	# * 0: quit the program when an output file exist;
	# * positive: force overwriting output files whenever possible;
	# * negative: notify at most this number of times for existing files.
VERBOSE         = 
	# Boolean option, indicate whether to show detailed outputs (unset: T).
'''
	fo.write(output)
	fo.close()

def mk_sbatch(batch_fn,conf_file,rootdir='/global/homes/d/desica/BRICKMASKcode/brickmask/'):
    fo = open(batch_fn,'w')
    #conf_file = rootdir+conf_file
    outs = '''#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH -J holiv3_1
#SBATCH --constraint=cpu
#SBATCH --account=desi
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=64
#SBATCH  --dependency=afterany:46580373
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main

module load cpu cray-fftw
export CFITSIO_DIR=/global/common/software/desi/users/naimgk/cfitsio
export LD_LIBRARY_PATH=$CFITSIO_DIR/lib:$LD_LIBRARY_PATH

srun -n 64 -c 4 --cpu-bind=cores '''+rootdir+'/BRICKMASK -c '''+conf_file+'\n'
    fo.close()


rundir = '/global/homes/d/desica/BRICKMASKcode/brickmask/'
qsomin = 100
qsomax = 125
fileroot = rundir+'holi_AJRrun1'
mk_inputandoutput_fn(file_root=fileroot,qsomin=qsomin,qsomax=qsomax,elgmin=0,elgmax=50)
mk_confile('.conf',fileroot)
mk_sbatch(fileroot+'.sbatch',fileroot+'.conf',rundir)
