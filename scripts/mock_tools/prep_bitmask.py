def mk_confile(conf_fn,input_fn,output_fn):
	fo = open(conf_fn,'w')
	output = """
	# Configuration file for BRICKMASK (default: `brickmask.conf').
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
	
	INPUT_FILES     = """+input_fn+'\n'
	+"""
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
	OUTPUT_FILES    = """+output_fn+'\n'+
	"""
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
	"""
	fo.write(output)
	fo.close()