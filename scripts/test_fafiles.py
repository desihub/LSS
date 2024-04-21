import numpy as np
import fitsio
from glob import glob

fls = glob('/dvs_ro/cfs/cdirs/desi/survey/fiberassign/main/*/*.fits')

for fl in fls:
	try:
		f = fitsio.read(fl,rows=1)
	except:
	    print(fl)