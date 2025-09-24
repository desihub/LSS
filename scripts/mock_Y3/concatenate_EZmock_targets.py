#!/usr/bin/env python
# coding: utf-8

# Check whether the input catalog has applied Y3 footprint mask. If not, apply it here.

import sys
import numpy as np
from astropy.table import Table, vstack
import fitsio


seed = int(sys.argv[1])
idir = sys.argv[2]
odir = sys.argv[3]

tracers = ['LRG', 'ELG', 'QSO']


#idir0 = "/pscratch/sd/j/jerryou/DESI_Y3/SecondGenMocks/EZmock/{tracer}/"
idir0 = idir+"/{tracer}/"
filename = "EZmock_{tracer}_DA2_c000_ph000_NScomb_{seed:04d}.fits.gz"



output = Table()
for tracer in tracers:
    ifile = (idir0 + filename).format(tracer=tracer, seed=seed)
    print(ifile)
    cat = Table.read(ifile)
    if np.sum((cat['STATUS']&2)>0) < len(cat):
        print("apply Y3 footprint mask.")
        y3_mask = (cat['STATUS']&2)>0
        cat = cat[y3_mask]
    output = vstack((output, cat))    


print("Targets:", np.unique(output['TARGETID']))


#odir = "/pscratch/sd/j/jerryou/DESI_Y3/SecondGenMocks/EZmock/forFA/"
ofile = odir + f"forFA{seed}.fits"
output.write(ofile, overwrite=True)




