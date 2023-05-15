import os
import pandas as pd
from astropy.table import Table
from astropy.io import fits
import numpy as np
import glob
#This script modifies SV3 products to meet the datamodel standards. It substitutes the files indicated in the specified directory, 
#Therefore, be careful. Suggestion is to work on a copy of the original files.

#This is the directory where we make the modifications. 
dir_ = '/global/cfs/cdirs/desi/survey/catalogs/edr_prepfor_public/fuji/sv3/v2/LSScats/clustering'

for root, dirs, files in os.walk(dir_):
    for file in files:
        if file.endswith(".fits") and file.startswith("BGS"):
            print(file)
            filename = os.path.join(root,file)
            hdul = fits.open(filename,mode='update')
            data_to_read = hdul[1]
            data_to_read.data['EQ_ALL_0P1'] = -0.97*(data_to_read.data['Z']-0.1)
            data_to_read.data['ABSMAG_R'] = data_to_read.data['ABSMAG_R'] + data_to_read.data['EQ_ALL_0P0']

            hdul.close()
        else:
            pass



