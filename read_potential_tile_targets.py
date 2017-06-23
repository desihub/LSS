import os
from astropy.io import fits
from astropy.table import Table, vstack
from astropy.table import Column, unique
import numpy as np
i=0
j=0
for filename in os.listdir(os.getcwd()):
    tile_data = fits.open(filename)[2].data
    add_rows_all = Table(tile_data)
    add_rows = unique(add_rows_all,keys='POTENTIALTARGETID')
    if i==0:
        obs2 = add_rows
    else:
        obs2 =(vstack([obs2, add_rows]))
    i+=1

rows_unique = unique(obs2,keys='POTENTIALTARGETID')
        

