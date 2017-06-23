import os
from astropy.io import fits
from astropy.table import Table, vstack
from astropy.table import Column
import numpy as np

i=0
j=0
inin=0
file_no = 8
for filename in os.listdir(os.getcwd()):
    inin +=1  

    tile_data = fits.open(filename)[1].data
    mask = tile_data['DESI_TARGET']==2
    add_rows = Table(tile_data[mask])
    if i==0:
        obs2 = add_rows
    else:
        obs2 =(vstack([obs2, add_rows]))
    i+=1
    if i %100 ==0:
        outfile = "ELG_%d_%d.csv" % (file_no, j)   
        obs2.write(outfile, format='csv')
        i=0
        j+=1
outfile = "ELG_%d_100.csv" % file_no
print(inin)
obs2.write(outfile, format='csv')

