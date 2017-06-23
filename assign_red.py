from astropy.table import Table, join
from astropy.io import fits, ascii
import numpy as np
from random import shuffle

file_index =6
infile = '/scratch1/scratchdirs/angela/rand%db/radec.csv' % file_index
outfile = '/scratch1/scratchdirs/angela/rand%db/ELG_FArands_mtl0_%d.fits' % (file_index, file_index)
data = ascii.read("/global/u1/a/angela/DESI/ELG_nz_zcat0.csv", format='csv')
radec= ascii.read(infile, format='no_header')
num_lines = len(radec)
# num_lines = sum(1 for line in open(infile))
number_0f_ELG = num_lines*data['number']
red_max = data['red_max']
red_min = data['red_min']
number_0f_ELG = np.around(number_0f_ELG)
z_vector = np.random.random(num_lines)
cum_no_ELG = np.cumsum(number_0f_ELG)
val = int(cum_no_ELG[0])
z_vector[:val] = z_vector [:val]*(red_max[0]-red_min[0])+ red_min[0]
for i in range(1,len(cum_no_ELG)):
    z_vector[int(cum_no_ELG[i-1]):int(cum_no_ELG[i])] = z_vector[int(cum_no_ELG[i-1]):int(cum_no_ELG[i])] *(red_max[i]-red_min[i])+ red_min[i]
j=0
int_val = int(cum_no_ELG[-1])
if cum_no_ELG[-1] < num_lines:
    extra_lines = num_lines-cum_no_ELG
    for i in range(int_val,num_lines):
        z_vector[i] = red_min[j]
        j+=1

new_z_vect = shuffle(z_vector)
t = Table ([radec['col1'], radec['col2'],new_z_vect], names=['RA', 'DEC', 'Z'])
hdu = fits.table_to_hdu(t)
hdu.writeto(outfile)

