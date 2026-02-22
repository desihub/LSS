'''
Tihis script will read the information about how number density varies across the
bricks as a function of various parameters and will create a random catalogue
that describes the angular mask.

At the moment the fit is just a multidimensional linear fit.

The parameters currently considered are: depth_r, galdepth_r, extinction.

Lado Samushia (colkhis@gmail.com), May 2017
'''

from astropy.io import fits
from astropy.table import Table
import numpy as np
import argparse
from itertools import izip
import os

parser = argparse.ArgumentParser(description='Construct angular mask for targeting')
parser.add_argument('in_file',help='Input file with bricks information')
parser.add_argument('out_file',help='Output file with ra and dec of random points')
parser.add_argument('Nran',help='How many random files?')
args = parser.parse_args()

if os.path.isfile(args.out_file) == True:
    print("The random files with this name already exist. Can not overwrite.")
    print("Choose an alternative name for the output random catalogue files.")
    quit()


# The columns are (current format):
# Ngal, depth_r, gdepth_r, ebv, area, ramin, ramax, decmin, decmax
bricks = np.load(args.in_file)
area = bricks[4,:]

# Compute angular number density 
bricks[0,:] /= area
# Convert radians to degrees
bricks[5:9,:] = np.degrees(bricks[5:9,:])

# Correlation between number density and other parameters
CC = np.cov(bricks[:4,:],rowvar=1)
print(CC)

# Mean number density and parameters
MM = np.mean(bricks[:4,:],axis=1)
print(MM)

# I am solving the equation nth - nave = sum(acoeff_i*(par - parave))
a_coeff = np.dot(CC[0,1:],np.linalg.inv(CC[1:,1:]))
print(np.shape(a_coeff))
print(a_coeff)
mean_par = np.outer(np.mean(bricks[1:4,:],axis=1),np.ones(67698))
nave_th = MM[0] + np.dot(a_coeff, bricks[1:4,:] - mean_par)
Nave_th = nave_th*area

print("a_coeff",a_coeff)
print("bricks",bricks[1:4,:])
print("mean_par",mean_par)
print("nave_th",nave_th)

for i in range(int(args.Nran)):
# Generate random ra, dec inside the brick based on nave_th
    print("generating ra")
    ra = []
    for (min,max,num) in izip(bricks[5,:],bricks[6,:],Nave_th):
        ra.extend(np.random.uniform(min,max,num).tolist())

    print("generating dec")
    dec = []
    for (min,max,num) in izip(bricks[7,:],bricks[8,:],Nave_th):
        dec.extend(np.random.uniform(min,max,num).tolist())

    np.asarray(ra, dtype=np.float32)
    np.asarray(dec, dtype=np.float32)

    radec = Table([ra,dec],names=('ra','dec'))
    radec.write(args.out_file + "_" + str(i) + '.fits')
