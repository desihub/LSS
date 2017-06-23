"""
Lado Samushia: lado@phys.ksu.edu February 2016

This script will take tiles files and a secret file resulting from a fiber
assignment algorithm and will create a large scale structure catalogue
containing all objects. 

Usage:
FAtoCat.py $SCRATCH/tiles*.txt $SCRATCH/secret.txt catalogue.txt

catalogue.txt is an ASCII file with columns:
target_id, ra, dec, redshift, target_type, observed

observed is either 1 or 0 depending on whether the target was observed or not.
Since at the moment the code is run on simulations the correct redshift is
printed anyway.

For the format of tiles*.txt and secret.txt see
https://github.com/desihub/fiberassign/blob/master/doc/new_instructions.txt
"""
import sys
import numpy as np

catalogue_file = sys.argv[-1]
secret_file = sys.argv[-2]
tile_files = sys.argv[1:-2]

fcat = open(catalogue_file,'w')
fsec = open(secret_file,'r')

ra = []
dec = []
red = []
type = []

for line in fsec:
    columns = line.split()
    ra.append(columns[2])
    dec.append(columns[3])
    red.append(columns[4])
    type.append(columns[5])
fsec.close()

tar_size = len(type)
print(tar_size, len(red), len(dec), len(ra))
#Is the target already in the catalogue?
duplicated = np.zeros(tar_size, dtype=bool)

for file in tile_files:
    #print(file)
    with open(file,'r') as f:
        for line in f:
            columns = line.split()
            #-1 means the fiber was not assigned
            if columns[-1] == "-1":
                continue
            id = int(columns[-5])
            tt = int(type[id])
            #If the target is not QSOa, QSO, LRG, or ELG, 
            #or it is already in the catalogue
            if not(tt == 0 or tt == 1 or tt == 2 or tt == 3) or duplicated[id]:
                continue
            duplicated[id] = True
            ra_id = columns[-4]
            dec_id = columns[-3]
            fcat.write('%d %s %s %s %d 1\n'%(id, ra_id, dec_id, red[id], tt))
        
#Targets that don't have observed redshifts
for i in range(tar_size):
#    print(i)
    if duplicated[i] != True:
#        print(i)
#        print(duplicated[i])
#        print(ra[i])
#        print(dec[i])
#        print(red[i])
#        print(type[i])
        fcat.write('%d %s %s %s %s 0\n'%(i, ra[i], dec[i], red[i], type[i]))

fcat.close()
