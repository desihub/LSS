from __future__ import print_function
from astropy.io import fits 
from astropy.coordinates import SkyCoord 
from astropy import units as u                               
from itertools import izip                                                      
import numpy as np                                                              
import sys                                                                      
import desispec.brick                                                           
                                                                                
# Assign a unique integer to each brick                                         
def bn_to_unique_int(name):                                                     
    hs = name[4]                                                                
    if hs == "p":                                                               
        return int(name[:4] + "1" + name[5:])                                   
    else:                                                                       
        return int(name[:4] + "0" + name[5:])                                   
                                                                                
targets = fits.open("../../quicksurvey2016/input/dark/targets.fits")            
print("Loading targets")                                                        
# Target information                                                            
tbname = targets[1].data.field("brickname")                                     
tdesi = targets[1].data.field("desi_target")                                    
tdepthr = targets[1].data.field("depth_r")                                      
tgdepthr = targets[1].data.field("galdepth_r")                                  
tra = targets[1].data.field("ra")                                               
tdec = targets[1].data.field("dec")                                             
targets.close()                                                                 
                                                                                
# Do this for ELGs only for now (2 means ELG)                                   
print("select ELGs")
tbname = tbname[tdesi == 2]                                                     
tdepthr = tdepthr[tdesi == 2]                                                    
tgdepthr = tgdepthr[tdesi == 2]                                                  
tra = tra[tdesi == 2]                                                           
tdec = tdec[tdesi == 2]                                                         
                                                                                
# Unique bricknames                                                             
print("Map bricknames to integers")
tbname = map(bn_to_unique_int,tbname)                                           
print("Set unique bricks")
ubname = set(tbname)                                                            
tbname = np.asarray(tbname,dtype=np.int32) 
                                                                   
# Finding EBV for each target                                                   
tebv = np.zeros((np.size(tbname)))
print("Find extinction for each brick")
ebv_north = fits.open("/project/projectdirs/desi/software/edison/dust/v0_1/maps/SFD_dust_4096_ngp.fits")
ebv_south = fits.open("/project/projectdirs/desi/software/edison/dust/v0_1/maps/SFD_dust_4096_sgp.fits")
SCALE = ebv_north[0].header["lam_scal"]                                         
ebv_N = ebv_north[0].data                                                       
ebv_S = ebv_south[0].data                             
# Transform to galactic coordinates
print("Transforming coordinates to Galactic")
c_icrs = SkyCoord(ra=tra*u.degree, dec=tdec*u.degree, frame='icrs')             
l = c_icrs.galactic.l.radian                                                    
l = np.asarray(l,dtype=np.float64)                                              
b = c_icrs.galactic.b.radian                                                    
b = np.asarray(b,dtype=np.float64)                                              
NorS = np.sign(b)                                                               
print("Finding indeces")
X = np.sqrt(1-NorS*np.sin(b))*np.cos(l)*SCALE + SCALE                           
X = X.astype(int)                                                               
Y = NorS*np.sqrt(1-NorS*np.sin(b))*np.sin(l)*(-SCALE) + SCALE                   
Y = Y.astype(int)                                                               
print("Setting EBV")
tebv[b > 0] = ebv_N[X[b>0],Y[b>0]]                                              
tebv[b < 0] = ebv_S[X[b<0],Y[b<0]]                                              
                                                                                
# Brick area, ra, dec (Based on Stephen's hack)                                 
print("Dictionary for brick properties")
b = desispec.brick.Bricks(0.5)                                                  
dbrick_info = dict()                                                            
for row in range(len(b._brickname)):                                            
    decmin, decmax = np.radians(b._edges_dec[row:row+2])                        
    for col in range(len(b._brickname[row])):                                   
        ramin, ramax = np.radians(b._edges_ra[row][col:col+2])                  
        area = (ramax-ramin)*(np.sin(decmax) - np.sin(decmin))                  
        area *= (180/np.pi)**2                                                  
        brick_info = (area, ramin, ramax, decmin, decmax)                       
        dbrick_info[bn_to_unique_int(b._brickname[row][col])] = brick_info      
                                                                                
#                                                                               
udepthr = np.zeros(len(ubname))                                             
ugdepthr = np.zeros(len(ubname))                                            
uebv = np.zeros(len(ubname))                                                
ungal = np.zeros(len(ubname))                                               
area = np.zeros(len(ubname))
ramin = np.zeros(len(ubname))
decmin = np.zeros(len(ubname))
ramax = np.zeros(len(ubname))
decmax = np.zeros(len(ubname))

print("Loop over bricks")                                                                                
counter = 0
for ind, name in enumerate(ubname):                                             
    mask = (tbname == name)                                                     
    udepthr[ind] = np.mean(tdepthr[mask])                                       
    ugdepthr[ind] = np.mean(tgdepthr[mask])                                     
    uebv[ind] = np.mean(tebv[mask])                                             
    ungal[ind] = np.sum(mask)                                                   
    ubrick_info = dbrick_info[name]                                             
    area[ind], ramin[ind], ramax[ind], decmin[ind], decmax[ind] = ubrick_info   

sbrick = np.array([ungal,udepthr,ugdepthr,uebv,area,ramin,ramax,decmin,decmax])

np.save("bm_imaging.npy",sbrick)
