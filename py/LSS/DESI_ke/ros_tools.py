import numpy as np


# https://github.com/desihub/LSS/blob/3d65b099a763a41179cb8e706f85482b3633a254/py/LSS/SV3/cattools.py#L21  
@np.vectorize
def tile2rosette(tile):
    if tile < 433:
        return (tile-1)//27
    
    else:
        if tile >= 433 and tile < 436:
            return 13
        if tile >= 436 and tile < 439:
            return 14
        if tile >= 439 and tile < 442:
            return 15
        if tile >= 442 and tile <=480:
            return (tile-442)//3
            
        if tile > 480:
            return tile//30    

    return -999999

roscen = {0:(150.100,2.182),\
          1:(179.6,0),\
          2:(183.1,0),\
          3:(189.9,61.8),\
          4:(194.75,28.2),\
          5:(210.0,5.0),\
          6:(215.5,52.5),\
          7:(217.8,34.4),\
          8:(216.3,-0.6),\
          9:(219.8,-0.6),\
          10:(218.05,2.43),\
          11:(242.75,54.98),\
          12:(241.05,43.45),\
          13:(245.88,43.45),\
          14:(252.5,34.5),\
          15:(269.73,66.02),\
          16:(194.75,24.7),\
          17:(212.8,-0.6),\
          18:(269.73,62.52),\
          19:(236.1,43.45)}

def calc_rosr(rosn, ra, dec):
    # Given rosetter number and ra,dec, calculate distance from center 
    ra       = ra*np.pi/180.
    dec      = dec*np.pi/180.
    rac,decc = roscen[rosn]
    rac      = rac*np.pi/180.
    decc     = decc*np.pi/180.
    cd       = np.sin(dec)*np.sin(decc)+np.cos(dec)*np.cos(decc)*np.cos(rac-ra)
    ad       = np.arccos(cd)*180./np.pi

    return  ad

def ros_limits(dryrun):
    if dryrun:
        limits                = [0.9, 1.10]
    else:
        limits                = [0.2, 1.75]
        
    return limits
