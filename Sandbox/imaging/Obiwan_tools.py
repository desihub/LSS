import fitsio
import numpy as np
from astropy.table import Table

truthf = '/project/projectdirs/desi/users/ajross/MCdata/seed.fits'

ranf = '/global/cscratch1/sd/adamyers/dr9m-sep26-2020/0.42.0/randoms/resolve/randoms-1-0.fits'
#file has 63283969 entries

#obiout = os.environ['obiwan_out'] #get this setup at some point

outdir = '/global/cscratch1/sd/ajross/Obiwan/dr9m/obiwan_out/test/divided_randoms/'

def getran_brick(brick,maxg=25):
    #columns to consider from imaging randoms
    kr = ['RA','DEC','TARGETID','BRICKNAME'] #minimum number needed for now
    #kr = ['RELEASE','BRICKID','BRICKNAME','RA','DEC','NOBS_G','NOBS_R','NOBS_Z','PSFDEPTH_G','PSFDEPTH_R','PSFDEPTH_Z','GALDEPTH_G','GALDEPTH_R',\
 #'GALDEPTH_Z','PSFDEPTH_W1','PSFDEPTH_W2','PSFSIZE_G','PSFSIZE_R','PSFSIZE_Z','APFLUX_G','APFLUX_R','APFLUX_Z','APFLUX_IVAR_G','APFLUX_IVAR_R',\
 #'APFLUX_IVAR_Z','MASKBITS','WISEMASK_W1','WISEMASK_W2','EBV','PHOTSYS','TARGETID','HPXPIXEL']
    #Load imaging randoms
    rall = fitsio.read(ranf,columns=kr)
    #select brick in question
    w = rall['BRICKNAME'] == brick
    rb = rall[w]
    del rall
    #columns to use from truth file
    knames = ['objid','type','g', 'r','z','w1','w2','hsc_mizuki_photoz_best','rhalf']
    tf = fitsio.read(truthf,columns=knames)
    #cut faint objects from truth
    w = tf['g'] < maxg
    tf = tf[w]
    #select a random selection form the truth data
    tb = np.random.choice(tf,len(rb))
    to = Table()
    to['id'] = rb['TARGETID']
    to['ra'] = rb['RA']
    to['dec'] = rb['DEC']
    to['id_sample'] = tb['objid']
    to['redshift'] = tb['hsc_mizuki_photoz_best']
    to['rhalf'] = tb['rhalf']
    to['g'] = tb['g']
    to['r'] = tb['r']
    to['z'] = tb['z']
    ntype = np.ones(len(tb['type']),dtype=int)
    wt = tb['type'] == 'DEV'
    ntype[wt] = 4
    to['n'] = ntype
    to['ba'] = np.random.uniform(0.2,1.,size=len(rb))
    to['pa'] = np.random.uniform(0,180.,size=len(rb))
    outf = outdir+'brick_'+brick+'.fits'
    to.write(outf,format='fits',overwrite=True)
    
     