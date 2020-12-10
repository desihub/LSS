import os
import fitsio
import numpy as np
from astropy.table import Column, Table
from astropy.coordinates import SkyCoord
from astropy import units as u

truthf = '/project/projectdirs/desi/users/ajross/MCdata/seed.fits'

ranf = '/global/cfs/cdirs/desi/target/catalogs/dr9m/0.44.0/randoms/resolve/randoms-1-0.fits'


#obiout = os.environ['obiwan_out'] #get this setup at some point


#assert(name_for_randoms is not None);assert(startid is not None);assert(nobj is not None)
topdir = '/global/cscratch1/sd/ajross/Obiwan/dr9m/obiwan_out/test/'
#topdir_tractor = os.environ['obiwan_out']+'/output/'
topdir_tractor = topdir + 'output/'
#sim_topdir = os.environ['obiwan_out']+'/divided_randoms/'
sim_topdir = '/global/cscratch1/sd/ajross/Obiwan/dr9m/obiwan_out/test/divided_randoms/'
matched_dir = topdir+'matched_obiwan/'


def mkbricklist_sampebv(nbrick=100,reg='N',ebvm=0.002,ebvx=0.15,fn='test'):
    bands = ['G','R','Z']
    kr = ['PHOTSYS','BRICKNAME','EBV','DEC']+ ['NOBS_%s' % b for b in bands]
    rall = fitsio.read(ranf,columns=kr)
    print('total # of randoms:')
    print(len(rall))
    mask = rall['PHOTSYS'] == reg
    if reg == 'S':
        mask &= (rall['DEC'] > -30)
    print(' # of randoms after restricting to ' +reg)
    print(len(rall[mask]))
    for b in bands: 
        mask &= rall['NOBS_' + b]>0
    print(' # of randoms after restricting to nobs > 0')
    print(len(rall[mask]))
    rall = rall[mask]
    print(len(rall))
    bl = []
    es = (ebvx-ebvm)/nbrick
    outf = 'bricklist_'+fn+'.txt'
    fo = open(outf,'w')
    for i in range(0,nbrick):
        we = rall['EBV'] > i*es+ebvm
        we &= rall['EBV'] < (i+1)*es+ebvm
        re = rall[we]
        if len(re) > 0:
        	ind = 0
        	bn = re[0]['BRICKNAME']
        	while np.isin(re[ind]['BRICKNAME'],bl):
        		ind += 1
        		bn = re[ind]['BRICKNAME']
        		if ind > 10:
        		    return('TOOK MORE THAN 10 interations, probably a bug')
        	print(bn)
        	bl.append(bn)
        	fo.write(bn+'\n')
        else:
            print(i*es)
    
    fo.close()   
        
        
    
    

def getran_brick(brick,maxg=24,ming=22):
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
    w &= tf['g'] > ming
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
    outf = sim_topdir+'brick_'+brick+'.fits'
    to.write(outf,format='fits',overwrite=True)
    print('wrote '+str(len(rb))+' imaging randoms in obiwan format to '+outf)
    
def SV_brick_match(brickname, name_for_randoms = 'matched_', startid = 0, nobj = 200, angle = 1.5/3600):
    #initial code copied from https://raw.githubusercontent.com/DriftingPig/Obi-Metallica/master/collect/SV_collect.py
    
    rs_type= 'rs'+str(startid)
    #print(brickname,rs_type)
    fn_tractor = os.path.join(topdir_tractor,'tractor',brickname[:3],brickname,rs_type,'tractor-%s.fits' %brickname)
    fn_sim = os.path.join(topdir_tractor,'obiwan',brickname[:3],brickname,rs_type,'simcat-elg-%s.fits' %brickname)
    fn_original_sim = sim_topdir+'/brick_'+brickname+'.fits'
    
    tractor = Table.read(fn_tractor)
    sim = Table.read(fn_sim)
    
    original_sim = Table.read(fn_original_sim)[startid:startid+nobj] 
    
    #import pdb;pdb.set_trace()
    c1 = SkyCoord(ra=sim['ra']*u.degree, dec=sim['dec']*u.degree)
    c2 = SkyCoord(ra=np.array(tractor['ra'])*u.degree, dec=np.array(tractor['dec'])*u.degree)
    c3 = SkyCoord(ra=original_sim['ra']*u.degree, dec=original_sim['dec']*u.degree)

    idx1, d2d, d3d = c1.match_to_catalog_sky(c2)
    idx2, d2d2, d3d2 = c1.match_to_catalog_sky(c3)

    matched = d2d.value <= angle
    distance = d2d.value
    tc = tractor[idx1]

    ors = original_sim[idx2]
    
    tc.add_column(sim['gflux'],name = 'sim_gflux')
    tc.add_column(sim['rflux'],name='sim_rflux')
    tc.add_column(sim['zflux'],name='sim_zflux')
    tc.add_column(ors['redshift'],name='sim_redshift')
    tc.add_column(ors['id'],name='TARGETID')
    tc.add_column(sim['rhalf'],name='sim_rhalf')
    tc.add_column(sim['e1'],name='sim_e1')
    tc.add_column(sim['e2'],name='sim_e2')
    tc.add_column(sim['x'],name='sim_bx')
    tc.add_column(sim['y'],name='sim_by')
    
    tc['detected'] = np.array(matched,dtype=np.bool)
    tc.add_column(sim['n'],name='sim_sersic_n')
    outf = matched_dir+name_for_randoms+brickname+'.fits'
    tc.write(outf,format='fits',overwrite=True)
    print('wrote matched output to '+outf)
    #return tc

    
     