import numpy as np
import fitsio
from astropy.table import Table,join,vstack
import datetime
import os
import sys
import logging

ext_coeff = {'G':3.214, 'R':2.165,'Z':1.211,'W1':0.184,'W2':0.113}

from LSS.tabulated_cosmo import TabulatedDESI
cosmo = TabulatedDESI()
dis_dc = cosmo.comoving_radial_distance

def dl(z):   # Luminosity distance from now to z
    return dis_dc(z)*(1.+z)

def dm(z):
    return 5.*np.log10(dl(z)) + 25.

def radec2thphi(ra,dec):
    return (-dec+90.)*np.pi/180.,ra*np.pi/180.
    
def thphi2radec(theta,phi):
    return 180./np.pi*phi,-(180./np.pi*theta-90)
    
def mask_bad_fibers_time_dependent(dz, badfibers_td):
	# Use the time-dependent bad fiber file to remove data
	# dz is the input data file, which must have FIBER and LASTNIGHT columns

	#badfibers_td = open('/global/cfs/cdirs/desi/survey/catalogs/DA2/LSS/loa-v1/unique_badfibers_time-dependent.txt','r').readlines()

	inds_to_remove = np.array([])
	for i in range(len(badfibers_td)):
		print(i)
	
		if len(badfibers_td[i].split()) == 1:
			inds_to_remove = np.concatenate((inds_to_remove, (np.where(dz['FIBER'] == int(badfibers_td[i]))[-1])))
		elif len(badfibers_td[i].split()) == 2:
			inds_to_remove = np.concatenate((inds_to_remove, (np.where((dz['FIBER'] == int(badfibers_td[i].split()[0]))
				& (dz['LASTNIGHT'] >= int(badfibers_td[i].split()[1]))
				)[-1])))
			print(np.where((dz['FIBER'] == int(badfibers_td[i].split()[0]))
				& (dz['LASTNIGHT'] >= int(badfibers_td[i].split()[1]))))
			#print(5/0)
		elif len(badfibers_td[i].split()) == 3:
			inds_to_remove = np.concatenate((inds_to_remove, (np.where((dz['FIBER'] == int(badfibers_td[i].split()[0]))
				& (dz['LASTNIGHT'] >= int(badfibers_td[i].split()[1])) & 
				(dz['LASTNIGHT'] < int(badfibers_td[i].split()[2]))
				)[-1])))
			print((np.where((dz['FIBER'] == int(badfibers_td[i].split()[0]))
				& (dz['LASTNIGHT'] >= int(badfibers_td[i].split()[1])) & 
				(dz['LASTNIGHT'] < int(badfibers_td[i].split()[2]))
				)))
		elif len(badfibers_td[i].split()) == 4:
			inds_to_remove = np.concatenate((inds_to_remove, (np.where((dz['FIBER'] == int(badfibers_td[i].split()[0]))
				& (((dz['LASTNIGHT'] >= int(badfibers_td[i].split()[1])) & 
				(dz['LASTNIGHT'] < int(badfibers_td[i].split()[2])))
				| (dz['LASTNIGHT'] >= int(badfibers_td[i].split()[3]))
				))[-1])))
			print(np.where((dz['FIBER'] == int(badfibers_td[i].split()[0]))
				& ((((dz['LASTNIGHT'] >= int(badfibers_td[i].split()[1])) & 
				(dz['LASTNIGHT'] < int(badfibers_td[i].split()[2])))
				| (dz['LASTNIGHT'] >= int(badfibers_td[i].split()[3]))))
				))
		elif len(badfibers_td[i].split()) == 5:
			inds_to_remove = np.concatenate((inds_to_remove, (np.where((dz['FIBER'] == int(badfibers_td[i].split()[0]))
				& (((dz['LASTNIGHT'] >= int(badfibers_td[i].split()[1])) & 
				(dz['LASTNIGHT'] < int(badfibers_td[i].split()[2])))
				| ((dz['LASTNIGHT'] >= int(badfibers_td[i].split()[3]))
				& (dz['LASTNIGHT'] < int(badfibers_td[i].split()[4])))
				))[-1])))
			print((np.where((dz['FIBER'] == int(badfibers_td[i].split()[0]))
				& ((dz['LASTNIGHT'] >= int(badfibers_td[i].split()[1])) & 
				(dz['LASTNIGHT'] < int(badfibers_td[i].split()[2])))
				| ((dz['LASTNIGHT'] >= int(badfibers_td[i].split()[3]))
				& (dz['LASTNIGHT'] < int(badfibers_td[i].split()[4])))
				)))
	dz_inds = np.ones(len(dz)).astype('int')
	dz_inds[inds_to_remove.astype('int')] = 0
	return dz[dz_inds == 1]



#functions that shouldn't have any dependence on survey go here

def cut_specdat(dz,badfib=None,tsnr_min=0,tsnr_col='TSNR2_ELG',logger=None,fibstatusbits=None,remove_badfiber_spike_nz=False):
    from desitarget.targetmask import zwarn_mask
    selz = dz['ZWARN'] != 999999
    selz &= dz['ZWARN']*0 == 0 #just in case of nans
    fs = dz[selz]

    #first, need to find locations to veto based data
    nodata = fs["ZWARN_MTL"] & zwarn_mask["NODATA"] != 0
    num_nod = np.sum(nodata)
    printlog('number with no data '+str(num_nod),logger)
    badqa = fs["ZWARN_MTL"] & zwarn_mask.mask("BAD_SPECQA|BAD_PETALQA") != 0
    num_badqa = np.sum(badqa)
    printlog('number with bad qa '+str(num_badqa),logger)
    nomtl = nodata | badqa
    wfqa = ~nomtl
    if tsnr_min > 0:
        low_tsnr = dz[tsnr_col] < tsnr_min
        wfqa &= ~low_tsnr
        printlog('number at low tsnr2 '+str(sum(low_tsnr)),logger)
    if fibstatusbits is not None:
        bfs = np.zeros(len(dz),dtype=bool)
        for bit in fibstatusbits:
            bfs |= (dz['COADD_FIBERSTATUS'] & 2**bit) > 0
        wfqa &= ~bfs

    #veto fibers later determined to have poor success rates
    if badfib is not None:
        #bad = np.isin(fs['FIBER'],badfib)
        #printlog('number at bad fibers '+str(sum(bad)),logger)
        #wfqa &= ~bad
        cat_out = mask_bad_fibers_time_dependent(fs[wfqa], badfib)
    else:
        cat_out = fs[wfqa]
    if remove_badfiber_spike_nz:
        badfib1 = np.loadtxt('bad_nz_fibers_ks_test.txt')
        badfib2 = np.loadtxt('elg_bad_nz_spike_fibers_1.498_1.499.txt')
        bad = np.isin(fs['FIBER'],np.concatenate((badfib1,badfib2)))
        cat_out = cat_out[~bad]
        return cat_out
    else:
        return cat_out
    
def mask_bad_petal_nights(dz, prog='dark'):
	if prog == 'dark':
		bad_petal_night_file = open('/global/cfs/cdirs/desi/survey/catalogs/DA2/LSS/loa-v1/lrg_bad_per_petal-night.txt','r')
	elif prog == 'bright':
		bad_petal_night_file = open('/global/cfs/cdirs/desi/survey/catalogs/DA2/LSS/loa-v1/bgs_bright_bad_per_petal-night.txt','r')

	inds_to_remove = np.array([])
	for line in bad_petal_night_file:
		night = int(line.split()[0])
		for petal in line.split()[1:]:
			inds_to_remove = np.concatenate(
			(inds_to_remove, np.where( 
			(dz['LASTNIGHT'] == night) 
			& (dz['FIBER'] >= 500 * int(petal))
			 & (dz['FIBER'] < 500 * (int(petal)+1)))[-1]))
	
	dz_inds = np.ones(len(dz)).astype('int')
	dz_inds[inds_to_remove.astype('int')] = 0
	return dz[dz_inds == 1]


def goodz_infull(tp,dz,zcol='Z_not4clus'):
    if tp == 'LRG':
        z_suc= dz['ZWARN']==0
        z_suc &= dz['DELTACHI2']>15
        z_suc &= dz[zcol]<1.5

    if tp == 'ELG':
        z_suc = dz['o2c'] > 0.9

    if tp == 'QSO':
        z_suc = dz[zcol]*0 == 0
        z_suc &= dz[zcol] != 999999
        z_suc &= dz[zcol] != 1.e20


    if tp == 'BGS':    
        z_suc = dz['ZWARN']==0
        z_suc &= dz['DELTACHI2']>40
    
    return z_suc

def make_hp(value, hpix, nside, fill_with=np.nan):
    """ A Function to create a HEALPix map
    """
    m_ = np.zeros(12*nside*nside)
    m_[:] = fill_with
    m_[hpix] = value

    return m_

def get_debv(mapname = '/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars_y3/v0.1/final_maps/lss/desi_ebv_lss_256.fits'):
    #DR1 map is in '/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/kp3_maps/'
    #DR1 map named like /global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/kp3_maps/v1_desi_ebv_256.fits
    #DR2 map is in /global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars_y3/v0.1/final_maps/lss/
    #DR2 map named like desi_ebv_lss_256.fits
    import healpy as hp

    #dirmap = '/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/v0/kp3_maps/'
    
    nside = 256#64
    nest = False
    eclrs = ['gr','rz']
    debv = Table()
    for ec in eclrs:
        #ebvn = fitsio.read(dirmap+'v0_desi_ebv_'+ec+'_'+str(nside)+'.fits')
        #ebvn = fitsio.read(dirmap+'v1_desi_ebv_'+str(nside)+'.fits')
        ebvn = fitsio.read(mapname)
        debv_a = ebvn['EBV_DESI_'+ec.upper()]-ebvn['EBV_SFD_'+ec.upper()]
        debv_a = hp.reorder(debv_a,r2n=True)
        debv['EBV_DIFF_'+ec.upper()] = debv_a
    
    return debv

def get_skyres():
    import healpy as hp
    sky_g = {'N':np.zeros(256*256*12),'S':np.zeros(256*256*12)}
    sky_r = {'N':np.zeros(256*256*12),'S':np.zeros(256*256*12)}
    sky_z = {'N':np.zeros(256*256*12),'S':np.zeros(256*256*12)}
    f = fitsio.read('/dvs_ro/cfs/cdirs/desi/users/rongpu/imaging_mc/ism_mask/sky_resid_map_256_north.fits')
    pixr = f['HPXPIXEL']
    pix_nest = hp.ring2nest(256,pixr)
    for i in range(0,len(f)):
        pix = pix_nest[i]#f['HPXPIXEL'][i]
        sky_g['N'][pix] = f['sky_median_g'][i]
        sky_r['N'][pix] = f['sky_median_r'][i]
        sky_z['N'][pix] = f['sky_median_z'][i]
    f = fitsio.read('/dvs_ro/cfs/cdirs/desi/users/rongpu/imaging_mc/ism_mask/sky_resid_map_256_south.fits')
    pix = f['HPXPIXEL']
    pix_nest = hp.ring2nest(256,pix)
    for i in range(0,len(f)):
        pix = pix_nest[i]#f['HPXPIXEL'][i]
        sky_g['S'][pix] = f['sky_median_g'][i]
        sky_r['S'][pix] = f['sky_median_r'][i]
        sky_z['S'][pix] = f['sky_median_z'][i]
    return sky_g,sky_r,sky_z

def cutphotmask(aa,bits=None,logger=None):
    printlog(str(len(aa)) +' before imaging veto' ,logger=logger)
    keep = (aa['NOBS_G']>0) & (aa['NOBS_R']>0) & (aa['NOBS_Z']>0)
    if bits is not None:
        for biti in bits:
            keep &= ((aa['MASKBITS'] & 2**biti)==0)
    aa = aa[keep]
    printlog(str(len(aa)) +' after imaging veto',logger=logger )
    return aa

def get_zcmbdipole(ra,dec):
    '''
    given input arrays for ra/dec in degrees, calculate redshift from CMB dipole
    '''
    """
    In https://arxiv.org/pdf/1807.06205
    Planck 2018 results. I. Overview and the cosmological legacy of Planck
    
    Table 3
    Sun–CMB :
    v = 369.82 +- 0.11 km/s
    long = 264.021 +- 0.011 deg
    lat  = 48.253 +- 0.005 deg
    
    Conversion to RA Dec with python:
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    dipole = SkyCoord(l=264.021*u.degree, b=48.253*u.degree, frame='galactic')
    print(dipole.fk5)
    167.94191028, -6.9442636
    """

    z_cmb=369.82/2.9979e5
    ra_cmb=167.942
    dec_cmb=-6.944

    sin_ra_cmb = np.sin(np.pi*ra_cmb/180)
    cos_ra_cmb = np.cos(np.pi*ra_cmb/180)
    sin_dec_cmb = np.sin(np.pi*dec_cmb/180)
    cos_dec_cmb = np.cos(np.pi*dec_cmb/180)

    sin_ra = np.sin(np.pi*ra/180)
    cos_ra = np.cos(np.pi*ra/180)
    sin_dec = np.sin(np.pi*dec/180)
    cos_dec = np.cos(np.pi*dec/180)
    
    cos_angdis2cmb = cos_dec*cos_dec_cmb*(cos_ra*cos_ra_cmb+sin_ra*sin_ra_cmb)+sin_dec*sin_dec_cmb
    
    return z_cmb*cos_angdis2cmb

def splitGC(input_array):
    import LSS.common_tools as common
    '''
    given and input array with RA, DEC
    return selection for NGC
    '''
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    c = SkyCoord(input_array['RA']* u.deg,input_array['DEC']* u.deg,frame='icrs')
    gc = c.transform_to('galactic')
    sel_ngc = gc.b > 0
    return sel_ngc

def select_regressis_DES(input_array,ra_col='RA',dec_col='DEC'):    
    '''
    input_array with RA, DEC given by ra_col,dec_col
    return selection for DES as defined by regressis
    '''

    from regressis import footprint
    import healpy as hp
    foot = footprint.DR9Footprint(256, mask_lmc=False, clear_south=True, mask_around_des=False, cut_desi=False)
    north, south, des = foot.get_imaging_surveys()
    th,phi = (-input_array[dec_col]+90.)*np.pi/180.,input_array[ra_col]*np.pi/180.
    pix = hp.ang2pix(256,th,phi,nest=True)
    sel_des = des[pix]
    return sel_des

def find_znotposs_tloc(dz,priority_thresh=10000):
    #dz should contain the potential targets of a given type, after cutting bad fibers
    #priority_thresh cuts repeat observations (so, e.g., 3000 work for dark time main survey)
    tileids = np.unique(dz['TILEID'])
    ual = []
    ufl = []
    for tile in tileids:
        sel = dz['TILEID'] == tile
        dzs = dz[sel]
        sela = dzs['ZWARN'] != 999999
        sela &= dzs['ZWARN']*0 == 0
        tlida = np.unique(dzs[sela]['TILELOCID']) #tilelocids with good assignments
        #print(tile,len(sela),len(dzs),np.sum(sela))
        tida = np.unique(dzs[sela]['TARGETID']) #targetids with good assignments
        ua = ~np.isin(dzs['TARGETID'],tida) #columns that were not assigned
        #ua &= dzs['NUMOBS'] == 0
        ua &= dzs['PRIORITY'] > priority_thresh #columns corresponding to targets with priority indicating observations unfinished
        ua &= ~np.isin(dzs['TILELOCID'],tlida) #columns corresponding to tilelocid that were not assigned
        #combination then gives the tilelocid of unassigned targets that have not finished observation; these must have been blocked by something
        uatlids = np.unique(dzs[ua]['TILELOCID'])
        ual.append(uatlids)
        selp = dzs['PRIORITY'] > priority_thresh
        tlids_gp = np.unique(dzs[selp]['TILELOCID'])
        tlids_all = np.unique(dzs['TILELOCID'])
        tlids_full = tlids_all[~np.isin(tlids_all,tlids_gp)]
        print('done with tile '+str(tile),str(len(uatlids)),str(len(tlids_full)))
        ufl.append(tlids_full)
    print('concatenating')
    ualt = np.concatenate(ual)
    uflt = np.concatenate(ufl)
    return ualt,uflt


def find_znotposs(dz,logname=None):

    dz.sort('TARGETID')
    tidnoz = []
    tids = np.unique(dz['TARGETID'])
    ti = 0
    i = 0
    message = 'finding targetids that were not observed'
    if logname is None:
        print(message)
    else:
        logger = logging.getLogger(logname)
        logger.info(message)
    while i < len(dz):
        za = 0

        while dz[i]['TARGETID'] == tids[ti]:
            if dz[i]['ZWARN'] != 999999:
                za = 1
                #break
            i += 1
            if i == len(dz):
                break
        if za == 0:
            tidnoz.append(tids[ti])

        if ti%30000 == 0:
            if logname is None:
                print(ti)
            else:
                logger.info('at index '+str(ti))
                
        ti += 1


    selnoz = np.isin(dz['TARGETID'],tidnoz)
    tidsb = np.unique(dz[selnoz]['TILELOCID'])
    #dz = dz[selnoz]
    dz.sort('TILELOCID')
    tids = np.unique(dz['TILELOCID'])
    message = 'number of targetids with no obs '+str(len(tidnoz))
    if logname is None:
        print(message)
    else:
        logger.info(message)
    tlidnoz = []
    lznposs = []

    ti = 0
    i = 0

    while i < len(dz):
        za = 0

        while dz[i]['TILELOCID'] == tids[ti]:
            if dz[i]['ZWARN'] != 999999:
                za = 1
                #break
            i += 1
            if i == len(dz):
                break
        if za == 0:
            tlidnoz.append(tids[ti])
            #if np.isin(tids[ti],tidsb):
            #    lznposs.append(tids[ti])

        if ti%30000 == 0:
            if logname is None:
                print(ti,len(tids))
            else:
                logger.info(str(ti)+' '+str(len(tids)))
        ti += 1
    #the ones to veto are now the join of the two
    wtbtlid = np.isin(tlidnoz,tidsb)
    tlidnoz = np.array(tlidnoz)
    lznposs = tlidnoz[wtbtlid]
    message = 'number of locations where assignment was not possible because of priorities '+str(len(lznposs))
    if logname is None:
        print(message)
    else:
        logger.info(message)
    return lznposs

def comp_tile(dz):
    compa = []
    tll = []
    ti = 0
    print('getting completenes')
    #sorting by tiles makes things quicker with while statements below
    dz.sort('TILES')
    nts = len(np.unique(dz['TILES']))
    tlsl = dz['TILES']
    tlslu = np.unique(tlsl)
    laa = dz['LOCATION_ASSIGNED']

    i = 0
    while i < len(dz):
        tls  = []
        tlis = []
        nli = 0
        nai = 0

        while tlsl[i] == tlslu[ti]:
            nli += 1 #counting unique targetids within the given TILES value
            nai += laa[i] #counting the number assigned
            i += 1
            if i == len(dz):
                break

        if ti%1000 == 0:
            print('at tiles '+str(ti)+' of '+str(nts))
        cp = nai/nli #completeness is number assigned over number total
        compa.append(cp)
        tll.append(tlslu[ti])
        ti += 1
    return tll,compa

def comp_tileloc(dz):

    locl,nlocl = np.unique(dz['TILELOCID'],return_counts=True)
    wz = dz['LOCATION_ASSIGNED'] == 1
    dzz = dz[wz]

    loclz,nloclz = np.unique(dzz['TILELOCID'],return_counts=True)

    print('getting fraction assigned for each tilelocid')
    #should be one (sometimes zero, though) assigned target at each tilelocid and we are now counting how many targets there are per tilelocid
    #probability of assignment is then estimated as 1/n_tilelocid
    nm = 0
    nmt =0

    loco = []
    fzo = []
    nlist = 0
    nlistg1 = 0
    for i in range(0,len(locl)):
        if i%10000 == 0:
            print('at row '+str(i))
        nt = nlocl[i]
        loc = locl[i]
        w = loclz == loc
        nz = 0
        if len(loclz[w]) == 1:
            nz = nloclz[w][0] #these are supposed all be 1...
        else:
            nm += 1.
            nmt += nt
        if len(loclz[w]) > 1:
            print('why is len(loclz[w]) > 1?') #this should never happen
        loco.append(loc)
        frac = nz/nt
        #if type(frac) != float:
        #    if len(frac) > 1:
        #        nlistg1 += 1
        #    frac = frac[0]
        #    nlist += 1

        fzo.append(frac)
    print(str(nlist)+ ' were type list for some reason; '+str(nlistg1)+ ' had length greater than 1')
    print('number of fibers with no observation, number targets on those fibers')
    print(nm,nmt)


    return loco,fzo


def mknz(fcd,fcr,fout,bs=0.01,zmin=0.01,zmax=1.6,randens=2500.,compmd='ran',wtmd='clus'):
    '''
    fcd is the full path to the catalog file in fits format with the data; requires columns Z and WEIGHT
    fcr is the full path to the random catalog meant to occupy the same area as the data; assumed to come from the imaging randoms that have a density of 2500/deg2
    fout is the full path to the file name
    bs is the bin width for the n(z) calculation
    zmin is the lower edge of the first bin
    zmax is the upper edge of the last bin
    '''
    #cd = distance(om,1-om)
    ranf = fitsio.read_header(fcr,ext=1) #should have originally had 2500/deg2 density, so can convert to area
    area = ranf['NAXIS2']/randens
    print('area is '+str(area))
    outf = open(fout,'w')
    outf.write('#area is '+str(area)+'square degrees\n')
    
    if compmd == 'ran':
        ranf = fitsio.read(fcr)
        area = np.sum(ranf['FRAC_TLOBS_TILES'])/randens
        outf.write('#effective area is '+str(area)+'square degrees\n')

    df = fitsio.read(fcd)

    nbin = int((zmax-zmin)*(1+bs/10)/bs)
    if wtmd == 'clus':
        #this is what should be used for clustering catalogs because 'WEIGHT' gets renormalized
        wts = df['WEIGHT_COMP']*df['WEIGHT_SYS']*df['WEIGHT_ZFAIL']
    zhist = np.histogram(df['Z'],bins=nbin,range=(zmin,zmax),weights=wts)
    outf.write('#zmid zlow zhigh n(z) Nbin Vol_bin\n')
    for i in range(0,nbin):
        zl = zhist[1][i]
        zh = zhist[1][i+1]
        zm = (zh+zl)/2.
        voli = area/(360.*360./np.pi)*4.*np.pi/3.*(dis_dc(zh)**3.-dis_dc(zl)**3.)
        nbarz =  zhist[0][i]/voli
        outf.write(str(zm)+' '+str(zl)+' '+str(zh)+' '+str(nbarz)+' '+str(zhist[0][i])+' '+str(voli)+'\n')
    outf.close()
    return True

def mknz_full(fcd,fcr,tp,bs=0.01,zmin=0.01,zmax=1.6,randens=2500.,write='n',md='data',zcol='Z_not4clus',reg=None):
    '''
    fcd is the full path to the catalog file in fits format with the data; requires columns Z and WEIGHT
    fcr is the full path to the random catalog meant to occupy the same area as the data; assumed to come from the imaging randoms that have a density of 2500/deg2
    bs is the bin width for the n(z) calculation
    zmin is the lower edge of the first bin
    zmax is the upper edge of the last bin
    returns array with n(z) values
    '''
    #cd = distance(om,1-om)
    if reg is None:
        ranf = fitsio.read_header(fcr,ext=1) #should have originally had 2500/deg2 density, so can convert to area
        area = ranf['NAXIS2']/randens
    else:
        print(reg)
        ranf = fitsio.read(fcr,columns=["PHOTSYS"])
        selreg = ranf['PHOTSYS'] == reg
        area = len(ranf[selreg])/randens
        del ranf
        
    print('area is '+str(area))

    df = fitsio.read(fcd)
    if md == 'data':
        gz = goodz_infull(tp,df)
    if md == 'mock':
        gz = df['ZWARN'] == 0
    df = df[gz]
    wo = ''
    if reg is not None:
        selreg = df['PHOTSYS'] == reg
        df = df[selreg]
        wo = '_'+reg
    nbin = int((zmax-zmin)/bs)
    cols = list(df.dtype.names)
    wts = 1/df['FRACZ_TILELOCID']
    if 'FRAC_TLOBS_TILES' in cols:
        wts *= 1/df['FRAC_TLOBS_TILES']
    #if 'WEIGHT_SYS' in cols:
    #    #wtnorm = np.mean()
    #    wts *= df['WEIGHT_SYS']
    selnan = wts*0 != 0
    print('number of nans in weights '+str(np.sum(selnan)))
    wts[selnan] = 1.
    zhist = np.histogram(df[zcol],bins=nbin,range=(zmin,zmax),weights=wts)
    zl = zhist[1][:-1]
    zh = zhist[1][1:]
    zm = (zl+zh)/2.
    vol = area/(360.*360./np.pi)*4.*np.pi/3.*(dis_dc(zh)**3.-dis_dc(zl)**3.)
    nz = zhist[0]/vol
    #print(nz)
    if write == 'y':
        fout = fcd.replace('.dat.fits','')+wo+'_nz.txt'
        fout = fout.replace('dvs_ro','global')
        outf = open(fout,'w')
        outf.write('#area is '+str(area)+'square degrees\n')
        outf.write('#zmid zlow zhigh n(z) Nbin Vol_bin\n')

        for i in range(0,len(nz)):
            outf.write(str(zm[i])+' '+str(zl[i])+' '+str(zh[i])+' '+str(nz[i])+' '+str(zhist[0][i])+' '+str(vol[i])+'\n')
        outf.close()
    return nz

def get_cols_rzdshuf(tabd,cols_rd,cols_z):
    #shuffle ra,dec and z from data to make randoms
    inds_rd = np.random.choice(len(tabd),len(tabd))
    inds_z = np.random.choice(len(tabd),len(tabd))
    dshuf_rd = tabd[inds_rd]
    dshuf_z = tabd[inds_z]
    tabr = Table()
    for col in cols_rd:
        tabr[col] =  dshuf_rd[col]
    for col in cols_z:
        tabr[col] =  dshuf_z[col]
    return tabr

def clusran_shufrd(flin,ran_sw='',P0=10000,zmin=0.01,zmax=1.6,dz=0.01):
    #take existing data/random clustering catalogs and re-sample redshift dependent quantities to assign to randoms
    import LSS.common_tools as common
    fcdn = Table.read(flin+'_clustering.dat.fits')
    cols_rd = ['RA','DEC','NTILE']
    cols_z = ['Z']


    def _resamp(selregd,fcdn):
        tabsr = []
        fcdnn = fcdn[selregd]
        fcdns = fcdn[~selregd]
        tabsd = [fcdnn,fcdns]
        rdl =[]
        for i in range(0,len(tabsd)):
            tabr = get_cols_rzdshuf(tabsd[i],cols_rd,cols_z)
            tabsr.append(tabr)
        ffr = vstack(tabsr)   
        return ffr

    if 'NGC' in flin:
        selregd = fcdn['DEC'] > 32.375
        ffr = _resamp(selregd,fcdn)
    elif 'SGC' in flin and 'QSO' in flin:
        print('resampling in DES region')
        from regressis import footprint
        foot = footprint.DR9Footprint(256, mask_lmc=False, clear_south=True, mask_around_des=False, cut_desi=False)
        north, south, des = foot.get_imaging_surveys()
        th_ran,phi_ran = (-ffr['DEC']+90.)*np.pi/180.,ffr['RA']*np.pi/180.
        th_dat,phi_dat = (-fcdn['DEC']+90.)*np.pi/180.,fcdn['RA']*np.pi/180.
        pixr = hp.ang2pix(256,th_ran,phi_ran,nest=True)
        selregr = des[pixr]
        pixd = hp.ang2pix(256,th_dat,phi_dat,nest=True)
        selregd = des[pixd]
        ffr = _resamp(selregr,selregd,ffr,fcdn)
    else:
        ffr = get_cols_rzdshuf(fcdn,cols_rd,cols_z)

    comp_ntl = get_comp(flin)

    nzd = np.loadtxt(flin+'_nz.txt').transpose()[3] #column with nbar values
    zl = ffr['Z']
    nl = np.zeros(len(zl))
    for ii in range(0,len(zl)):
        z = zl[ii]
        zind = int((z-zmin)/dz)
        if z > zmin and z < zmax:
            nl[ii] = nzd[zind]
    ffr['NX'] = nl*comp_ntl[ffr['NTILE']-1]

    fkpl = 1/(1+ffr['NX']*P0)
    ffr['WEIGHT_FKP'] = fkpl
    ffr['WEIGHT'] = np.ones(len(ffr))
    return ffr
    
    
    #outfn =  flin+'_'+str(rann)+'_rdshuf_clustering.ran.fits'
    
    #if write_cat == 'y':
        #comments = ["'clustering' LSS catalog for random number "+str(rann)+", BASS/MzLS region","entries are only for data with good redshifts"]
    #    common.write_LSS(ffr,outfn)

def get_comp(fb,ran_sw=''):
    fn = fb.replace(ran_sw,'')+'_clustering.dat.fits'
    fd = Table(fitsio.read(fn))
    mean_comp = len(fd)/np.sum(fd['WEIGHT_COMP'])
    print('mean completeness '+str(mean_comp))
    ntl = np.unique(fd['NTILE'])
    comp_ntl = np.zeros(len(ntl))
    weight_ntl = np.zeros(len(ntl))
    for i in range(0,len(ntl)):
        sel = fd['NTILE'] == ntl[i]
        mean_ntweight = np.mean(fd['WEIGHT_COMP'][sel])        
        weight_ntl[i] = mean_ntweight
        comp_ntl[i] = 1/mean_ntweight#*mean_fracobs_tiles
    fran = fitsio.read(fb+'_0_clustering.ran.fits',columns=['NTILE','FRAC_TLOBS_TILES'])
    fttl = np.zeros(len(ntl))
    for i in range(0,len(ntl)): 
        sel = fran['NTILE'] == ntl[i]
        mean_fracobs_tiles = np.mean(fran[sel]['FRAC_TLOBS_TILES'])
        fttl[i] = mean_fracobs_tiles
    print(comp_ntl,fttl)
    comp_ntl = comp_ntl*fttl
    print('completeness per ntile:')
    print(comp_ntl)
    return comp_ntl

def get_weight_ntile(in_data):
    '''
    in_data should be an array of the clustering catalog data
    returns the mean weight as a function of NTILE, dividing by the completeness weight
    this is for eventual use in angular upweighting
    '''
    ntl = np.unique(in_data['NTILE'])
    weight_ntl = np.ones(len(ntl))
    fkp_ntl = np.ones(len(ntl))
    for i in range(0,len(ntl)):
        sel = in_data['NTILE'] == ntl[i]
        mean_ntweight = np.mean(in_data['WEIGHT'][sel]/in_data['WEIGHT_COMP'][sel]) #we want the mean weight used without WEIGHT_COMP       
        weight_ntl[i] = mean_ntweight
        fkp_ntl[i] = np.mean(in_data['WEIGHT_FKP'][sel])
    return weight_ntl,fkp_ntl

def add_weight_ntile(fb,logger=None,ranmin=0,nran=18,par='n',extradir='',tp='',nproc=9):
    from desitarget.internal import sharedmem
    fn = fb+'_NGC_clustering.dat.fits'
    clus_dn = fitsio.read(fn.replace('global','dvs_ro').replace(tp,extradir+tp))
    fs = fb+'_SGC_clustering.dat.fits'
    clus_ds = fitsio.read(fs.replace('global','dvs_ro').replace(tp,extradir+tp))
    clus_d = np.concatenate((clus_dn,clus_ds))
    weight_ntl,fkp_ntl = get_weight_ntile(clus_d)
    
    full_name = fb+'_full_HPmapcut.dat.fits'
    full_d = Table(fitsio.read(full_name))
    full_d['WEIGHT_NTILE'] = weight_ntl[full_d['NTILE']-1]
    full_d['WEIGHT_FKP_NTILE'] = fkp_ntl[full_d['NTILE']-1]
    write_LSS_scratchcp(full_d,full_name,logger=logger)
    
    def _parfun(rann):
        fn = fb+'_NGC_'+str(rann)+'_clustering.ran.fits'
        clus_rn = fitsio.read(fn.replace('global','dvs_ro').replace(tp,extradir+tp) )
        fs = fb+'_SGC_'+str(rann)+'_clustering.ran.fits'
        clus_rs = fitsio.read(fs.replace('global','dvs_ro').replace(tp,extradir+tp) )
        clus_r = np.concatenate((clus_rn,clus_rs))
        #weight_ntl,fkp_ntl = get_weight_ntile(clus_r) #just base it on data so it is equivalent for angular upweighting
    
        full_name = fb+'_'+str(rann)+'_full_HPmapcut.ran.fits'
        full_r = Table(fitsio.read(full_name))
        full_r['WEIGHT_NTILE'] = weight_ntl[full_r['NTILE']-1]
        full_r['WEIGHT_FKP_NTILE'] = fkp_ntl[full_r['NTILE']-1]
        write_LSS_scratchcp(full_r,full_name,logger=logger)

    if par == 'n':
        for rann in range(ranmin,nran):
            _parfun(rann)
            #print('done with random number '+str(rann))
    else:
        inds = np.arange(ranmin,nran)
        from multiprocessing import Pool
    
        #nproc = 9 #try this so doesn't run out of memory
        pool = sharedmem.MapReduce(np=nproc)
        #with Pool() as pool:#Pool(processes=nproc) as pool:
        with pool:
            res = pool.map(_parfun, inds)

    return True
    
    
    


def addnbar(fb,nran=18,bs=0.01,zmin=0.01,zmax=1.6,P0=10000,add_data=True,ran_sw='',ranmin=0,compmd='ran',par='n',nproc=18,logger=None):
    '''
    fb is the root of the file name, including the path
    nran is the number of random files to add the nz to
    bs is the bin size of the nz file (read this from file in future)
    zmin is the lower edge of the minimum bin (read this from file in future)
    zmax is the upper edge of the maximum bin (read this from file in the future)
    '''

    from desitarget.internal import sharedmem
    nzd = np.loadtxt(fb.replace(ran_sw,'')+'_nz.txt').transpose()[3] #column with nbar values
    fn = fb.replace(ran_sw,'')+'_clustering.dat.fits'
    #ff = fitsio.FITS(fn,'rw')
    #fd = Table(ff['LSS'].read())
    #fd = fitsio.read(fn) #reading in data with fitsio because it is much faster to loop through than table
    fd = Table(fitsio.read(fn.replace('global','dvs_ro')))
    zl = fd['Z']
    nl = np.zeros(len(zl))
    for ii in range(0,len(zl)):
        z = zl[ii]
        zind = int((z-zmin)/bs)
        if z > zmin and z < zmax:
            nl[ii] = nzd[zind]
    mean_comp = len(fd)/np.sum(fd['WEIGHT_COMP'])
    if logger is None:
        print('mean completeness '+str(mean_comp))
    else:
        logger.info('mean completeness '+str(mean_comp))
    nont = 0
    if 'NTILE' not in list(fd.dtype.names):
        fd['NTILE'] = np.ones(len(fd),dtype=int)
        if logger is None:
            print('added NTILE = 1 column because column did not exist')
        else:
            logger.info('added NTILE = 1 column because column did not exist')
        nont = 1
    ntl = np.unique(fd['NTILE'])
    comp_ntl = np.ones(len(ntl))
    weight_ntl = np.ones(len(ntl))
    for i in range(0,len(ntl)):
        sel = fd['NTILE'] == ntl[i]
        mean_ntweight = np.mean(fd['WEIGHT_COMP'][sel])        
        weight_ntl[i] = mean_ntweight
        comp_ntl[i] = 1/mean_ntweight#*mean_fracobs_tiles
    
    if compmd == 'ran':
        fran = fitsio.read(fb.replace('global','dvs_ro')+'_0_clustering.ran.fits',columns=['NTILE','FRAC_TLOBS_TILES'])
        fttl = np.zeros(len(ntl))
        for i in range(0,len(ntl)): 
            sel = fran['NTILE'] == ntl[i]
            mean_fracobs_tiles = np.mean(fran[sel]['FRAC_TLOBS_TILES'])
            fttl[i] = mean_fracobs_tiles
    else:
        fttl = np.ones(len(ntl))
    print(comp_ntl,fttl)
    comp_ntl = comp_ntl*fttl
    print('completeness per ntile:')
    print(comp_ntl)
    #del fd
    #ft = Table.read(fn)
    #ft['NZ'] = nl
    #fd['NZ'] = nl
    fd['NX'] = nl*comp_ntl[fd['NTILE']-1]
    #AJR is unsure why multiply these three here rather than use the original 'WEIGHT'...
    fd['WEIGHT'] = fd['WEIGHT_COMP']*fd['WEIGHT_SYS']*fd['WEIGHT_ZFAIL']/weight_ntl[fd['NTILE']-1]
    cols = list(fd.dtype.names)
    if 'WEIGHT_BLIND' in cols:
        fd['WEIGHT'] *= fd['WEIGHT_BLIND']
    #ff['LSS'].insert_column('NZ',nl)
    print(np.min(nl),np.max(nl))

    #fkpl = 1./(1+nl*P0*mean_comp)
    #fkpl = comp_ntl[fd['NTILE']-1]/(1+nl*P0*comp_ntl[fd['NTILE']-1])
    fkpl = 1/(1+fd['NX']*P0)
    #ft['WEIGHT_FKP'] = 1./(1+ft['NZ']*P0)
    if add_data:
        fd['WEIGHT_FKP'] = fkpl
        write_LSS_scratchcp(fd,fn,logger=logger)
    #fd = np.array(fd)
    #ff['LSS'].insert_column('WEIGHT_FKP',fkpl)
    #ff['LSS'].write(fd)
    #ff['LSS'].write_history("added NZ and WEIGHT_FKP columns on "+datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
    #ff.close()
    #ft.write(fn,format='fits',overwrite=True)
    print('done with data')
    def _parfun(rann):
        fn = fb+'_'+str(rann)+'_clustering.ran.fits'
        #ff = fitsio.FITS(fn,'rw')
        #fd = ff['LSS'].read()
        fd = Table(fitsio.read(fn.replace('global','dvs_ro') ))
        #fd = fitsio.read(fn) #reading in data with fitsio because it is much faster to loop through than table
        zl = fd['Z']
        nl = np.zeros(len(zl))
        for ii in range(0,len(zl)):
            z = zl[ii]
            zind = int((z-zmin)/bs)
            if z > zmin and z < zmax:
                nl[ii] = nzd[zind]
        #del fd
        #ft = Table.read(fn)
        #ft['NZ'] = nl
        #ff['LSS'].insert_column('NZ',nl)
        #fd['NZ'] = nl
        print('mean completeness '+str(mean_comp))
        if nont == 1:
            fd['NTILE'] = np.ones(len(fd),dtype=int)
            print('added NTILE = 1 column because column did not exist')

        fd['NX'] = nl*comp_ntl[fd['NTILE']-1]
        #the following lines should keep the relative weighting and the region normalization in place
        wt = fd['WEIGHT_COMP']*fd['WEIGHT_SYS']*fd['WEIGHT_ZFAIL']
        if compmd == 'ran':
            wt *= fd['FRAC_TLOBS_TILES']
        cols = list(fd.dtype.names)
        if 'WEIGHT_BLIND' in cols:
            wt *= fd['WEIGHT_BLIND']

        wtfac = np.ones(len(fd))
        sel = wt > 0
        wtfac[sel] = fd['WEIGHT'][sel]/wt[sel]
        print(np.mean(wtfac))
        fd['WEIGHT'] = wtfac*wt/weight_ntl[fd['NTILE']-1] #this should keep, e.g., N/S normalization in place
        #fkpl = 1./(1+nl*P0*mean_comp)
        #fkpl = comp_ntl[fd['NTILE']-1]/(1+nl*P0*comp_ntl[fd['NTILE']-1])
        fkpl = 1/(1+fd['NX']*P0)
        fd['WEIGHT_FKP'] = fkpl
        write_LSS_scratchcp(fd,fn,logger=logger)
        #ff['LSS'].insert_column('WEIGHT_FKP',fkpl)
        #fd = np.array(fd)
        #ff['LSS'].write(fd)
        #ff['LSS'].write_history("added NZ and WEIGHT_FKP columns on "+datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        #ff.close()
        #ft['WEIGHT_FKP'] = 1./(1+ft['NZ']*P0)
        #ft.write(fn,format='fits',overwrite=True)
    if par == 'n':
        for rann in range(ranmin,nran):
            _parfun(rann)
            #print('done with random number '+str(rann))
    else:
        inds = np.arange(ranmin,nran)
        from multiprocessing import Pool
    
        #nproc = 9 #try this so doesn't run out of memory
        pool = sharedmem.MapReduce(np=nproc)
        #with Pool(processes=nproc) as pool:
        with pool:
            res = pool.map(_parfun, inds)

    return True

def addFKPfull(fb,nz,tp,bs=0.01,zmin=0.01,zmax=1.6,P0=10000,add_data=True,md='data',zcol='Z_not4clus',logger=None):
    '''
    fb is the file name, including the path
    nran is the number of random files to add the nz to
    bs is the bin size of the nz file (read this from file in future)
    zmin is the lower edge of the minimum bin (read this from file in future)
    zmax is the upper edge of the maximum bin (read this from file in the future)
    '''

    fd = Table(fitsio.read(fb))
    
    zl = fd[zcol]
    zind = ((zl-zmin)/bs).astype(int)
    gz = fd['ZWARN'] != 999999
    if md == 'data':
        gz &= goodz_infull(tp,fd)
    gz &= zl > zmin
    gz &= zl < zmax
    
    print(np.min(fd[gz]['FRACZ_TILELOCID']),np.max(fd[gz]['FRACZ_TILELOCID']))
    nl = np.zeros(len(fd))
    nl[gz] = nz[zind[gz]]
    mean_comp = len(fd[gz])/np.sum(1./fd[gz]['FRACZ_TILELOCID'])
    printlog('mean completeness '+str(mean_comp),logger)

    fkpl = 1./(1+nl*P0*mean_comp)
    #ft['WEIGHT_FKP'] = 1./(1+ft['NZ']*P0)
    if add_data:
        fd['WEIGHT_FKP'] = fkpl
        write_LSS_scratchcp(fd,fb.replace('dvs_ro','global'),logger=logger)
    return True

def join_with_fastspec(infn,fscols=['TARGETID','ABSMAG01_SDSS_G','ABSMAG01_SDSS_R'],inroot='/dvs_ro/cfs/cdirs/',\
outroot='/global/cfs/cdirs/',fsver='v2.0',fsrel='dr1',specrel='iron',prog='bright'):
    #v2.0 used in v1 through v1.2 of catalogs
    indata = fitsio.read(inroot+infn)
    print(len(indata))
    #fsfn = '/dvs_ro/cfs/cdirs/desi/public/'+fsrel+'/vac/'+fsrel+'/fastspecfit/'+specrel+'/'+fsver+'/catalogs/fastspec-'+specrel+'-main-'+prog+'.fits'
    if fsver == 'latest':
        fsfn = '/dvs_ro/cfs/cdirs/desi/spectro/fastspecfit/'+specrel+'/catalogs/fastspec-'+specrel+'-main-'+prog+'.fits'
    else:
        fsfn = '/dvs_ro/cfs/cdirs/desi/public/'+fsrel+'/vac/'+fsrel+'/fastspecfit/'+specrel+'/'+fsver+'/catalogs/fastspec-'+specrel+'-main-'+prog+'.fits'
    fastspecdata = fitsio.read(fsfn,columns=fscols)
    dcols = list(indata.dtype.names)
    indata = Table(indata)
    
    for col in fscols:
        if col != 'TARGETID':
            if col in dcols:
                indata.remove_column(col)
                print('removed '+col+' from input data')
                
    jt = join(indata,Table(fastspecdata),keys=['TARGETID'],join_type='left')
    jt = Table(jt, masked=True, copy=False)
    print(len(jt),'length of joined table, should agree with above')
    for col in fscols:
        jt[col] = jt[col].filled(999999)
        sel = jt[col] == 999999
        print(col,len(jt[sel]),'number with 999999')
    outfn = outroot+infn
    write_LSS(jt,outfn)

def add_dered_flux(data,fcols=['G','R','Z','W1','W2']):
    #data should be table with fcols flux columns existing
    for col in fcols:
        data['flux_'+col.lower()+'_dered'] = data['FLUX_'+col]/data['MW_TRANSMISSION_'+col]
    return data

def add_ke(dat,zcol='Z'):#,n_processes=100):
    from multiprocessing import Pool
    #dat should be table with flux_g_dered and flux_r_dered
    #from kcorr package https://github.com/SgmAstro/DESI, needs to be added to path
    #
    #ke_code_root = '/global/homes/a/ajross/desicode/DESI_ke'
    ke_code_root = os.environ['LSSCODE']+'/LSS/py/LSS/DESI_ke'
    sys.path.append(ke_code_root)
    os.environ['CODE_ROOT'] = ke_code_root
    print(os.environ['CODE_ROOT'])
    from   smith_kcorr     import GAMA_KCorrection
    from   rest_gmr        import smith_rest_gmr,rest_gmr
    from   tmr_ecorr       import tmr_ecorr, tmr_q

    kcorr_r   = GAMA_KCorrection(band='R')
    kcorr_g   = GAMA_KCorrection(band='G')

    r_dered = 22.5 - 2.5*np.log10(dat['flux_r_dered'])
    g_dered = 22.5 - 2.5*np.log10(dat['flux_g_dered'])
    gmr = g_dered-r_dered

#     chunk_size = len(dat)//n_processes
#     list = []
#     for i in range(0,n_processes):
#         list.append(0)
#     def _wrapper(N):
#         mini = N*chunk_size
#         maxi = mini+chunk_size
#         if maxi > len(dat):
#             maxi = len(dat)
#         idx = np.arange(mini,maxi)
#         data = Table()
#         data['idx'] = idx
#         data['REST_GMR_0P1'], rest_gmr_0p1_warn = smith_rest_gmr(dat[zcol][mini:maxi], gmr[mini:maxi])
#         list[N] = data
#         #return data
#
#     with Pool(processes=n_processes+1) as pool:
#         #res = pool.map(_wrapper, np.arange(n_processes))
#         pool.map(_wrapper, np.arange(n_processes))
#
#     res = vstack(list)#vstack(res)
#     res.sort('idx')
#     res.remove_column('idx')
#     print(len(res),len(dat))
    selz = dat[zcol] > 0
    selz &= dat[zcol] < 4.5
    cols = ['REST_GMR_0P1','KCORR_R0P1','KCORR_G0P1','KCORR_R0P0','KCORR_G0P0','REST_GMR_0P0','EQ_ALL_0P0','EQ_ALL_0P1','ABSMAG_RP1','ABSMAG_RP0']
    for col in cols:
        dat[col] = np.zeros(len(dat))
    dat['REST_GMR_0P1'][selz], rest_gmr_0p1_warn = smith_rest_gmr(dat[zcol][selz], gmr[selz])
    dat['KCORR_R0P1'][selz] = kcorr_r.k(dat[zcol][selz], dat['REST_GMR_0P1'][selz])
    dat['KCORR_G0P1'][selz] = kcorr_g.k(dat[zcol][selz], dat['REST_GMR_0P1'][selz])
    dat['KCORR_R0P0'][selz] = kcorr_r.k_nonnative_zref(0.0, dat[zcol][selz], dat['REST_GMR_0P1'][selz])
    dat['KCORR_G0P0'][selz] = kcorr_g.k_nonnative_zref(0.0, dat[zcol][selz], dat['REST_GMR_0P1'][selz])
    dat['REST_GMR_0P0'][selz] = gmr[selz] - (dat['KCORR_G0P0'][selz] - dat['KCORR_R0P0'][selz])
    dat['EQ_ALL_0P0'][selz]   = tmr_ecorr(dat[zcol][selz], dat['REST_GMR_0P0'][selz], aall=True)
    dat['EQ_ALL_0P1'][selz]   = tmr_ecorr(dat[zcol][selz], dat['REST_GMR_0P1'][selz], aall=True)
    dat['ABSMAG_RP1'][selz] = r_dered[selz] -dm(dat[zcol][selz])-dat['KCORR_R0P1'][selz]-dat['EQ_ALL_0P1'][selz]
    dat['ABSMAG_RP0'][selz] = r_dered[selz] -dm(dat[zcol][selz])-dat['KCORR_R0P0'][selz]-dat['EQ_ALL_0P0'][selz]
    return dat
    #abg = g_dered -dm(data['Z'])


def join_etar(fn,tracer,tarver='1.1.1'):
    tr = tracer
    if tr == 'BGS_BRIGHT':
        tr = 'BGS_ANY'
    etar_fn = '/dvs_ro/cfs/cdirs/desi/survey/catalogs/main/LSS/'+tr+'targets_pixelDR9v'+tarver+'.fits'
    ef = Table(fitsio.read(etar_fn))
    ef.remove_columns(['BRICKNAME','RA','DEC','PHOTSYS','DESI_TARGET'])
    df = fitsio.read(fn)
    print(len(df))
    df = join(df,ef,keys=['TARGETID'])
    print(len(df),'should match above')
    comments = ['Adding imaging mask column']
    write_LSS(df,fn,comments)


def add_map_cols(fn,rann,logger=None,new_cols=['HALPHA', 'HALPHA_ERROR', 'CALIB_G', 'CALIB_R', 'CALIB_Z', 'EBV_MPF_Mean_FW15', 'EBV_MPF_Mean_ZptCorr_FW15', 'EBV_MPF_Var_FW15', 'EBV_MPF_VarCorr_FW15', 'EBV_MPF_Mean_FW6P1', 'EBV_MPF_Mean_ZptCorr_FW6P1', 'EBV_MPF_Var_FW6P1', 'EBV_MPF_VarCorr_FW6P1', 'EBV_SGF14', 'BETA_ML', 'BETA_MEAN', 'BETA_RMS', 'HI', 'KAPPA_PLANCK'],fid_cols=['EBV','PSFDEPTH_G','PSFDEPTH_R','PSFDEPTH_Z','GALDEPTH_G','GALDEPTH_R','GALDEPTH_Z','PSFDEPTH_W1','PSFDEPTH_W2','PSFSIZE_G','PSFSIZE_R','PSFSIZE_Z'],redo=True):
    fid_fn = '/dvs_ro/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/randoms-1-'+str(rann)+'.fits'
    new_fn = '/dvs_ro/cfs/cdirs/desi/survey/catalogs/external_input_maps/mapvalues/randoms-1-'+str(rann)+'-skymapvalues.fits'
    mask_fn = '/dvs_ro/cfs/cdirs/desi/survey/catalogs/external_input_maps/maskvalues/randoms-1-'+str(rann)+'-skymapmask.fits'
   
    printlog('reading '+fn.replace('global','dvs_ro'),logger)
    df = Table(fitsio.read(fn.replace('global','dvs_ro')))
    printlog('read '+fn.replace('global','dvs_ro'),logger)
    
    col = 'SKYMAP_MASK'
    domask = True
    if np.isin(col,list(df.dtype.names)):
        printlog(col+' already in '+fn,logger)
        if redo:
            df.remove_columns([col])
            printlog('will replace '+col,logger) 
        else:
            printlog('not replacing '+col,logger)
            domask = False
    if domask:
        mask = fitsio.read(mask_fn)
        df = join(df,mask,keys=['TARGETID'])
        del mask
        
    cols2read_new = ['TARGETID']
    for col in new_cols:
        if np.isin(col,list(df.dtype.names)):
            printlog(col+' already in '+fn,logger)
            if redo:
                df.remove_columns([col])
                printlog('will replace '+col,logger)
                cols2read_new.append(col)
            else:
                printlog('not replacing '+col,logger)
        else:
            cols2read_new.append(col)

    rannew = fitsio.read(new_fn,columns=cols2read_new)
    printlog('read '+new_fn,logger)
    df = join(df,rannew,keys=['TARGETID'])
    printlog('joined '+new_fn,logger)
    del rannew
    cols2read_fid = ['TARGETID']
    for col in fid_cols:
        if np.isin(col,list(df.dtype.names)):
            printlog(col+' already in '+fn,logger)
            if redo:
                df.remove_columns([col])
                printlog('will replace '+col,logger)
                cols2read_fid.append(col)
            else:
                printlog('not replacing '+col,logger)
        else:
            cols2read_fid.append(col)

    ranfid = fitsio.read(fid_fn,columns=cols2read_fid)
    printlog('read '+fid_fn,logger)
    df = join(df,ranfid,keys=['TARGETID'])

    #print(len(df))
    #comments = ['Adding map columns']
    printlog('about to write to '+fn,logger)
    write_LSS_scratchcp(df,fn,logger=logger)#,comments)
    printlog('wrote to '+fn,logger)
    return

def add_veto_col(fn,ran=False,tracer_mask='lrg',rann=0,tarver='targetsDR9v1.1.1',redo=False,logger=None,return_array=False):
    mask_fn = '/dvs_ro/cfs/cdirs/desi/survey/catalogs/main/LSS/'+tracer_mask.upper()+tarver+'_'+tracer_mask+'imask.fits'
    if ran:
        mask_fn = '/dvs_ro/cfs/cdirs/desi/survey/catalogs/main/LSS/randoms-1-'+str(rann)+tracer_mask+'imask.fits'
    maskf = fitsio.read(mask_fn)
    df = fitsio.read(fn.replace('global','dvs_ro'))
    if np.isin(tracer_mask+'_mask',list(df.dtype.names)):
        printlog('mask column already in '+fn,logger)
        if redo:
            df = Table(df)
            df.remove_columns([tracer_mask+'_mask'])
            print('will replace '+tracer_mask)
        elif return_array:
            return df
        else:    
            return True
    else:
        printlog('adding '+tracer_mask,logger)
    #print(len(df))
    sel = np.isin(maskf['TARGETID'],df['TARGETID'])
    maskf = maskf[sel]
    #print(len(maskf))
    maskf = Table(maskf)
    maskf.sort('TARGETID')
    df = Table(df)
    df.sort('TARGETID')
    if np.array_equal(df['TARGETID'],maskf['TARGETID']):
        df[tracer_mask+'_mask'] = maskf[tracer_mask+'_mask']
    else:
        return('TARGETIDs do not match! exiting')
    #df = join(df,maskf,keys=['TARGETID'])
    #print(len(df),'should match above')
    #comments = ['Adding imaging mask column']
    if return_array:
        return df
    write_LSS_scratchcp(df,fn,logger=logger)#,comments)
    del df

def parse_circandrec_mask(custom_mask_fn):
    '''
    Parse the custom mask file and return astropy tables.
​
    Args:
        custom_mask_fn: full path to the custom mask file
    Returns:
        circ_mask and rect_mask: circular and rectangular masks in astropy table format
    '''
    with open(custom_mask_fn, 'r') as f:
        lines = list(map(str.strip, f.readlines()))
    circ_mask_arr = []
    rect_mask_arr = []
    for line in lines:
        if line!='' and line[0]!='#':
            line = line[:line.find('#')]
            line = list(map(float, line.split(',')))
            if len(line)==3:
                circ_mask_arr.append(line)
            elif len(line)==4:
                rect_mask_arr.append(line)
            else:
                raise ValueError
    
    circ_mask_arr = np.array(circ_mask_arr)
    rect_mask_arr = np.array(rect_mask_arr)
    
    circ_mask = Table(circ_mask_arr, names=['RA', 'DEC', 'radius'])
    rect_mask = Table(rect_mask_arr, names=['ramin', 'ramax', 'decmin', 'decmax'])
    
    return circ_mask, rect_mask

def maskcircandrec(indata,maskfn,logger=None):
    '''
    indata should have RA,DEC columns
    mask should be path to text file with entries for circular and rectangular masks
    outputs indices to mask
    '''
    circ_mask,rect_mask = parse_circandrec_mask(maskfn)
    mask_rect = np.zeros(len(indata),dtype='bool')
    for radec in rect_mask:
        ramin, ramax, decmin, decmax = radec
        mask_rect |= (indata['RA']>ramin) & (indata['RA']<ramax) & (indata['DEC']>decmin) & (indata['DEC']<decmax)
    printlog('comparison of data removed by rectangular mask: '+str(len(indata))+','+str(len(indata[mask_rect])),logger)
    
    dat_cx = np.cos(np.radians(indata['RA']))*np.cos(np.radians(indata['DEC']))
    dat_cy = np.sin(np.radians(indata['RA']))*np.cos(np.radians(indata['DEC']))
    dat_cz = np.sin(np.radians(indata['DEC']))

    mask_cx = np.cos(np.radians(circ_mask['RA']))*np.cos(np.radians(circ_mask['DEC']))
    mask_cy = np.sin(np.radians(circ_mask['RA']))*np.cos(np.radians(circ_mask['DEC']))
    mask_cz = np.sin(np.radians(circ_mask['DEC']))
    
    mat1 = np.array([dat_cx, dat_cy, dat_cz]).T
    mat2 = np.array([mask_cx, mask_cy, mask_cz])
    
    dist = np.dot(mat1, mat2)
    
    mask = np.any(dist>np.cos(np.radians(circ_mask['radius']/3600.)), axis=1)
    
    printlog('comparison of data removed by circular mask: '+str(len(indata))+','+str(len(indata[mask])),logger)
    #print(len(indata),len(indata[mask]))
    
    mask |= mask_rect
    
    printlog('comparison of data removed by both rectangular and circular mask: '+str(len(indata))+','+str(len(indata[mask])),logger)
    #print(len(indata),len(indata[mask]))
    
    return mask
    
    
    

def apply_veto(fin,fout=None,ebits=None,zmask=False,maxp=3400,comp_only=False,reccircmasks=None,wo='y',mapveto='',logger=None):
    '''
    fl is a string with the path to the file name to load
    fout is a string with the path to the outpur file
    ebits are the new imaging mask bits to apply
    zmask is whether or not to apply any zmask
    maxp is the maximum priority to keep in the data files
    '''
    if isinstance(fin, str):
        ff = Table(fitsio.read(fin.replace('global','dvs_ro')))#+'full_noveto.'+dr+'.fits')
    else:
        ff = fin
        del fin
    printlog('length of input '+str(len(ff)),logger)
    seld = ff['GOODHARDLOC'] == 1
    printlog('length after cutting to good locations '+str(len(ff[seld])),logger)
    if '.dat' in fout:
        #seld &= ff['PRIORITY_INIT'] <= maxp
        seld &= ff['PRIORITY_ASSIGNED'] <= maxp
        printlog('length after cutting locations with priority_assigned > '+str(maxp)+': '+str(len(ff[seld])),logger)
    if '.ran' in fout:
        #seld &= ff['ZPOSSLOC'] == 1
        #print('length after cutting locations where target type could not be observed: '+str(len(ff[seld])))
        #selp9 = ff['PRIORITY'] == 999999
        #ff['PRIORITY'][selp9] = -999999
        seld &= ff['PRIORITY'] <= maxp
        printlog('length after cutting locations with priority > '+str(maxp)+': '+str(len(ff[seld])),logger)


    ff = ff[seld]
    
    if reccircmasks is not None:
        for maskfn in reccircmasks:
            mask = maskcircandrec(ff,maskfn,logger=logger)
            ff = ff[~mask]
    
    if ebits is not None:
        printlog('number before imaging mask '+str(len(ff)),logger)
        if ebits == 'lrg_mask':
            sel = ff['lrg_mask'] == 0
            ff = ff[sel]
        else:
            ff = cutphotmask(ff,ebits,logger=logger)
        printlog('number after imaging mask '+str(len(ff)),logger)

    if zmask:
        whz = ff['Z'] < 1.6
        ff = ff[whz]

        fzm = fitsio.read('/global/homes/m/mjwilson/desi/DX2DROPOUT/radial_mask.fits')
        zma = []
        for z in ff['Z']:
            zind = int(z/1e-6)
            zma.append(fzm[zind]['RADIAL_MASK'])
        zma = np.array(zma)
        wm = zma == 0
        ff = ff[wm]

    
    if '.dat' in fout:
        ff['Z'].name = 'Z_not4clus'
        printlog('updating completeness',logger)
        compa = []
        fractl = []
        tll = []
        ti = 0
        ff.sort('TILES')
        nts = len(np.unique(ff['TILES']))
        tlsl = ff['TILES']
        tlslu = np.unique(tlsl)
        laa = ff['LOCATION_ASSIGNED']
        lta = ff['TILELOCID_ASSIGNED']
        #print('TILELOCID_ASSIGNED',np.unique(ff['TILELOCID_ASSIGNED'],return_counts=True),len(ff))

        # for tls in np.unique(dz['TILES']): #this is really slow now, need to figure out a better way
        i = 0
        tot = 0
        atot = 0
        tltot = 0
        while i < len(ff):
            tls = []
            tlis = []
            nli = 0 #initialize total available per tile group
            nai = 0 #initialize total assigned
            nti = 0 #initialize total at location where something of the same type was assigned

            while tlsl[i] == tlslu[ti]:
                nli += 1
                nai += laa[i] #laa is true/false assigned
                nti += lta[i] #lta is true/false something of the same type was assigned
                i += 1
                if i == len(ff):
                    break

            if ti % 10000 == 0:
                printlog('at tiles ' + str(ti) + ' of ' + str(nts),logger)

            tot += nli
            atot += nai
            tltot += nti
            cp = nai / nli #
            fract = nti/nli
            # print(tls,cp,no,nt)
            compa.append(cp)
            fractl.append(fract)
            tll.append(tlslu[ti])
            ti += 1
        #print(tot,atot,tltot)
        comp_dicta = dict(zip(tll, compa))
        fract_dicta = dict(zip(tll, fractl))
        tlobs_fn = fout.replace('full'+mapveto+'.dat.fits','frac_tlobs.fits')
        tlobs = Table()
        tlobs['TILES'] = tll
        tlobs['FRAC_TLOBS_TILES'] = fractl
        write_LSS_scratchcp(tlobs,tlobs_fn,logger=logger)
        del tlobs
        fcompa = []
        fracta = []
        for tl in ff['TILES']:
            fcompa.append(comp_dicta[tl])
            fracta.append(fract_dicta[tl])
        ff['COMP_TILE'] = np.array(fcompa)
        ff['FRAC_TLOBS_TILES'] = np.array(fracta)
        printlog('data quantities measured, moving to write-out phase',logger)
        #print(np.sum(ff['FRAC_TLOBS_TILES']),len(ff))
        #if comp_only:
        #    return True
    
    if '.ran' in fout:
        printlog('area is ' + str(len(ff) / 2500),logger)
    #comments = ["'full' LSS catalog without after vetos for priority, good hardware and imaging quality","entries are for targetid that showed up in POTENTIAL_ASSIGNMENTS"]
    if wo == 'y':
        write_LSS_scratchcp(ff, fout,logger=logger)#, comments)
    if '.dat' in fout:
        wz = ff['ZWARN'] != 999999
        wz &= ff['ZWARN'] * 0 == 0
        wz &= ff['ZWARN'] != 1.e20
        comp = len(ff[wz])/len(ff)
        printlog('assignment completeness is '+str(comp),logger)
        printlog('sum of 1/(FRACZ_TILELOCID*FRAC_TLOBS_TILES), 1/COMP_TILE, and length of input; should approximately match',logger)
        printlog(str(np.sum(1. / (ff[wz]['FRACZ_TILELOCID']*ff[wz]['FRAC_TLOBS_TILES'])))+','+ str(np.sum(1. / ff[wz]['COMP_TILE']))+','+str( len(ff)),logger)
    if fout is None or wo == 'n':
        return ff
    del ff
#     tmpfn = fout+'.tmp'
#     if os.path.isfile(tmpfn):
#         os.system('rm '+tmpfn)
#     fd = fitsio.FITS(tmpfn, "rw")
#     fd.write(np.array(ff),extname='LSS')
#     fd['LSS'].write_history("created (or over-written) on "+datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
#     fd.close()
#     os.system('mv '+tmpfn+' '+fout)
    #ff.write(fout,overwrite=True,format='fits')

def apply_map_veto(fin,fout,mapn,maps,mapcuts,nside=256,logger=None):
    din = fitsio.read(fin)
    din = apply_map_veto_arrays(din,mapn,maps,mapcuts,nside=nside,logger=logger)
#     mask = np.ones(len(din),dtype='bool')
#     if 'PHOTSYS' not in list(din.dtype.names):
#         din = addNS(Table(din))
#     seln = din['PHOTSYS'] == 'N'
#     
#         
#     import healpy as hp
#     th,phi = radec2thphi(din['RA'],din['DEC'])
#     pix = hp.ang2pix(nside,th,phi,nest=True)
#     maps2cut = list(mapcuts.keys())
#     inlen = len(din)
#     print('initial',inlen)
#     for mp in maps2cut:
#         mvals = np.zeros(len(din))
#         if 'DEPTH' in mp:
#             bnd = mp.split('_')[-1]
#             mvals[seln] = mapn[mp][pix[seln]]*10**(-0.4*ext_coeff[bnd]*mapn['EBV'][pix[seln]])
#             mvals[~seln] = maps[mp][pix[~seln]]*10**(-0.4*ext_coeff[bnd]*maps['EBV'][pix[~seln]])
#             mask &= mvals > mapcuts[mp]
#             print(mp,len(din[mask]),len(din[mask])/inlen)
#             
#         else:
#             mvals[seln] = mapn[mp][pix[seln]]
#             if len(mvals[seln]) > 0:
#                 print(np.min(mvals[seln]),np.max(mvals[seln]))
#             mvals[~seln] = maps[mp][pix[~seln]]
#             print(np.min(mvals[~seln]),np.max(mvals[~seln]))
#             if mp == 'STARDENS':
#                 mvals = np.log10(mvals)
#             mask &= mvals < mapcuts[mp]   
#             print(mp,len(din[mask]),len(din[mask])/inlen)
#    write_LSS(din[mask],fout) 
    write_LSS_scratchcp(din,fout,logger=logger) 
 
def apply_map_veto_arrays(din,mapn,maps,mapcuts,nside=256,logger=None):
    mask = np.ones(len(din),dtype='bool')
    if 'PHOTSYS' not in list(din.dtype.names):
        din = addNS(Table(din))
    seln = din['PHOTSYS'] == 'N'
    
        
    import healpy as hp
    th,phi = radec2thphi(din['RA'],din['DEC'])
    pix = hp.ang2pix(nside,th,phi,nest=True)
    maps2cut = list(mapcuts.keys())
    inlen = len(din)
    printlog('initial '+str(inlen),logger)
    for mp in maps2cut:
        mvals = np.zeros(len(din))
        if 'DEPTH' in mp:
            bnd = mp.split('_')[-1]
            mvals[seln] = mapn[mp][pix[seln]]*10**(-0.4*ext_coeff[bnd]*mapn['EBV'][pix[seln]])
            mvals[~seln] = maps[mp][pix[~seln]]*10**(-0.4*ext_coeff[bnd]*maps['EBV'][pix[~seln]])
            mask &= mvals > mapcuts[mp]
            printlog(mp+','+str(len(din[mask]))+','+str(len(din[mask])/inlen),logger)
            
        else:
            mvals[seln] = mapn[mp][pix[seln]]
            if len(mvals[seln]) > 0:
                printlog(str(np.min(mvals[seln]))+','+str(np.max(mvals[seln])),logger=logger)
            mvals[~seln] = maps[mp][pix[~seln]]
            printlog(str(np.min(mvals[~seln]))+','+str(np.max(mvals[~seln])),logger)
            if mp == 'STARDENS':
                mvals = np.log10(mvals)
            mask &= mvals < mapcuts[mp]   
            printlog(mp+','+str(len(din[mask]))+','+str(len(din[mask])/inlen),logger)
    return din[mask]

            

def get_tlcomp(fin):
    '''
    fin is the full path of the catalog to use 
    '''
    ff = Table(fitsio.read(fin))#+'full_noveto.'+dr+'.fits')
    print('getting completeness')
    compa = []
    fractl = []
    tll = []
    ti = 0
    ff.sort('TILES')
    nts = len(np.unique(ff['TILES']))
    tlsl = ff['TILES']
    tlslu = np.unique(tlsl)
    laa = ff['LOCATION_ASSIGNED']
    lta = ff['TILELOCID_ASSIGNED']
    print('TILELOCID_ASSIGNED',np.unique(ff['TILELOCID_ASSIGNED'],return_counts=True),len(ff))

    # for tls in np.unique(dz['TILES']): #this is really slow now, need to figure out a better way
    i = 0
    tot = 0
    atot = 0
    tltot = 0
    while i < len(ff):
        tls = []
        tlis = []
        nli = 0 #initialize total available per tile group
        nai = 0 #initialize total assigned
        nti = 0 #initialize total at location where something of the same type was assigned

        while tlsl[i] == tlslu[ti]:
            nli += 1
            nai += laa[i] #laa is true/false assigned
            nti += lta[i] #lta is true/false something of the same type was assigned
            i += 1
            if i == len(ff):
                break

        if ti % 10000 == 0:
            print('at tiles ' + str(ti) + ' of ' + str(nts))

        tot += nli
        atot += nai
        tltot += nti
        cp = nai / nli #
        fract = nti/nli
        # print(tls,cp,no,nt)
        compa.append(cp)
        fractl.append(fract)
        tll.append(tlslu[ti])
        ti += 1
    #print(tot,atot,tltot)
    comp_dicta = dict(zip(tll, compa))
    fract_dicta = dict(zip(tll, fractl))
    tlobs_fn = fin.replace('full.dat.fits','frac_tlobs.fits')
    tlobs = Table()
    tlobs['TILES'] = tll
    tlobs['FRAC_TLOBS_TILES'] = fractl
    write_LSS(tlobs,tlobs_fn)



def write_LSS(ff, outf, comments=None,extname='LSS'):
    '''
    ff is the structured array/Table to be written out as an LSS catalog
    outf is the full path to write out
    comments is a list of comments to include in the header
    '''
    import shutil
    ranstring = int(np.random.random()*1e10)
    tmpfn = outf+'.tmp'#os.getenv('SCRATCH')+'/'+outf.split('/')[-1] + '.tmp'+str(ranstring)
    if os.path.isfile(tmpfn):
        os.system('rm ' + tmpfn)
    fd = fitsio.FITS(tmpfn, "rw")
    fd.write(np.array(ff), extname=extname)
    if comments is not None:
        for comment in comments:
            fd[extname].write_comment(comment)
    #fd[extname].write_history("updated on " + datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
    fd.close()
    print('closed fits file')
    #shutil.move(tmpfn, outf)
    #os.rename(tmpfn, outf)
    testcol = list(ff.dtype.names)[0]
    try:
        fitsio.read(tmpfn,columns=(testcol))
    except:
        print('read failed, output corrupted?!')
        return 'FAILED'    
    os.system('mv ' + tmpfn + ' ' + outf) #for some reason shutil is giving people permission issues but mv does not for actually getting the file in place
    #os.system('chmod 775 ' + outf) #this should fix permissions for the group
    print('moved output to ' + outf)
    return True

def write_LSS_scratchcp(ff, outf, comments=None,extname='LSS',logger=None):
    '''
    ff is the structured array/Table to be written out as an LSS catalog
    outf is the full path to write out
    comments is a list of comments to include in the header
    this will write to a temporary file on scratch and then copy it, then delete the temporary file once verify a successful copy
    '''
    import shutil
    printlog('will write to '+outf,logger)
    rng = np.random.default_rng()#seed=rann)
    ranstring = int(rng.random()*1e10)
    tmpfn = os.getenv('SCRATCH')+'/'+outf.split('/')[-1] + '.tmp'+str(ranstring)
    if os.path.isfile(tmpfn):
        os.system('rm ' + tmpfn)
    fd = fitsio.FITS(tmpfn, "rw")
    fd.write(np.array(ff), extname=extname)
    if comments is not None:
        for comment in comments:
            fd[extname].write_comment(comment)
    #fd[extname].write_history("updated on " + datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
    fd.close()
    if logger is None:
        print('closed fits file '+tmpfn)
    else:
        logger.info('closed fits file '+tmpfn)
    #shutil.move(tmpfn, outf)
    #os.rename(tmpfn, outf)
    testcol = list(ff.dtype.names)[0]
    try:
        fitsio.read(tmpfn,columns=(testcol))
    except:
        if logger is None:
            print('read failed, output corrupted?! '+tmpfn)
        else:
            logger.info('read failed, output corrupted?! '+tmpfn)
        return 'FAILED'    
    os.system('cp ' + tmpfn + ' ' + outf) 
    os.system('chmod 775 ' + outf) #this should fix permissions for the group
    if logger is None:
        print('moved output to ' + outf)
    else:
        logger.info('moved output to ' + outf)
    df = 0
    if logger is None:
        print('checking read of column ' + testcol)
    else:
        logger.info('checking read of column ' + testcol)

    try:
        fitsio.read(outf.replace('global','dvs_ro'),columns=(testcol))
        df = 1
    except:
        print('read failed, copy failed?! check temporary file '+tmpfn)
        return 'FAILED'    
    if logger is None:
        print('removing temp file')
    else:
        logger.info('removing temp file')

    if df == 1:
        os.system('rm '+tmpfn)
    return True


def create_sky_targets(dirname, columns=None, format_output='fits', release='1.1.1', version='main', program='dark', dr='dr9', nfiles=10, mpicomm=None):
    """
    It is a real nightmare to work with all the fits file in hpdirname_list more than 800 for each directory. Rewrite it in a smaller number of files.

    To create the sky target files:

    /* #my_code.py
    from mpi4py import MPI
    mpicomm = MPI.COMM_WORLD
    create_sky_targets('somewhere/on/nersc', nfiles=1, mpicomm=mpicomm)
    */

    On cori:
    salloc -N 1 -C haswell -t 02:00:00 --qos interactive -L SCRATCH,project
    srun -n 64 python my_code.py

    """
    from glob import glob
    import mpytools as mpy
    from mpi4py import MPI
    from fiberassign.fba_launch_io import get_desitarget_paths

    if mpicomm is None:
        print('This code works under MPI paradigm throught cosmodesi/mpytools Software. Please pass a mpi communicator as input: mpicomm=MPI.COMM_WORLD')
        sys.exit(2)

    if columns is None:
        columns = ["RA", "DEC", "TARGETID", "DESI_TARGET", "BGS_TARGET", "MWS_TARGET", "SUBPRIORITY", "OBSCONDITIONS", "PRIORITY_INIT", "NUMOBS_INIT"]

    # Now load the sky target files.  These are main-survey files that we will
    # force to be treated as the survey type of the other target files.
    mydirs = get_desitarget_paths(release, version, program, dr=dr)
    skydirs = [mydirs["sky"]]
    if os.path.isdir(mydirs["skysupp"]):
        skydirs.append(mydirs["skysupp"])

    for hpdirname in skydirs:
        fns = []
        if mpicomm.rank == 0:
            fns = glob(os.path.join(hpdirname, "*fits"))
        fns = list(mpy.bcast(fns, mpiroot=0))

        # Create mlpytools.Catalog
        targets = mpy.Catalog.read(fns, mpicomm=mpicomm)

        if format_output == 'fits':
            # first save it in nfiles fits files.
            # Create filenames to write targets in less files.
            basename = '-'.join(os.path.basename(fns[0]).split('-')[:-2])
            fns_to_save = [os.path.join(dirname, f'{basename}-{i}.fits') for i in range(nfiles)]

            start = MPI.Wtime()
            targets[columns].write(fns_to_save, filetype='fits')
            mpicomm.Barrier()
            if mpicomm.rank == 0: print(f'Write {basename} in {nfiles} done in {MPI.Wtime() - start:2.2f} s.')
        elif format_output == 'bigfile':
            # Save it on one bigfile (more efficient ?)
            start = MPI.Wtime()
            targets[columns].write(os.path.join(dirname, f'{basename}'), filetype='bigfile')
            mpicomm.Barrier()
            if mpicomm.rank == 0: print(f'Write {basename} in {nfiles} done in {MPI.Wtime() - start:2.2f} s.')
        else:
            print('Only fits or bigfile format are extected...')
            sys.exit(2)


def combtiles_pa_wdup(tiles, fbadir, outdir, tarf, addcols=['TARGETID', 'RA', 'DEC'], fba=True, tp='dark', ran='ran'):
    if ran == 'dat':
        # addcols.append('PRIORITY')
        addcols.append('PRIORITY_INIT')
        addcols.append('DESI_TARGET')
    s = 0
    td = 0
    # tiles.sort('ZDATE')
    print(len(tiles))
    outf = outdir + '/' + ran + 'comb_' + tp + 'wdup.fits'
    if fba:
        pa_hdu = 'FAVAIL'
    tl = []
    for tile in tiles['TILEID']:
        if fba:
            ffa = fbadir + '/fba-' + str(tile).zfill(6) + '.fits'
        if os.path.isfile(ffa):
            fa = Table(fitsio.read(ffa, ext=pa_hdu))
            sel = fa['TARGETID'] >= 0
            fa = fa[sel]
            td += 1
            fa['TILEID'] = int(tile)
            tl.append(fa)
            print(td, len(tiles))
        else:
            print('did not find ' + ffa)
    dat_comb = vstack(tl)
    print(len(dat_comb))
    tar_in = fitsio.read(tarf, columns=addcols)
    dat_comb = join(dat_comb, tar_in, keys=['TARGETID'])
    print(len(dat_comb))

    dat_comb.write(outf, format='fits', overwrite=True)
    print('wrote ' + outf)
    return dat_comb

def combtiles_wdup_altmtl(pa_hdu, tiles, fbadir, outf, tarf, addcols=['TARGETID', 'RA', 'DEC'],par='n',logger=None):
    s = 0
    td = 0
    print('size of tiles', len(tiles))
    tl = []
    
    tids = fitsio.read(tarf,columns=['TARGETID'])['TARGETID']
    
    def _get_fa(tile):
        fadate = return_altmtl_fba_fadate(tile)
        ffa = os.path.join(fbadir, fadate, 'fba-'+str(tile).zfill(6)+'.fits')
        if pa_hdu == 'FAVAIL':
            fa = Table(fitsio.read(ffa, ext=pa_hdu))
            sel = np.isin(fa['TARGETID'],tids)
            fa = fa[sel] #for targets, we only want science targets
        else:
            tar_hdu = 'FTARGETS'
            fa = Table(fitsio.read(ffa,ext=pa_hdu,columns=['TARGETID','LOCATION']))
            ft = Table(fitsio.read(ffa,ext=tar_hdu,columns=['TARGETID','PRIORITY','SUBPRIORITY']))
            sel = fa['TARGETID'] >= 0
            fa = fa[sel]
            lb4join = len(fa)
            #td += 1
            fa['TILEID'] = int(tile)
            fa = join(fa,ft,keys=['TARGETID'])
            if len(fa) != lb4join:
                print(tile,lb4join,len(fa))
        sel = fa['TARGETID'] >= 0
        fa = fa[sel]
        #td += 1
        fa['TILEID'] = int(tile)
        return fa
    tls = tiles['TILEID']
    if par == 'n':
        for tile in tiles['TILEID']:
            fa = _get_fa(tile)
            tl.append(fa)
    if par == 'y':
        #doesn't seem to work within function
        from concurrent.futures import ProcessPoolExecutor
        
        with ProcessPoolExecutor() as executor:
            for fa in executor.map(_get_fa, list(tls)):
                tl.append(fa)
        
    dat_comb = vstack(tl)
    printlog('size combitles for ' + pa_hdu+' , '+str(len(dat_comb)),logger=logger)
    tar_in = fitsio.read(tarf, columns=addcols)
    dat_comb = join(dat_comb, tar_in, keys=['TARGETID'],join_type='left')
    print(len(dat_comb))

    write_LSS_scratchcp(dat_comb,outf,logger=logger)
    #dat_comb.write(outf, format='fits', overwrite=True)
    #print('wrote ' + outf)
    return dat_comb


def combtiles_assign_wdup(tiles,fbadir,outdir,tarf,addcols=['TARGETID','RSDZ','TRUEZ','ZWARN'],fba=True,tp='dark'):

    s = 0
    td = 0
    #tiles.sort('ZDATE')
    print(len(tiles))
    outf = outdir+'/datcomb_'+tp+'assignwdup.fits'
    if fba:
        pa_hdu = 'FASSIGN'
        tar_hdu = 'FTARGETS'
    tl = []
    for tile in tiles['TILEID']:
        if fba:
            ffa = fbadir+'/fba-'+str(tile).zfill(6)+'.fits'
        if os.path.isfile(ffa):
            fa = Table(fitsio.read(ffa,ext=pa_hdu,columns=['TARGETID','LOCATION']))
            ft = Table(fitsio.read(ffa,ext=pa_hdu,columns=['TARGETID','PRIORITY','SUBPRIORITY']))
            sel = fa['TARGETID'] >= 0
            fa = fa[sel]
            lb4join = len(fa)
            td += 1
            fa['TILEID'] = int(tile)
            fa = join(fa,ft,keys=['TARGETID'])
            if len(fa) != lb4join:
                print(tile,lb4join,len(fa))
            tl.append(fa)
            print(td,len(tiles))
        else:
            print('did not find '+ffa)
    dat_comb = vstack(tl)
    print(len(dat_comb))
    tar_in = fitsio.read(tarf,columns=addcols)
    dat_comb = join(dat_comb,tar_in,keys=['TARGETID'])
    print(len(dat_comb))

    dat_comb.write(outf,format='fits', overwrite=True)
    print('wrote '+outf)
    return dat_comb

def addNS(tab):
    '''
    given a table that already includes RA,DEC, add PHOTSYS column denoting whether
    the data is in the DECaLS ('S') or BASS/MzLS ('N') photometric region
    '''
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    c = SkyCoord(tab['RA']* u.deg,tab['DEC']* u.deg,frame='icrs')
    gc = c.transform_to('galactic')
    sel_ngc = gc.b > 0

    #wra = (tab['RA'] > 100-tab['DEC'])
    #wra &= (tab['RA'] < 280 +tab['DEC'])
    tab['PHOTSYS'] = 'S'
    seln = tab['DEC'] > 32.375
    seln &= sel_ngc#wra
    tab['PHOTSYS'][seln] = 'N'
    return tab

def return_altmtl_fba_fadate(tileid):
    ts = str(tileid).zfill(6)
    FAOrigName = '/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'+ts[:3]+'/fiberassign-'+ts+'.fits.gz'
    fhtOrig = fitsio.read_header(FAOrigName)
    fadate = fhtOrig['RUNDATE']
    return ''.join(fadate.split('T')[0].split('-'))

def return_hp_givenradec(nside, ra, dec):
    theta, phi = np.radians(90-dec), np.radians(ra)
    return hp.ang2pix(nside, theta, phi, nest=True)

def printlog(message,logger):
    if logger is not None:
        logger.info(message)
    else:
        print(message)


def compute_wntmp(bw, prob_obs, loc_assgn, ntile, ntile_range=[0,15]):
    """
    Takes as an input the "full" sample (bitweights, prob_obs, location_assigned, ntile)
    and returns the ntile-missing-power (NTMP) factor as function of ntile.
    The routine also returns a similar factor, obtained by counting the zero-probability objects,
    which has not been used so far, as NTMP seems to give better performaces.
    In the current implementation the bitweights are almost redundant (they are only used to obtain
    the total number of bits). With minimal modification (e.g. comment/uncomment a few lines below)
    the bitewights can be used directly to compute IIPs, recurrences, etc, thus making prob_obs redundant.
    """

    nbits = 64 * np.shape(bw)[1]
    # aboo = unpack_bitweights(bw)
    # print(np.shape(aboo))
    # recurr = np.sum(aboo_full, axis=1)    
    recurr = prob_obs*nbits
    wiip = (nbits+1)/(recurr+1)
    zerop_msk = (recurr==0) & (~loc_assgn)
    
    #print(np.sum(zerop_msk))
    
    idx_nt = [] 
    idx_nt_zp = []
    idx_nt_la = []
    for n in range(ntile_range[0], ntile_range[1]+1):
        idx_nt.append(np.where(ntile == n))
        idx_nt_zp.append(np.where((ntile == n) & (zerop_msk)))  
        idx_nt_la.append(np.where((ntile == n) & (loc_assgn)))
        
    f_ntzp = np.ones(ntile_range[1]+1)
    f_ntmisspw = np.ones(ntile_range[1]+1)
    for n in range(ntile_range[0], ntile_range[1]+1):
        n_nt = len(idx_nt[n][0])
        n_ntzp = len(idx_nt_zp[n][0])
        n_ntwiip = np.sum(wiip[idx_nt_la[n][0]])
        n_ntmisspw = n_nt - n_ntwiip    
        #print(n, n_nt, n_ntzp, n_ntmisspw)
        if n_nt > 0: f_ntzp[n] = n_ntzp / n_nt
        if n_nt > 0: f_ntmisspw[n] = n_ntmisspw / n_nt   

    return f_ntmisspw, f_ntzp


def apply_wntmp(ntile, f_ntmisspw, f_ntzp, ntile_range=[0,15], randoms=True):
    """
    Takes as an input ntile from sample of interest ("full" or "clustering")
    plus the f_NTMP (and f_NTZP) factors produced by the "compute_wntmp" routine and
    returns the corresponding w_NTMP weights.
    In addition the routine returns ntile-zero-probability (NTZP) weights, not used so far,
    as NTMP seems to give better performaces.
    The the logical variable "randoms" (set to True by default) determines wheater the weights
    are for the randoms or the data by putting the factor f either at the numerator or the denominator)
    """

    idx_nt = [] 
    for n in range(ntile_range[0], ntile_range[1]+1):
        idx_nt.append(np.where(ntile == n))
        
    w_ntzp = np.ones(len(ntile))
    w_ntmisspw = np.ones(len(ntile))

    if randoms:
        for n in range(ntile_range[0], ntile_range[1]+1):
            w_ntzp[idx_nt[n][0]] = 1-f_ntzp[n]
            w_ntmisspw[idx_nt[n][0]] = 1-f_ntmisspw[n]
    else:  
        for n in range(ntile_range[0], ntile_range[1]+1):
            if 1-f_ntzp[n] > 0: w_ntzp[idx_nt[n][0]] = 1/(1-f_ntzp[n])
            if 1-f_ntmisspw[n] > 0: w_ntmisspw[idx_nt[n][0]] = 1/(1-f_ntmisspw[n])
        
    return w_ntmisspw, w_ntzp
