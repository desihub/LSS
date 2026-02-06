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

def expand_ran(in_ran_fn,parent_ran_fn,in_clus_fileroot,rancols=['TARGETID','RA','DEC'],datacols=['TARGETID','Z'],logger=None):
    #function to add columns to randoms, most useful for mock randoms where the same column values are used
    #assumes data is saved in the LSS h5 format; could edit to allow functionality for fits or other formats
    '''
    Parameters
    ----------
    in_ran_fn : string 
        full path to input randoms to expand, in LSS h5 format, containing at least the columns 'TARGETID','TARGETID_DATA','WEIGHT','NX'
    parent_ran_fn : string 
        full path to randoms that are a superset of data to expand and contain at least TARGETID, RA, DEC, in LSS h5 format
    in_clus_fileroot : string
        directory + tracer for clustering catalogs (allowing NGC/SGC to be loaded by simply adding +'_'+reg+'_clustering.dat.h5'), in LSS h5 format
    rancols : list 
        a list of the column names to add to the input table from the parent random catalog via TARGETID match
    datacols : list 
        a list of the column names to add to the input table from the data catalog via a TARGETID_DATA to TARGETID match
    Returns
    ----------
    A table in astropy format for the randoms containing the desired columns
    '''
    
    #t0 = time.time()
    parent_ran = read_hdf5_blosc(parent_ran_fn,columns=rancols)
    in_table = read_hdf5_blosc(in_ran_fn,columns=['TARGETID','TARGETID_DATA','WEIGHT','NX'])
    #tran = time.time()
    #print(str(rann)+' read original randoms;'+str(tran-t0))
    
    olen = len(in_table)
    tids, in_ind, orig_ind = np.intersect1d(in_table['TARGETID'], parent_ran['TARGETID'], return_indices=True)
    in_table = in_table[in_ind]
    #print(np.array_equal(tids,in_table['TARGETID']))
    parent_ran = parent_ran[orig_ind]
    for col in rancols:
        if col != 'TARGETID':
            in_table[col] = parent_ran[col]
    #in_table = join(in_table,in_ran,keys=['TARGETID']) #astropy join is much slower
    #t1 = time.time()
    #print(str(rann)+' joined to original randoms;'+str(t1-t0))
    del parent_ran
    regl = ['NGC','SGC']
    datal = []
    for reg in regl:
        datal.append(read_hdf5_blosc(in_clus_fileroot+'_'+reg+'_clustering.dat.h5',columns=datacols))
    in_data = vstack(datal)
    #t2 = time.time()
    #print(str(rann)+' stacked data;'+str(t2-t0))
    del datal    
    in_data.rename_column('TARGETID', 'TARGETID_DATA')
    #in_table = join(in_table,in_data,keys=['TARGETID_DATA']) #astropy join is much slower
    
    if in_data['TARGETID_DATA'].max() < 1000000000:
        lookup = np.arange(1 + in_data['TARGETID_DATA'].max())
        lookup[in_data['TARGETID_DATA']] = np.arange(len(in_data))
        indices = lookup[in_table['TARGETID_DATA']]
    else:
        printlog('TARGETID_DATA max is too large, using slower method',None)
        sorted_idx = np.argsort(in_data['TARGETID_DATA'])
        idx_in_sorted = np.searchsorted(in_data['TARGETID_DATA'], in_table['TARGETID_DATA'], sorter=sorted_idx)
        indices = sorted_idx[idx_in_sorted]
    for col in datacols:
        if col != 'TARGETID':
            in_table[col] = in_data[col][indices]

    #t3 = time.time()
    #print(str(rann)+' done;'+str(t3-t0))
    #print(olen,len(in_table))
    return in_table

def mask_bad_fibers_time_dependent(dz, badfibers_td,logger=None):
    # Use the time-dependent bad fiber file to remove data
    # dz is the input data file, which must have FIBER and LASTNIGHT columns

    #badfibers_td = open('/global/cfs/cdirs/desi/survey/catalogs/DA2/LSS/loa-v1/unique_badfibers_time-dependent.txt','r').readlines()

    inds_to_remove = np.array([])
    for i in range(len(badfibers_td)):
        #print(i)
    
        if len(badfibers_td[i].split()) == 1:
            inds_to_remove = np.concatenate((inds_to_remove, (np.where(dz['FIBER'] == int(badfibers_td[i]))[-1])))
        elif len(badfibers_td[i].split()) == 2:
            inds_to_remove = np.concatenate((inds_to_remove, (np.where((dz['FIBER'] == int(badfibers_td[i].split()[0]))
                & (dz['LASTNIGHT'] >= int(badfibers_td[i].split()[1]))
                )[-1])))
            #print(np.where((dz['FIBER'] == int(badfibers_td[i].split()[0])
            #   & (dz['LASTNIGHT'] >= int(badfibers_td[i].split()[1]))))
            #print(5/0)
        elif len(badfibers_td[i].split()) == 3:
            inds_to_remove = np.concatenate((inds_to_remove, (np.where((dz['FIBER'] == int(badfibers_td[i].split()[0]))
                & (dz['LASTNIGHT'] >= int(badfibers_td[i].split()[1])) & 
                (dz['LASTNIGHT'] < int(badfibers_td[i].split()[2]))
                )[-1])))
            #print((np.where((dz['FIBER'] == int(badfibers_td[i].split()[0])
            #   & (dz['LASTNIGHT'] >= int(badfibers_td[i].split()[1])) & 
            #   (dz['LASTNIGHT'] < int(badfibers_td[i].split()[2]))
            #   ))
        elif len(badfibers_td[i].split()) == 4:
            inds_to_remove = np.concatenate((inds_to_remove, (np.where((dz['FIBER'] == int(badfibers_td[i].split()[0]))
                & (((dz['LASTNIGHT'] >= int(badfibers_td[i].split()[1])) & 
                (dz['LASTNIGHT'] < int(badfibers_td[i].split()[2])))
                | (dz['LASTNIGHT'] >= int(badfibers_td[i].split()[3]))
                ))[-1])))
            #print(np.where((dz['FIBER'] == int(badfibers_td[i].split()[0])
            #   & ((((dz['LASTNIGHT'] >= int(badfibers_td[i].split()[1])) & 
            #   (dz['LASTNIGHT'] < int(badfibers_td[i].split()[2])))
            #   | (dz['LASTNIGHT'] >= int(badfibers_td[i].split()[3]))))
            #   ))
        elif len(badfibers_td[i].split()) == 5:
            inds_to_remove = np.concatenate((inds_to_remove, (np.where((dz['FIBER'] == int(badfibers_td[i].split()[0]))
                & (((dz['LASTNIGHT'] >= int(badfibers_td[i].split()[1])) & 
                (dz['LASTNIGHT'] < int(badfibers_td[i].split()[2])))
                | ((dz['LASTNIGHT'] >= int(badfibers_td[i].split()[3]))
                & (dz['LASTNIGHT'] < int(badfibers_td[i].split()[4])))
                ))[-1])))
            #print((np.where((dz['FIBER'] == int(badfibers_td[i].split()[0])
            #   & ((dz['LASTNIGHT'] >= int(badfibers_td[i].split()[1])) & 
            #   (dz['LASTNIGHT'] < int(badfibers_td[i].split()[2])))
            #   | ((dz['LASTNIGHT'] >= int(badfibers_td[i].split()[3]))
            #   & (dz['LASTNIGHT'] < int(badfibers_td[i].split()[4])))
            #   )))
    dz_inds = np.ones(len(dz)).astype('int')
    dz_inds[inds_to_remove.astype('int')] = 0
    nremove = len(dz_inds)-np.sum(dz_inds)
    printlog('number removed from time dependent mask is '+str(nremove),logger)
    return dz[dz_inds == 1]



#functions that shouldn't have any dependence on survey go here

def cut_specdat(dz,badfib=None,tsnr_min=0,tsnr_col='TSNR2_ELG',logger=None,fibstatusbits=None,remove_badfiber_spike_nz=False,mask_petal_nights=False):
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
        if type(badfib).__name__ == 'list':
            printlog('will mask fibers with time dependence',logger)
            cat_out = mask_bad_fibers_time_dependent(fs[wfqa], badfib,logger=logger)
        else:
            bad = np.isin(fs['FIBER'],badfib)
            printlog('number at bad fibers '+str(sum(bad)),logger)
            cat_out = fs[wfqa&~bad]
    else:
        cat_out = fs[wfqa]
    if remove_badfiber_spike_nz:
        badfib1 = np.loadtxt('/global/cfs/cdirs/desi/survey/catalogs/DA2/LSS/loa-v1/bad_nz_fibers_ks_test.txt')
        badfib2 = np.loadtxt('/global/cfs/cdirs/desi/survey/catalogs/DA2/LSS/loa-v1/elg_bad_nz_spike_fibers_1.498_1.499.txt')
        bad = np.isin(cat_out['FIBER'],np.concatenate((badfib1,badfib2)))
        printlog('number removed from spike mask is '+str(np.sum(bad)),logger)
        cat_out = cat_out[~bad]
        
    if mask_petal_nights:
        cat_out = mask_bad_petal_nights(cat_out,logger=logger)
    
    return cat_out
    
def mask_bad_petal_nights(dz, prog='dark',logger=None):
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
    nremove = len(dz_inds)-np.sum(dz_inds)
    printlog('number removed from petal night mask is '+str(nremove),logger)

    return dz[dz_inds == 1]


def goodz_infull(tp,dz,zcol='Z_not4clus'):
    if (tp == 'LRG') or (tp == 'LGE'):
        z_suc= dz['ZWARN']==0
        #print(np.sum(z_suc))
        z_suc &= dz['DELTACHI2']>15
        #print(np.sum(z_suc))
        z_suc &= dz[zcol]<1.5
        #print(np.sum(z_suc))

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

def mk_zcmbmap(nside=256,nest=True):
    import healpy as hp
    pix = np.arange(12*nside*nside)
    th,phi = hp.pix2ang(nside,pix,nest=nest)
    ra,dec = thphi2radec(th,phi)
    zcmb = get_zcmbdipole(ra,dec)
    return zcmb

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
    if '.fits' in fcr:
        ranf = fitsio.read_header(fcr,ext=1) #should have originally had 2500/deg2 density, so can convert to area
        area = ranf['NAXIS2']/randens
    
        print('area is '+str(area))
        outf = open(fout,'w')
        outf.write('#area is '+str(area)+'square degrees\n')
    
        if compmd == 'ran':
            ranf = fitsio.read(fcr)
            area = np.sum(ranf['FRAC_TLOBS_TILES'])/randens
            outf.write('#effective area is '+str(area)+'square degrees\n')
    
    
    if '.h5' in fcr:
        ranf = read_hdf5_blosc(fcr) #should have originally had 2500/deg2 density, so can convert to area
        area = len(ranf)/randens
    
        print('area is '+str(area))
        outf = open(fout,'w')
        outf.write('#area is '+str(area)+'square degrees\n')
    
        if compmd == 'ran':
            area = np.sum(ranf['FRAC_TLOBS_TILES'])/randens
            outf.write('#effective area is '+str(area)+'square degrees\n')

    del ranf
    if '.fits' in fcd:
        df = fitsio.read(fcd)
    if '.h5' in fcd:
        df = read_hdf5_blosc(fcd)

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

