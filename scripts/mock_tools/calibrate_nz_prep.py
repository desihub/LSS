import numpy as np
import os
from numpy.random import Generator, PCG64
rng = Generator(PCG64())

path_to_nz = '/pscratch/sd/a/acarnero/codes/desi-cutsky-mock/nz_files/HighFidelity'

def get_nz(z_cat, tracer_type, ns=None):
    ''' The function where the n(z) is read and the NZ column is computed for the given redshifts.   '''
    if tracer_type == 'LRG':
        nzfile = 'nz_lrg_v4.txt'
    elif tracer_type == 'QSO':
        nzfile = 'NZ_QSO_v3.txt'
    elif tracer_type == 'BGS':
        nzfile = 'NZ_QSO_v3.txt'
    elif tracer_type == 'ELG':
        ns = True
        north_file = 'nz_elg_N_v5.txt'
        south_file = 'nz_elg_S_v5.txt'
    
    if ns is None:
        print('reading nz function', os.path.join(path_to_nz, nzfile))
        z, nz = np.loadtxt(os.path.join(path_to_nz, nzfile), usecols=([0,1]), unpack=True)
        return np.interp(z_cat, z, nz, left=0, right=0)

    else:
        print('Reading North and South calibration separately from', os.path.join(path_to_nz, north_file), os.path.join(path_to_nz, south_file))
        z, nz = np.loadtxt(os.path.join(path_to_nz, north_file), usecols=([0,1]), unpack=True)
        north_interp = np.interp(z_cat, z, nz, left=0, right=0)

        z, nz = np.loadtxt(os.path.join(path_to_nz, south_file), usecols=([0,1]), unpack=True)
        south_interp = np.interp(z_cat, z, nz, left=0, right=0)

        return north_interp, south_interp

def downsample_aux(z_cat, ran, n_mean, tracer_type, val=1):
    """ downsample galaxies following n(z) model specified in galtype"""

    nz = get_nz(z_cat, tracer_type)
    if isinstance(n_mean, list) or isinstance(n_mean, np.ndarray):
        n_mean_interp = np.interp(z_cat, n_mean[0], n_mean[1], left=0, right=0)
        print('n_mean is list or array')

        # downsample
        nz_selected = ran < (nz) / (n_mean_interp)
    else:
        nz_selected = ran < (nz) / (n_mean)
    ###nz_selected = ran < (1+nz) / (1+n_mean)
    
    idx = np.where(nz_selected)

    #idx is the ones to select
    print("DOWNSAMPLE: Selected {} out of {} galaxies.".format(len(idx[0]), len(z_cat)), flush=True)

    newbits = np.zeros(len(z_cat), dtype=np.int32)
    newbits[idx] = val

    return newbits, nz
    #return idx, nz

def return_north(ra, dec):
    '''
    given a table that already includes RA,DEC, add PHOTSYS column denoting whether
    the data is in the DECaLS ('S') or BASS/MzLS ('N') photometric region
    '''
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    c = SkyCoord(ra* u.deg, dec* u.deg,frame='icrs')
    gc = c.transform_to('galactic')
    sel_ngc = gc.b > 0

    seln = dec > 32.375

    sel = seln&sel_ngc
    return sel


def downsample_aux_NS(z_cat, ran, n_mean, radec, tracer_type):
    """ downsample galaxies following n(z) model specified in galtype"""


    nz_N, nz_S = get_nz(z_cat, tracer_type, ns = True)
    ra = radec[0]
    dec = radec[1]

    is_north = return_north(ra, dec)
    is_south = ~return_north(ra, dec)

        # downsample
    print(n_mean[0])
    print(n_mean[1])
    n_mean_interpN = np.interp(z_cat, n_mean[0][0], n_mean[0][1], left=0, right=0)
    n_mean_interpS = np.interp(z_cat, n_mean[1][0], n_mean[1][1], left=0, right=0)

    nz_selected_n = (ran < nz_N / n_mean_interpN) & is_north
    nz_selected_s = (ran < nz_S / n_mean_interpS) & is_south
#        n_N = nz_N / n_mean        
#        n_S = nz_S / n_mean        
    idx = np.where(nz_selected_n|nz_selected_s)

    print("DOWNSAMPLE North and South: Selected {} out of {} galaxies.".format(len(idx[0]), len(z_cat)), flush=True)

    #bitval = bits(ask=ask)

    newbits = np.zeros(len(z_cat), dtype=np.int32)
    newbits[idx] = 1
    return newbits, nz_N, nz_S

def get_lop(z_cat, ran):
    
    lop_file = 'nz_lop_v5.txt'
    z,nz = np.loadtxt(os.path.join(path_to_nz, lop_file), usecols=([0,1]), unpack=True)
    frac_lop = np.interp(z_cat, z, nz, left=0, right=0)
    nz_selected = ran < frac_lop

    idx = np.where(nz_selected)
    newbits = np.zeros(len(z_cat), dtype=np.int32)
    print("ELG_LOP selection: Selected {} out of {} galaxies.".format(len(idx[0]), len(z_cat)), flush=True)
    newbits[idx] = 4
    return newbits



def downsample(z_cat, n_mean, tracer_type, radec = None):
    """ downsample galaxies following n(z) model specified in galtype"""

    ran_i = rng.random(len(z_cat))

    outbits = []

    if tracer_type == "LRG" or tracer_type == "QSO":
        outbits, nzgood = downsample_aux(z_cat, ran_i, n_mean, tracer_type)
        ran = [ran_i]

    elif tracer_type == "ELG":
        newbits, nz_N, nz_S = downsample_aux_NS(z_cat, ran_i, n_mean, radec, tracer_type)
        ran_n               = rng.random(len(z_cat))
        ran_n[newbits == 0] = np.inf
        newbits_LOP      = get_lop(z_cat, ran_n)

        outbits = np.bitwise_or(newbits, newbits_LOP)
        ran = [ran_i, ran_n]

    else:
        print("Wrong galaxy type.")
        os._exit(1)

    return outbits, ran



def calibrate_nz(input_data, redshift_column = 'Z_RSD', tracer_type='LRG', n_mean=None, survey='DA2'):
    print('Entering NZ calibration')

    limits = {'LRG':[0., 1.6], 'ELG':[0.,2.1], 'QSO':[0.,3.5]}
    areas = {'DA2':13176.892}
    areas_split = {'DA2': [3212.9032, 9963.9888]}  #primero north y luego south
    ra = input_data['RA'][()]
    dec = input_data['DEC'][()]
    z = input_data[redshift_column][()]

    if tracer_type != 'ELG':

#    limits[tracer_type] = [np.min(z), np.max(z)]
        if n_mean is None:
            print('n_mean not given, estimate as a function of z')
        #vol = get_volume(tracer_type, limits[tracer_type], areas[survey])

            ztarget, n_mean = mknz(input_data, areas[survey], zmax=limits[tracer_type][1], zcol=redshift_column)
            n_mean = [ztarget, n_mean]
        #n_mean = np.array(n_mean)
        #print('volume per tracer is', vol)
        #n_mean = len(ra[(z>=limits[tracer_type][0])*(z<=limits[tracer_type][1])])/vol
        #print('mean density in working mock is', n_mean)    


        down_bit, ran_arr = downsample(z, n_mean, tracer_type, radec=[ra, dec])
    else:
        print('IS ELG')
        if n_mean is None:
            maskN = return_north(input_data['RA'], input_data['DEC'])

            ztargetN, n_meanN = mknz(input_data[maskN], areas_split[survey][0], zmax=limits[tracer_type][1], zcol=redshift_column)
            n_meanN = [ztargetN, n_meanN]
            print(n_meanN)
            ztargetS, n_meanS = mknz(input_data[~maskN], areas_split[survey][1], zmax=limits[tracer_type][1], zcol=redshift_column)
            n_meanS = [ztargetS, n_meanS]
            print(n_meanS)
            down_bit, ran_arr = downsample(z, [n_meanN, n_meanS], tracer_type, radec=[ra, dec])
        else:
            down_bit, ran_arr = downsample(z, n_mean, tracer_type, radec=[ra, dec])



    out_arr = down_bit.astype(np.int32)
    input_data['STATUS'] = out_arr
    return input_data

def mask_abacusHF(nz=0, foot=None, nz_lop=0):
    if foot == 'Y1':
        Y5 = 0
        Y1 = 1
        Y3 = 0
    elif foot == 'DA2':
        Y5 = 0
        Y1 = 0
        Y3 = 1
    else:
        Y5 = 1
        Y1 = 0
        Y3 = 0

    return nz * (2**0) + Y5 * (2**1) + nz_lop * (2**2) + Y1 * (2**3) + Y3 * (2**5)


def get_volume(tracer_type, zlims, area):
    from LSS.tabulated_cosmo import TabulatedDESI 
    cosmo = TabulatedDESI() 
    dis_dc = cosmo.comoving_radial_distance

    
    return area/(360.*360./np.pi)*4.*np.pi/3.*(dis_dc(zlims[1])**3.-dis_dc(zlims[0])**3.)


def mknz(df, area, bs = 0.01, zmin = 0.01, zmax = 1.6, randens = 2500., zcol='Z'):
    from LSS.tabulated_cosmo import TabulatedDESI
    cosmo = TabulatedDESI()
    dis_dc = cosmo.comoving_radial_distance
    print('area is '+str(area))


    nbin = int((zmax-zmin)*(1+bs/10)/bs)
    zhist = np.histogram(df[zcol],bins=nbin,range=(zmin,zmax)) #,weights=wts)

    zreturn, nzreturn = [], []
    for i in range(0,nbin):
        zl = zhist[1][i]
        zh = zhist[1][i+1]
        zm = (zh+zl)/2.
        voli = area/(360.*360./np.pi)*4.*np.pi/3.*(dis_dc(zh)**3.-dis_dc(zl)**3.)
        nbarz =  zhist[0][i]/voli
        zreturn.append(zm)
        nzreturn.append(nbarz)
    return zreturn, nzreturn

