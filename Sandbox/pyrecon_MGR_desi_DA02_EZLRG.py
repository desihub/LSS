from pyrecon import MultiGridReconstruction, utils
from astropy import units as u
from astropy.coordinates import SkyCoord
from LSS.tabulated_cosmo import TabulatedDESI
from cosmoprimo.fiducial import DESI
import numpy as np
from datetime import datetime
from astropy.table import vstack,hstack
from astropy.table import Table
from astropy.io import fits
import os


for i in range(27,28) :

    start = datetime.now()
    print('Start at %s'%start)

    print('Data Catalog')

    xt,yt,zt,wt = Table(),Table(),Table(),Table()

    data_fn = '/global/project/projectdirs/desi/users/dvalcin/EZMOCKS/LRG/Mocks/LRG_z0.800_cutsky_seed%d_data.fits'%(i)
    open_data = fits.open(data_fn)
    table = Table.read(data_fn,format='fits')
    temp_xt = table['RA']
    temp_yt = table['DEC']
    temp_zt = table['Z']
    w_dt = table['WEIGHT_FKP']

    temp_x = open_data[1].data['RA']
    temp_y = open_data[1].data['DEC']
    temp_z = open_data[1].data['Z']
    w_d = open_data[1].data['WEIGHT_FKP']

    cosmo = DESI()
    cosmo = cosmo.get_background(engine='class')

    dist_d = cosmo.comoving_radial_distance(temp_z)

    c_d = SkyCoord(ra=temp_x*u.degree, dec=temp_y*u.degree, distance=dist_d*u.Mpc)

    xx = c_d.cartesian.x.value
    yy = c_d.cartesian.y.value
    zz = c_d.cartesian.z.value

    xt = vstack([xt,xx])
    yt = vstack([yt,yy])
    zt = vstack([zt,zz])
    
    pos_dat = np.vstack((xx,yy,zz)).T

    table_dat = hstack([xt,yt,zt,w_dt])

    print('Random Catalog')

    xx,yy,zz,w_r = [],[],[],[]

    xt,yt,zt,wt = Table(),Table(),Table(),Table()

    for n in range(10):
        data_fn = '/global/project/projectdirs/desi/users/dvalcin/EZMOCKS/LRG/Mocks/LRG_z0.800_cutsky_S%d00_random.fits'%(n+1)
        open_data = fits.open(data_fn)
        table = Table.read(data_fn,format='fits')
        temp_xt = table['RA']
        temp_yt = table['DEC']
        temp_zt = table['Z']
        temp_wt = table['WEIGHT_FKP']
        xt = vstack([xt,temp_xt])
        yt = vstack([yt,temp_yt])
        zt = vstack([zt,temp_zt])
        wt = vstack([wt,temp_wt])
        temp_x = open_data[1].data['RA']
        temp_y = open_data[1].data['DEC']
        temp_z = open_data[1].data['Z']
        temp_w_r = open_data[1].data['WEIGHT_FKP']

        cosmo = DESI()
        cosmo = cosmo.get_background(engine='class')

        dist_d = cosmo.comoving_radial_distance(temp_z)

        c_d = SkyCoord(ra=temp_x*u.degree, dec=temp_y*u.degree, distance=dist_d*u.Mpc)

        xx += list(c_d.cartesian.x.value)
        yy += list(c_d.cartesian.y.value)
        zz += list(c_d.cartesian.z.value)
        w_r += list(temp_w_r)

    pos_ran = np.vstack((xx,yy,zz)).T

    w_r = np.array(w_r)

    table_r = hstack([xt,yt,zt,wt])

    print('Start recon')

    recon = MultiGridReconstruction(f=0.8,bias=1.99,los='x',positions=pos_dat,fft_engine='fftw',nmesh=None,cellsize=7,dtype='f4')
    recon.assign_data(pos_dat,w_d)
    recon.assign_randoms(pos_ran,w_r)
    recon.set_density_contrast()
    recon.run()

    pos_dat_rec = recon.read_shifted_positions(pos_dat,field = 'disp+rsd')
    pos_ran_rec = recon.read_shifted_positions(pos_ran,field = 'disp')

    dist, ra, dec = utils.cartesian_to_sky(pos_dat_rec)
    dist_r, ra_r, dec_r = utils.cartesian_to_sky(pos_ran_rec)

    print('Start to write fits')

    distance_to_redshift = utils.DistanceToRedshift(TabulatedDESI().comoving_radial_distance)

    table_dat['RA'],table_dat['DEC'],table_dat['Z'] = ra, dec, distance_to_redshift(dist)
    table_r['RA'],table_r['DEC'],table_r['Z'] = ra_r, dec_r, distance_to_redshift(dist_r)
    table_dat.write(os.environ['CSCRATCH']+'/DA02_recon/reciso_MGR_LRG_z0.800_fkp_seed%d_data.fits'%(i),format='fits',overwrite=True)
    table_r.write(os.environ['CSCRATCH']+'/DA02_recon/reciso_MGR_LRG_z0.800_fkp_seed%d_random_S100_1000.fits'%(i),format='fits',overwrite=True)

    end = datetime.now()
    print('Mock %s done'%i)
    print('Took %s'%(end-start))
