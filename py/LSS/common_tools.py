import numpy as np
import fitsio

from LSS.tabulated_cosmo import TabulatedDESI
cosmo = TabulatedDESI()
dis_dc = cosmo.comoving_radial_distance



#functions that shouldn't have any dependence on survey go here

def mknz(fcd,fcr,fout,bs=0.01,zmin=0.01,zmax=1.6):
    '''
    fcd is the full path to the catalog file in fits format with the data; requires columns Z and WEIGHT
    fcr is the full path to the random catalog meant to occupy the same area as the data; assumed to come from the imaging randoms that have a density of 2500/deg2
    fout is the full path to the file name
    bs is the bin width for the n(z) calculation
    zmin is the lower edge of the first bin
    zmax is the upper edge of the last bin
    '''
    #cd = distance(om,1-om)
    ranf = fitsio.read(fcr) #should have originally had 2500/deg2 density, so can convert to area
    area = len(ranf)/2500.
    print('area is '+str(area))
    
    df = fitsio.read(fcd)
    
    nbin = int((zmax-zmin)/bs)
    zhist = np.histogram(df['Z'],bins=nbin,range=(zmin,zmax),weights=df['WEIGHT'])
    outf = open(fout,'w')
    outf.write('#area is '+str(area)+'square degrees\n')
    outf.write('#zmid zlow zhigh n(z) Nbin Vol_bin\n')
    for i in range(0,nbin):
        zl = zhist[1][i]
        zh = zhist[1][i+1]
        zm = (zh+zl)/2.
        voli = area/(360.*360./np.pi)*4.*np.pi/3.*(dis_dc(zh)**3.-dis_dc(zl)**3.)
        nbarz =  zhist[0][i]/voli
        outf.write(str(zm)+' '+str(zl)+' '+str(zh)+' '+str(nbarz)+' '+str(zhist[0][i])+' '+str(voli)+'\n')
    outf.close()
