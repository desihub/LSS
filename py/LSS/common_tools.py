import numpy as np
import fitsio
from astropy.table import Table,join

from LSS.tabulated_cosmo import TabulatedDESI
cosmo = TabulatedDESI()
dis_dc = cosmo.comoving_radial_distance



#functions that shouldn't have any dependence on survey go here

def find_znotposs(dz):

    dz.sort('TARGETID')
    tidnoz = []
    tids = np.unique(dz['TARGETID'])
    ti = 0
    i = 0
    
    print('finding targetids that were not observed')
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
            print(ti)
        ti += 1 

    
    selnoz = np.isin(dz['TARGETID'],tidnoz)
    tidsb = np.unique(dz[selnoz]['TILELOCID'])
    #dz = dz[selnoz]
    dz.sort('TILELOCID')
    tids = np.unique(dz['TILELOCID'])
    print('number of targetids with no obs '+str(len(tidnoz)))
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
            print(ti,len(tids))
        ti += 1 
    #the ones to veto are now the join of the two
    wtbtlid = np.isin(tlidnoz,tidsb)
    tlidnoz = np.array(tlidnoz)
    lznposs = tlidnoz[wtbtlid]
    print('number of locations where assignment was not possible because of priorities '+str(len(lznposs)))
    return lznposs


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

def addnbar(fb,nran=18,bs=0.01,zmin=0.01,zmax=1.6):
    '''
    fb is the root of the file name, including the path
    nran is the number of random files to add the nz to 
    bs is the bin size of the nz file (read this from file in future)
    zmin is the lower edge of the minimum bin (read this from file in future)
    zmax is the upper edge of the maximum bin (read this from file in the future)
    '''
    
    nzd = np.loadtxt(fb+'_nz.dat').transpose()[3] #column with nbar values
    fn = fb+'_clustering.dat.fits'
    fd = fitsio.read(fn) #reading in data with fitsio because it is much faster to loop through than table
    zl = fd['Z']
    nl = np.zeros(len(zl))
    for ii in range(0,len(zl)):
        z = zl[ii]
        zind = int((z-zmin)/bs)
        if z > zmin and z < zmax:
            nl[ii] = nzd[zind]
    del fd
    ft = Table.read(fn)
    ft['NZ'] = nl
    ft.write(fn,format='fits',overwrite=True)        
    print('done with data')
    for rann in range(0,nran):
        fn = fb+'_'+str(rann)+'_clustering.ran.fits'
        fd = fitsio.read(fn) #reading in data with fitsio because it is much faster to loop through than table
        zl = fd['Z']
        nl = np.zeros(len(zl))
        for ii in range(0,len(zl)):
            z = zl[ii]
            zind = int((z-zmin)/bs)
            if z > zmin and z < zmax:
                nl[ii] = nzd[zind]
        del fd
        ft = Table.read(fn)
        ft['NZ'] = nl
        ft.write(fn,format='fits',overwrite=True)      
        print('done with random number '+str(rann))  
    return True        
