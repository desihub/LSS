import fitsio
import numpy as np
import sys
import os
try:
    dirpcadw = os.environ['CSCRATCH']+'/pcadw/'
    dirpc = os.environ['CSCRATCH']+'/paircounts/'
except:
    print('NEED TO BE ON NERSC, ARE YOU?')

from matplotlib import pyplot as plt

def P2(mu):
    return .5*(3.*mu**2.-1.)
    
def P4(mu):
    return .125*(35.*mu**4.-30.*mu**2.+3.)

def P6(mu):
    return 1./16.*(231.*mu**6.-315.*mu**4.+105.*mu**2.-5.)

def P8(mu):
    return 1./128.*(6435.*mu**8.-12012.*mu**6.+6930.*mu**4.-1260.*mu**2.+35.)

dirxi = os.environ['CSCRATCH']+'/SV1xi/'

om = 0.31

def prep4czxi(type,zmin,zmax,nran=10,indir='',ver='test',outdir=os.environ['CSCRATCH']+'/cz/',tile='alltiles',subset='deep',fkp=False,ranwt1=False,subt=None):
    '''
    prepare catalogs to be used by Cheng Zhao's paircount code
    '''
    fkpw = ''
    if fkp:
        fkpw = 'fkp'

    df = fitsio.read(indir+'/'+ver+'/'+type+tile+'_'+subset+'_clustering.dat.fits')
    if subt is not None:
        from desitarget.sv1 import sv1_targetmask
        tb = sv1_targetmask.desi_mask[subt]
        sel = (df['SV1_DESI_TARGET'] & tb) > 0
        df = df[sel]
        
    so = 'SV1_'+ver+type+fkpw+str(zmin)+str(zmax)
    ifiled = outdir+'g'+so+'4xi.dat'
    fo = open(ifiled,'w')
    w = (df['Z'] > zmin) & (df['Z'] < zmax) #& (df['NTILE'] > mintile)
    df = df[w]
    wt = df['WEIGHT']
    if fkp:
        wt *= df['WEIGHT_FKP']

    for i in range(0,len(df)):
        fo.write(str(df['RA'][i])+' '+str(df['DEC'][i])+' '+str(df['Z'][i])+' '+str(wt[i])+'\n')
    fo.close()
    
    ifiler = outdir+'r'+so+'4xi.dat'
    fo = open(ifiler,'w')
    for nr in range(0,nran):
        df = fitsio.read(indir+'/'+ver+'/'+type+tile+'_'+subset+'_'+str(nr)+'_clustering.ran.fits')
        
        if subt is not None:
            sel = (df['SV1_DESI_TARGET'] & tb) > 0
            df = df[sel]
        
        w = (df['Z'] > zmin) & (df['Z'] < zmax) #& (df['NTILE'] > mintile)

        df = df[w]
        if ranwt1:
            wt = np.ones(len(df))
        else:    
            wt = df['WEIGHT']
        if fkp:
            wt *= df['WEIGHT_FKP']
        print('maximum random weight is '+str(np.max(wt)))  

        for i in range(0,len(df)):
            fo.write(str(df['RA'][i])+' '+str(df['DEC'][i])+' '+str(df['Z'][i])+' '+str(wt[i])+'\n')
    fo.close()
    dirczpc = outdir + 'paircounts/'
    froot = dirczpc+so
    cf = '../Sandbox/czxi/fcfc_smu.conf'
    ddf = froot+'.dd'
    drf = froot+'.dr'
    rrf = froot+'.rr'
    fo = open('czpc.sh','w')
    fo.write('#!/bin/bash\n')
    fo.write('/global/u2/z/zhaoc/programs/FCFC_2D/2pcf -c '+cf+' -d '+ifiled+' -r '+ifiler+' --data-z-min='+str(zmin)+' --data-z-max='+str(zmax)+' --rand-z-min='+str(zmin)+' --rand-z-max='+str(zmax)+' --dd='+ddf+' --dr='+drf+' --rr='+rrf+' -p 7 -f')
    fo.close()

def calcxi_dataCZ(type,zmin,zmax,dirczpc = os.environ['CSCRATCH']+'/cz/paircounts/',fa='',bs=5,start=0,rec='',mumin=0,mumax=1,mupow=0,ver='test',fkp=False,rxp=50):
    fkpw = ''
    if fkp:
        fkpw = 'fkp'

    so = 'SV1_'+ver+type+fkpw+str(zmin)+str(zmax)

    froot = dirczpc+so
    if rec == '':

        dd = np.loadtxt(froot+'.dd').transpose()[-1]#*ddnorm
        dr = np.loadtxt(froot+'.dr').transpose()[-1]#*drnorm
        rr = np.loadtxt(froot+'.rr').transpose()[-1]

    if rec == '_rec':

        #fn += '_rec'
        dd = np.loadtxt(indir+fn+'.dd').transpose()[-1]#*ddnorm
        dr = np.loadtxt(indir+fn+'.ds').transpose()[-1]#*drnorm
        ss = np.loadtxt(indir+fn+'.ss').transpose()[-1]#*rrnorm 
        rr = np.loadtxt(indir+fnnorec+'.rr').transpose()[-1]    

    nb = (200-start)//bs
    xil = np.zeros(nb)
    xil2 = np.zeros(nb)
    xil4 = np.zeros(nb)

    nmub = 120
    dmu = 1./float(nmub)
    mubm = 0
    if mumin != 0:
        mubm = int(mumin*nmub)
    mubx = nmub
    if mumax != 1:
        mubx = int(mumax*nmub)
    for i in range(start,nb*bs+start,bs):
        xib = 0
        xib2 = 0
        xib4 = 0
        ddt = 0
        drt = 0
        rrt = 0
        sst = 0
        w = 0
        w2 = 0
        w4 = 0
        mut = 0
        rmin = i
        rmax = rmin+bs
        for m in range(mubm,mubx):
            ddb = 0
            drb = 0
            rrb = 0
            ssb = 0
            mu = m/float(nmub) + 0.5/float(nmub)
            for b in range(0,bs):
                bin = nmub*(i+b)+m
                if bin < 24000:
                    ddb += dd[bin]
                    drb += dr[bin]
                    rrb += rr[bin]
                    if rec == '_rec' or rec == 'shuff':
                        ssb += ss[bin]
                        sst += ss[bin]
                ddt += dd[bin]
                drt += dr[bin]
                rrt += rr[bin]
            xi = 0
            if rrb > 0:
                if rec == '_rec' or rec == 'shuff':
                    xi = (ddb-2.*drb+ssb)/rrb
                else:       
                    xi = (ddb-2.*drb+rrb)/rrb
            else:
                print('rrb=0 at mu '+str(mu)+' s '+str((rmin+rmax)/2.))     

            xib += xi*dmu*(mu**mupow)
            xib2 += xi*dmu*P2(mu)*5.
            xib4 += xi*dmu*P4(mu)*9.        
        xil[i//bs] = xib
        xil2[i//bs] = xib2
        xil4[i//bs] = xib4
    muw = ''
    fo = open(dirxi+'xi024'+so+rec+muw+str(bs)+'st'+str(start)+fa+'.dat','w')
    rl = []
    for i in range(0,len(xil)):
        r = bs/2.+i*bs+start
        rl.append(r)
        fo.write(str(r)+' '+str(xil[i])+' '+str(xil2[i])+' '+str(xil4[i])+'\n')
    fo.close()
    indx = int(rxp/bs)
    rl = np.array(rl)
    plt.plot(rl[:indx],rl[:indx]**2.*xil[:indx],'k-')
    plt.xlabel(r'$s$ (Mpc/h)')
    plt.ylabel(r'$s^2\xi$')
    plt.show()
    return True



def createSourcesrd_ad(sample,tile,date,ii=0,zmin=.5,zmax=1.1,datadir=''):
    '''
    prepare sv files for paircounts
    '''
    #from healpix import healpix,radec2thphi
    from random import random
    from LSS.Cosmo import distance
    d = distance(om,1.-om) #cosmology assumed in final BOSS analyses, make sure this conforms to current
    #h = healpix()
    file = sample+tile+'_'+date

    fd = fitsio.read(datadir+file+'_clustering.dat.fits')
    wz = (fd['Z'] > zmin) & (fd['Z'] < zmax)
    zw = '_zm'+str(zmin)+'zx'+str(zmax)
    zl = fd[wz]['Z']

    cdl = np.zeros(len(zl))
    for i in range(0,len(zl)):
        cdl[i] = d.dc(zl[i])    
    sral = np.sin(np.radians(fd[wz]['RA']))
    cral = np.cos(np.radians(fd[wz]['RA']))
    sdecl = np.sin(np.radians(fd[wz]['DEC']))
    cdecl = np.cos(np.radians(fd[wz]['DEC']))
    wl = fd[wz]['WEIGHT']
    print(str(len(cdl))+' data objects going out for paircounts')
    gf ='g'+file+zw
    fdo = open(dirpcadw+gf +'pcadw.dat','w')
    for i in range(0,len(cdl)):
        fdo.write(str(sral[i])+' '+str(cral[i])+' '+str(sdecl[i])+' '+str(cdecl[i])+' '+str(cdl[i])+' '+str(wl[i])+'\n')

    fdo.close()

    fr = fitsio.read(datadir+file+'_'+str(ii)+'_clustering.ran.fits')
    wz = (fr['Z'] > zmin) & (fr['Z'] < zmax)

    zl = fr[wz]['Z']

    cdl = np.zeros(len(zl))
    for i in range(0,len(zl)):
        cdl[i] = d.dc(zl[i])    
    sral = np.sin(np.radians(fr[wz]['RA']))
    cral = np.cos(np.radians(fr[wz]['RA']))
    sdecl = np.sin(np.radians(fr[wz]['DEC']))
    cdecl = np.cos(np.radians(fr[wz]['DEC']))
    wl = np.ones(len(cdl))
    print(str(len(cdl))+' random objects going out for paircounts')
    rf = 'r'+file+str(ii)+zw
    fdo = open(dirpcadw+rf +'pcadw.dat','w')
    for i in range(0,len(cdl)):
        fdo.write(str(sral[i])+' '+str(cral[i])+' '+str(sdecl[i])+' '+str(cdecl[i])+' '+str(cdl[i])+' '+str(wl[i])+'\n')

    fdo.close()
    fo = open('dopc'+gf+'.sh','w')
    fo.write('#!/bin/bash\n')
    fo.write('/global/homes/a/ajross/pp2pt_Dmufb '+gf +' '+gf +' \n')#while in AJR's home directory, I think permissions allow it...?
    fo.write('/global/homes/a/ajross/pp2pt_Dmufb '+gf +' '+rf +' \n')
    fo.write('/global/homes/a/ajross/pp2pt_Dmufb '+rf +' '+rf +' \n')
    fo.close()

    
    return gf

def createSourcesrd_ari(sample,tile,date,ii,zmin=.5,zmax=1.1,datadir=''):
    '''
    prepare sv files for paircounts
    '''
    #from healpix import healpix,radec2thphi
    from random import random
    from LSS.Cosmo import distance
    d = distance(om,1.-om) #cosmology assumed in final BOSS analyses, make sure this conforms to current
    #h = healpix()
    file = sample+tile+'_'+date

    fr = fitsio.read(datadir+file+'_'+str(ii)+'_clustering.ran.fits')
    wz = (fr['Z'] > zmin) & (fr['Z'] < zmax)

    zl = fr[wz]['Z']

    cdl = np.zeros(len(zl))
    for i in range(0,len(zl)):
        cdl[i] = d.dc(zl[i])    
    sral = np.sin(np.radians(fr[wz]['RA']))
    cral = np.cos(np.radians(fr[wz]['RA']))
    sdecl = np.sin(np.radians(fr[wz]['DEC']))
    cdecl = np.cos(np.radians(fr[wz]['DEC']))
    wl = np.ones(len(cdl))
    print(str(len(cdl))+' random objects going out for paircounts')
    zw = '_zm'+str(zmin)+'zx'+str(zmax)
    gf ='g'+file+zw
    rf = 'r'+file+str(ii)+zw
    fdo = open(dirpcadw+rf +'pcadw.dat','w')
    for i in range(0,len(cdl)):
        fdo.write(str(sral[i])+' '+str(cral[i])+' '+str(sdecl[i])+' '+str(cdecl[i])+' '+str(cdl[i])+' '+str(wl[i])+'\n')

    fdo.close()
    fo = open('dopc'+rf+'.sh','w')
    fo.write('#!/bin/bash\n')
    fo.write('/global/homes/a/ajross/pp2pt_Dmufb '+gf +' '+rf +' \n')
    fo.write('/global/homes/a/ajross/pp2pt_Dmufb '+rf +' '+rf +' \n')
    fo.close()

    return rf

def ppxilcalc_LSDfjack_bs(sample,tile,date,zmin=.5,zmax=1.1,bs=1,start=0,rmaxf=250,rmax=50,mumin=0,mumax=1.,wmu='',mom=0,nran=1,vis=False):
    fl = sample+tile+'_'+date+'_zm'+str(zmin)+'zx'+str(zmax)
    flr0 = sample+tile+'_'+date+'0_zm'+str(zmin)+'zx'+str(zmax)
    DDnl = []   
    DDnorml = 0
    DDnormt = 0
    DRnl = []
    DRnorml = 0
    DRnormt = 0
    RRnl = []
    RRnl0 = []
    nmubin = 100
    nbin = rmax/bs
    for i in range(0,rmaxf*nmubin):
        DDnl.append(0)
        DRnl.append(0)
        RRnl.append(0)
        RRnl0.append(0)
    RRnorml = 0
    RRnormt = 0
    pl = []
    nmut = 0
    for i in range(0,nmubin):
        mu = i/float(nmubin)+.5/float(nmubin)
        mub = int(mu*nmubin)
        ##print mu
        if mu > mumin and mu < mumax:
            pl.append((1.,P2(mu),P4(mu),P6(mu),P8(mu),mub))
            nmut += 1.
        else:
            pl.append((0,0,0,0,0,0))    
    fdp = open(dirpc+'g'+fl+'g'+fl+'2ptdmu.dat').readlines()
    DDnormt += float(fdp[0])
    fdnp = open(dirpc+'g'+fl+'r'+flr0+'2ptdmu.dat').readlines()
    fr = open(dirpc+'r'+flr0+'r'+flr0+'2ptdmu.dat').readlines()
    DRnormt += float(fdnp[0])
    RRnormt += float(fr[0])
    #print(DRnormt,RRnormt)
    for k in range(1,len(fdp)):
        dp = float(fdp[k])
        dr = float(fdnp[k])
        rp = float(fr[k])
        DDnl[k-1] += dp
        DRnl[k-1] += dr
        RRnl[k-1] += rp
    for n in range(1,nran):
        flrn = sample+tile+'_'+date+str(n)+'_zm'+str(zmin)+'zx'+str(zmax)
        fdnp = open(dirpc+'g'+fl+'r'+flrn+'2ptdmu.dat').readlines()
        fr = open(dirpc+'r'+flrn+'r'+flrn+'2ptdmu.dat').readlines()
        DRnormt += float(fdnp[0])
        RRnormt += float(fr[0])
        #print(DRnormt,RRnormt)
        for k in range(1,len(fdp)):
            dr = float(fdnp[k])
            rp = float(fr[k])
            DRnl[k-1] += dr
            RRnl[k-1] += rp
    #print(DDnormt,DRnormt,RRnormt)                
    xil = np.zeros(int(nbin),'f')
    for i in range(start,rmax,bs):
        xi = 0
        dd = 0
        dr = 0
        rr = 0

        ddt = 0
        drt = 0
        rrt = 0
        nmunz = 0
        for j in range(0,nmubin):
            if wmu != 'counts':
                dd = 0
                dr = 0
                rr = 0
            for k in range(0,bs):
                bin = nmubin*(i+k)+j            
                if bin < len(RRnl):
                    #if RRnl[bin] == 0:
                    #   pass
        
                    #else:
                    dd += DDnl[bin]
                    rr += RRnl[bin]
                    dr += DRnl[bin]
                    ddt +=dd
                    rrt += rr
                    drt += dr
        
            #if rr != 0 and wm == 'muw':            
            if wmu != 'counts' and rr != 0:
                xi += pl[j][mom]/float(nmut)*(dd/DDnormt-2*dr/DRnormt+rr/RRnormt)*RRnormt/rr
                nmunz += 1.
        if wmu == 'counts':
            xi = (dd/DDnormt-2*dr/DRnormt+rr/RRnormt)*RRnormt/rr        
        else:
            xi *= 100./nmunz #correct for empty mu bins
            if nmunz != 100.:
                print('there were only '+str(nmunz)+' mu bins with rr counts in bin '+str(i)+' to '+str(i+bs)+' mpc/h')
        if i/bs < nbin:
            xil[i//bs] = xi
        #print(ddt/DDnormt,drt/DRnormt,rrt/RRnormt)
    rl = []
    for i in range(0,len(xil)):
        rl.append(start+bs/2.+bs*i)
    rl = np.array(rl)
    if vis:
        plt.plot(rl,xil)
        plt.show()
    bsst = str(bs)+'st'+str(start)
    if wmu == 'counts':
        outf = dirxi+'xi'+fl+bsst+'.dat'
        
    else:
        outf = dirxi+'xi'+str(2*mom)+fl+bsst+'.dat'
    fo = open(outf,'w')    
    for i in range(0,len(rl)):
        fo.write(str(rl[i])+' '+str(xil[i])+'\n')
    fo.close()  
    print('wrote results to '+outf)  
    return xil

def plotxi(xidir=''):
    dl = np.loadtxt(xidir+'xiLRG70003_20200219_zm0.5zx1.1bsc.dat').transpose() 
    de = np.loadtxt(xidir+'xiELG70004_20200219_zm0.8zx1.6bsc.dat').transpose()
    ml = (dl[0]/7.78)**-1.98
    plt.loglog(dl[0],dl[1],'r-',label='MINI_SV_LRG, Tile 70003')
    plt.loglog(de[0],de[1],'b-',label='MINI_SV_ELG, Tile 70004')
    plt.loglog(dl[0],ml,'r:',label=r'$(r/7.78)^{-1.98}$ (Kitanidis et al.)')
    plt.legend()
    plt.xlabel(r'$r$ ($h^{-1}$Mpc)')
    plt.ylabel(r'$\xi$')
    plt.title('From 202200219')
    plt.savefig(xidir+'miniSVxi.png')
    plt.show()

def plot3ELG(xidir=''):
    d1 = np.loadtxt(xidir+'minisvxi/xiELG70004_20200219_zm0.8zx1.61st0.dat').transpose() 
    d2 = np.loadtxt(xidir+'xiELG70005_20200228_zm0.8zx1.61st0.dat').transpose()
    d3 = np.loadtxt(xidir+'xiELG70006_20200303_zm0.8zx1.61st0.dat').transpose()
    plt.loglog(d1[0],d1[1],label='MINI_SV_ELG, Tile 70004')
    plt.loglog(d2[0],d2[1],label='MINI_SV_ELG, Tile 70005')
    plt.loglog(d3[0],d3[1],label='MINI_SV_ELG, Tile 70006')
    dm = (d1[1]+d2[1]+d3[1])/3.
    plt.loglog(d3[0],dm,'k-',label='MINI_SV_ELG, mean 70004,5,6')
    plt.legend()
    plt.ylim(1.e-3,30)
    plt.show()

def plotELG0(xidir=''):
    d1 = np.loadtxt(xidir+'xiELG67142_20200315_zm0.8zx1.61st0.dat').transpose() 
    d2 = np.loadtxt(xidir+'xiELG67230_20200315_zm0.8zx1.61st0.dat').transpose()
    plt.loglog(d1[0],d1[1],label='SV0_ELG, Tile 67142')
    plt.loglog(d2[0],d2[1],label='SV0_ELG, Tile 67230')
    dm = (d1[1]+d2[1])/2.
    plt.loglog(d2[0],dm,'k-',label='SV0_ELG, mean')
    plt.legend()
    plt.ylim(1.e-3,30)
    plt.show()


def plot2LRG(xidir=''):
    d1 = np.loadtxt(xidir+'xiLRG70002_20200304_zm0.5zx1.11st0.dat').transpose() 
    d2 = np.loadtxt(xidir+'xiLRG70003_20200228_zm0.5zx1.11st0.dat').transpose()
    plt.loglog(d1[0],d1[1],label='MINI_SV_LRG, Tile 70002')
    plt.loglog(d2[0],d2[1],label='MINI_SV_LRG, Tile 70003')
    dm = (d1[1]+d2[1])/2.
    plt.loglog(d2[0],dm,'k-',label='MINI_SV_LRG, mean 70002,3')
    plt.legend()
    plt.ylim(1.e-3,30)
    plt.show()

def plotLRG0(xidir=''):
    d1 = np.loadtxt(xidir+'xiLRG68001_20200315_zm0.5zx1.11st0.dat').transpose() 
    d2 = np.loadtxt(xidir+'xiLRG68002_20200315_zm0.5zx1.11st0.dat').transpose()
    plt.loglog(d1[0],d1[1],label='SV0_LRG, Tile 68001')
    plt.loglog(d2[0],d2[1],label='SV0_LRG, Tile 68002')
    dm = (d1[1]+d2[1])/2.
    plt.loglog(d2[0],dm,'k-',label='SV0_LRG, mean')
    plt.legend()
    plt.ylim(1.e-3,30)
    plt.show()

def plotQSO0(xidir=''):
    d1 = np.loadtxt(xidir+'xiQSO68001_20200315_zm0.8zx2.25st0.dat').transpose() 
    d2 = np.loadtxt(xidir+'xiQSO68002_20200315_zm0.8zx2.25st0.dat').transpose()
    plt.loglog(d1[0],d1[1],label='SV0_QSO, Tile 68001')
    plt.loglog(d2[0],d2[1],label='SV0_QSO, Tile 68002')
    dm = (d1[1]+d2[1])/2.
    plt.loglog(d2[0],dm,'k-',label='SV0_QSO, mean')
    plt.legend()
    plt.ylim(1.e-3,30)
    plt.show()


def plotxicomb(xidir=''):
    xilin = np.loadtxt(xidir+'xi0Challenge_matterpower0.42.04915.00.dat').transpose()
    d1 = np.loadtxt(xidir+'xiELG70004_20200219_zm0.8zx1.6bsc.dat').transpose() 
    d2 = np.loadtxt(xidir+'xiELG70005_20200228_zm0.8zx1.6bsc.dat').transpose()
    d3 = np.loadtxt(xidir+'xiELG70006_20200303_zm0.8zx1.6bsc.dat').transpose()
    dme = (d1[1]+d2[1]+d3[1])/3.

    
    plt.loglog(d1[0],dme,'b-',label='MINI_SV_ELG, mean 70004,5,6')
    d1 = np.loadtxt(xidir+'xiLRG70002_20200304_zm0.5zx1.1bsc.dat').transpose() 
    d2 = np.loadtxt(xidir+'xiLRG70003_20200228_zm0.5zx1.1bsc.dat').transpose()
    dml = (d1[1]+d2[1])/2.
    ml = (d1[0]/7.78)**-1.98
    plt.loglog(d1[0],dml,'r-',label='MINI_SV_LRG, mean 70002,3')
    plt.loglog(d1[0],ml,'r:',label=r'$(r/7.78)^{-1.98}$ (Kitanidis et al.)')
    plt.legend()
    plt.xlabel(r'$r$ ($h^{-1}$Mpc)')
    plt.ylabel(r'$\xi$')
    plt.savefig(xidir+'miniSVxicomb.png')
    plt.show()

def plotxicomb_dec(bs=1,xidir=''):
    xilin = np.loadtxt(xidir+'xi0Challenge_matterpower0.42.04915.00.dat').transpose()
    d1 = np.loadtxt(xidir+'xiELG80606_20201216_zm0.8zx1.6'+str(bs)+'st0.dat').transpose() 
    d2 = np.loadtxt(xidir+'xiELG80608_20201215_zm0.8zx1.6'+str(bs)+'st0.dat').transpose()
    d3 = np.loadtxt(xidir+'xiELG80610_20201216_zm0.8zx1.6'+str(bs)+'st0.dat').transpose()
    dme = (d1[1]+d2[1]+d3[1])/3.    
    plt.loglog(d1[0],dme,'b-',label='DEC_SV1_ELG, mean 80606,08,10')

    d1 = np.loadtxt(xidir+'xiLRG80605_20201215_zm0.5zx1.1'+str(bs)+'st0.dat').transpose() 
    d2 = np.loadtxt(xidir+'xiLRG80607_20201215_zm0.5zx1.1'+str(bs)+'st0.dat').transpose()
    d3 = np.loadtxt(xidir+'xiLRG80609_20201216_zm0.5zx1.1'+str(bs)+'st0.dat').transpose()
    dml = (d1[1]+d2[1]+d3[1])/3.
    ml = (d1[0]/7.78)**-1.98
    plt.loglog(d1[0],dml,'r-',label='DEC_SV1_LRG, mean 80605,7,9')

    plt.loglog(d1[0],ml,'r--',label=r'$(r/7.78)^{-1.98}$ (Kitanidis et al.)')
    plt.loglog(xilin[0],xilin[1]*1.4,'r:',label=r'3.1$\xi_{\rm lin}(z=0.8)$')
    plt.loglog(xilin[0],xilin[1]*.7,'b:',label=r'2$\xi_{\rm lin}(z=1.1)$')
    plt.legend()
    plt.xlabel(r'$r$ ($h^{-1}$Mpc)')
    plt.ylabel(r'$\xi$')
    plt.ylim(1.e-2,70)
    plt.xlim(0.25,60)
    plt.savefig(xidir+'xicomb'+str(bs)+'.png')
    plt.show()

def plotxicomb_deep(bs=1,xidir=''):
    xilin = np.loadtxt(os.environ['HOME']+'/BAOtemplates/xi0Challenge_matterpower0.42.04915.00.dat').transpose()
    d1 = np.loadtxt(xidir+'xiELG80606_deep_zm0.8zx1.6'+str(bs)+'st0.dat').transpose() 
    d2 = np.loadtxt(xidir+'xiELG80608_deep_zm0.8zx1.6'+str(bs)+'st0.dat').transpose()
    d3 = np.loadtxt(xidir+'xiELG80610_deep_zm0.8zx1.6'+str(bs)+'st0.dat').transpose()
    d4 = np.loadtxt(xidir+'xiELG80621_deep_zm0.8zx1.6'+str(bs)+'st0.dat').transpose()
    d5 = np.loadtxt(xidir+'xiELG80623_deep_zm0.8zx1.6'+str(bs)+'st0.dat').transpose()
    dme = (d1[1]+d2[1]+d3[1]+d4[1]+d5[1])/5.    
    plt.loglog(d1[0],dme,'b-',label='DEEP_SV1_ELG, mean 80606,08,10,21,23')

    d1 = np.loadtxt(xidir+'xiLRG80605_deep_zm0.5zx1.1'+str(bs)+'st0.dat').transpose() 
    d2 = np.loadtxt(xidir+'xiLRG80607_deep_zm0.5zx1.1'+str(bs)+'st0.dat').transpose()
    d3 = np.loadtxt(xidir+'xiLRG80609_deep_zm0.5zx1.1'+str(bs)+'st0.dat').transpose()
    d4 = np.loadtxt(xidir+'xiLRG80620_deep_zm0.5zx1.1'+str(bs)+'st0.dat').transpose()
    d5 = np.loadtxt(xidir+'xiLRG80622_deep_zm0.5zx1.1'+str(bs)+'st0.dat').transpose()
    dml = (d1[1]+d2[1]+d3[1]+d4[1]+d5[1])/5.
    ml = (d1[0]/7.78)**-1.98
    plt.loglog(d1[0],dml,'r-',label='DEC_SV1_LRG, mean 80605,7,9,20,22')

    plt.loglog(d1[0],ml,'r--',label=r'$(r/7.78)^{-1.98}$ (Kitanidis et al.)')
    plt.loglog(xilin[0],xilin[1]*1.4,'r:',label=r'3.1$\xi_{\rm lin}(z=0.8)$')
    plt.loglog(xilin[0],xilin[1]*.7,'b:',label=r'2$\xi_{\rm lin}(z=1.1)$')
    plt.legend()
    plt.xlabel(r'$r$ ($h^{-1}$Mpc)')
    plt.ylabel(r'$\xi$')
    plt.ylim(1.e-2,70)
    plt.xlim(0.25,60)
    plt.savefig(xidir+'xicomb'+str(bs)+'.png')
    plt.show()


def plotxicomb0(xidir=''):
    xilin = np.loadtxt(xidir+'xi0Challenge_matterpower0.42.04915.00.dat').transpose()
    d1 = np.loadtxt(xidir+'xiELG67142_20200315_zm0.8zx1.61st0.dat').transpose() 
    d2 = np.loadtxt(xidir+'xiELG67230_20200315_zm0.8zx1.61st0.dat').transpose()
    dme = (d1[1]+d2[1])/2.  
    plt.loglog(d1[0],dme,'b-',label='SV0_ELG, mean 67142,67230')

    d1 = np.loadtxt(xidir+'xiLRG68001_20200315_zm0.5zx1.11st0fid.dat').transpose() 
    d2 = np.loadtxt(xidir+'xiLRG68002_20200315_zm0.5zx1.11st0fid.dat').transpose()
    dml = (d1[1]+d2[1])/2.
    ml = (d1[0]/7.78)**-1.98
    plt.loglog(d1[0],dml,'r-',label='SV0_LRG, mean 68001,68002')
    plt.loglog(d1[0],ml,'r--',label=r'$(r/7.78)^{-1.98}$ (Kitanidis et al.)')
    plt.loglog(xilin[0],xilin[1]*1.4,'r:',label=r'3.1$\xi_{\rm lin}(z=0.8)$')
    plt.loglog(xilin[0],xilin[1]*.7,'b:',label=r'2$\xi_{\rm lin}(z=1.1)$')
    plt.legend()
    plt.xlabel(r'$r$ ($h^{-1}$Mpc)')
    plt.ylabel(r'$\xi$')
    plt.ylim(1.e-2,70)
    plt.xlim(0.25,60)
    plt.savefig(xidir+'miniSV0xicomb.png')
    plt.show()


if __name__ == '__main__':
    import subprocess
    minisvdir = '/project/projectdirs/desi/users/ajross/catalogs/minisv2/'
    datadir = minisvdir+'LSScats/'

    night = '20200315'
    type = 'LRG'
    tile = '68001'
#   gf = createSourcesrd_ad(type,tile,night)
#   subprocess.run(['chmod','+x','dopc'+gf+'.sh'])
#   subprocess.run('./dopc'+gf+'.sh')
#   ppxilcalc_LSDfjack_bs(type,tile,night,zmin=.5,zmax=1.1)
#   ppxilcalc_LSDfjack_bs(type,tile,night,zmin=.5,zmax=1.1,bs=5)

    type = 'QSO'
    gf = createSourcesrd_ad(type,tile,night,zmin=.8,zmax=2.2,datadir=datadir)
    subprocess.run(['chmod','+x','dopc'+gf+'.sh'])
    subprocess.run('./dopc'+gf+'.sh')
    ppxilcalc_LSDfjack_bs(type,tile,night,zmin=.8,zmax=2.2)
    ppxilcalc_LSDfjack_bs(type,tile,night,zmin=.8,zmax=2.2,bs=5)
    tile = '68002'

#   type = 'LRG'
#   gf = createSourcesrd_ad(type,tile,night)
#   subprocess.run(['chmod','+x','dopc'+gf+'.sh'])
#   subprocess.run('./dopc'+gf+'.sh')
#   ppxilcalc_LSDfjack_bs(type,tile,night,zmin=.5,zmax=1.1)
#   ppxilcalc_LSDfjack_bs(type,tile,night,zmin=.5,zmax=1.1,bs=5)

    type = 'QSO'
    gf = createSourcesrd_ad(type,tile,night,zmin=.8,zmax=2.2,datadir=datadir)
    subprocess.run(['chmod','+x','dopc'+gf+'.sh'])
    subprocess.run('./dopc'+gf+'.sh')
    ppxilcalc_LSDfjack_bs(type,tile,night,zmin=.8,zmax=2.2)
    ppxilcalc_LSDfjack_bs(type,tile,night,zmin=.8,zmax=2.2,bs=5)

#   type = 'ELG'
#   tile = '67230'
#   gf = createSourcesrd_ad(type,tile,night,zmin=.8,zmax=1.6)
#   subprocess.run(['chmod','+x','dopc'+gf+'.sh'])
#   subprocess.run('./dopc'+gf+'.sh')
#   ppxilcalc_LSDfjack_bs(type,tile,night,zmin=.8,zmax=1.6)
#   ppxilcalc_LSDfjack_bs(type,tile,night,zmin=.8,zmax=1.6,bs=5)
#   tile = '67142'
#   gf = createSourcesrd_ad(type,tile,night,zmin=.8,zmax=1.6)
#   subprocess.run(['chmod','+x','dopc'+gf+'.sh'])
#   subprocess.run('./dopc'+gf+'.sh')
#   ppxilcalc_LSDfjack_bs(type,tile,night,zmin=.8,zmax=1.6)
#   ppxilcalc_LSDfjack_bs(type,tile,night,zmin=.8,zmax=1.6,bs=5)
