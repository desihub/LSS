dirpcadw = '/global/cscratch1/sd/ajross/pcadw/'
dirpc = '/global/cscratch1/sd/ajross/paircounts/'
dirczpc = '/global/cscratch1/sd/ajross/cz/paircounts/'
dircz = '/global/cscratch1/sd/ajross/cz/'
dirxi = '/project/projectdirs/desi/users/ajross/e2exi/'
import fitsio
import numpy as np

from matplotlib import pyplot as plt

def P2(mu):
	return .5*(3.*mu**2.-1.)
	
def P4(mu):
	return .125*(35.*mu**4.-30.*mu**2.+3.)

def P6(mu):
	return 1./16.*(231.*mu**6.-315.*mu**4.+105.*mu**2.-5.)

def P8(mu):
	return 1./128.*(6435.*mu**8.-12012.*mu**6.+6930.*mu**4.-1260.*mu**2.+35.)


om = 0.31

def createSourcesrd_ad(sample,tile,date,zmin=.5,zmax=1.1):
	'''
	prepare minisv files for paircounts
	'''
	#from healpix import healpix,radec2thphi
	from random import random
	from Cosmo import distance
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

	fr = fitsio.read(datadir+file+'_clustering.ran.fits')
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
	print(str(len(cdl))+' randdom objects going out for paircounts')
	rf = 'r'+file+zw
	fdo = open(dirpcadw+rf +'pcadw.dat','w')
	for i in range(0,len(cdl)):
		fdo.write(str(sral[i])+' '+str(cral[i])+' '+str(sdecl[i])+' '+str(cdecl[i])+' '+str(cdl[i])+' '+str(wl[i])+'\n')

	fdo.close()
	fo = open('dopc'+gf+'.sh','w')
	fo.write('#!/bin/bash\n')
	fo.write('./pp2pt_Dmufb '+gf +' '+gf +' \n')
	fo.write('./pp2pt_Dmufb '+gf +' '+rf +' \n')
	fo.write('./pp2pt_Dmufb '+rf +' '+rf +' \n')
	fo.close()

	
	return gf

def ppxilcalc_LSDfjack_bs(sample,tile,date,zmin=.5,zmax=1.1,bs=1,start=0,rmaxf=250,rmax=50,mumin=0,mumax=1.,wmu='counts',mom=0):
	fl = sample+tile+'_'+date+'_zm'+str(zmin)+'zx'+str(zmax)
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
	fdnp = open(dirpc+'g'+fl+'r'+fl+'2ptdmu.dat').readlines()
	fr = open(dirpc+'r'+fl+'r'+fl+'2ptdmu.dat').readlines()
	DRnormt += float(fdnp[0])
	RRnormt += float(fr[0])
	for k in range(1,len(fdp)):
		dp = float(fdp[k])
		dr = float(fdnp[k])
		rp = float(fr[k])
		DDnl[k-1] += dp
		DRnl[k-1] += dr
		RRnl[k-1] += rp
					
	xil = np.zeros(int(nbin),'f')
	for i in range(start,rmax,bs):
		xi = 0
		dd = 0
		dr = 0
		rr = 0
	
		ddt = 0
		drt = 0
		rrt = 0
		for j in range(0,nmubin):
			if wmu != 'counts':
				dd = 0
				dr = 0
				rr = 0
			for k in range(0,bs):
				bin = nmubin*(i+k)+j			
				if bin < len(RRnl):
					#if RRnl[bin] == 0:
					#	pass
			
					#else:
					dd += DDnl[bin]
					rr += RRnl[bin]
					dr += DRnl[bin]
					ddt +=dd
					rrt += rr
					drt += dr
			
			#if rr != 0 and wm == 'muw':			
			if wmu != 'counts':
				xi += pl[j][mom]/float(nmut)*(dd/DDnormt-2*dr/DRnormt+rr/RRnormt)*RRnormt/rr
		if wmu == 'counts':
			xi = (dd/DDnormt-2*dr/DRnormt+rr/RRnormt)*RRnormt/rr		
		if i/bs < nbin:
			xil[i//bs] = xi
		#print ddt/DDnormt,drt/DRnormt,rrt/RRnormt
	rl = []
	for i in range(0,len(xil)):
		rl.append(start+bs/2.+bs*i)
	rl = np.array(rl)
	#plt.plot(rl,xil)
	#plt.show()
	bsst = str(bs)+'st'+str(start)
	fo = open('xi'+fl+bsst+'.dat','w')
	for i in range(0,len(rl)):
		fo.write(str(rl[i])+' '+str(xil[i])+'\n')
	fo.close()		
	return xil

def prep4czxi(type,zmin,zmax,program='dark',truez='',ver='g',fkp=True,test='',mintile=0,pg=-1):
	if test == 'test':
		e2edir = '/project/projectdirs/desi/users/ajross/catalogs/test/'
	else:
		e2edir ='/global/cfs/cdirs/desi/datachallenge/onepercent/catalogs/'
		#e2edir = '/global/homes/m/mjwilson/desi/survey-validation/svdc-spring2020'+ver+'-onepercent/run/catalogs/'
	fkpw = ''
	if fkp:
		fkpw = 'fkp'
	pgw = ''
	if pg != -1:
		pgw += str(pg)	
	df = fitsio.read(e2edir+program+'/'+type+'_oneper'+truez+'_clus.dat.fits')
	so = 'e2e_oneper'+ver+type+truez+fkpw+str(zmin)+str(zmax)+str(mintile)+pgw
	ifiled = dircz+'g'+so+'4xi.dat'
	fo = open(ifiled,'w')
	w = (df['Z'] > zmin) & (df['Z'] < zmax) #& (df['NTILE'] > mintile)
	if pg != -1:
		w &= (df['PROGRAM'] & 2**pg > 0)
	df = df[w]
	wt = df['WEIGHT']
	if fkp:
		wt *= df['WEIGHT_FKP']

	for i in range(0,len(df)):
		fo.write(str(df['RA'][i])+' '+str(df['DEC'][i])+' '+str(df['Z'][i])+' '+str(wt[i])+'\n')
	fo.close()
	df = fitsio.read(e2edir+program+'/'+type+'_oneper'+truez+'_clus.ran.fits')
	ifiler = dircz+'r'+so+'4xi.dat'
	fo = open(ifiler,'w')
	w = (df['Z'] > zmin) & (df['Z'] < zmax) #& (df['NTILE'] > mintile)
	if pg != -1:
		w &= (df['PROGRAM'] & 2**pg > 0)

	df = df[w]
	wt = df['WEIGHT']
	if fkp:
		wt *= df['WEIGHT_FKP']

	for i in range(0,len(df)):
		fo.write(str(df['RA'][i])+' '+str(df['DEC'][i])+' '+str(df['Z'][i])+' '+str(wt[i])+'\n')
	fo.close()
	print(dirczpc)
	froot = dirczpc+so
	cf = 'czxi/fcfc_smu.conf'
	ddf = froot+'.dd'
	drf = froot+'.dr'
	rrf = froot+'.rr'
	fo = open('czpc.sh','w')
	fo.write('#!/bin/bash\n')
	fo.write('/global/u2/z/zhaoc/programs/FCFC_2D/2pcf -c '+cf+' -d '+ifiled+' -r '+ifiler+' --data-z-min='+str(zmin)+' --data-z-max='+str(zmax)+' --rand-z-min='+str(zmin)+' --rand-z-max='+str(zmax)+' --dd='+ddf+' --dr='+drf+' --rr='+rrf+' -p 7 -f')
	fo.close()

def calcxi_dataCZ(type,zmin,zmax,truez='',bs=5,start=0,rec='',mumin=0,mumax=1,mupow=0,ver='g',fkp=True,test='',mintile=0,pg=-1):
    fkpw = ''
    if fkp:
        fkpw = 'fkp'
    pgw = ''
    if pg != -1:
        pgw += str(pg)	

    so = 'e2e_oneper'+ver+type+truez+fkpw+str(zmin)+str(zmax)+str(mintile)+pgw

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
            if rec == '_rec' or rec == 'shuff':
                xi = (ddb-2.*drb+ssb)/rrb
            else:		
                xi = (ddb-2.*drb+rrb)/rrb

            xib += xi*dmu*(mu**mupow)
            xib2 += xi*dmu*P2(mu)*5.
            xib4 += xi*dmu*P4(mu)*9.		
        xil[i//bs] = xib
        xil2[i//bs] = xib2
        xil4[i//bs] = xib4
    muw = ''
    fo = open(dirxi+'xi024'+so+rec+muw+str(bs)+'st'+str(start)+'.dat','w')
    rl = []
    for i in range(0,len(xil)):
        r = bs/2.+i*bs+start
        rl.append(r)
        fo.write(str(r)+' '+str(xil[i])+' '+str(xil2[i])+' '+str(xil4[i])+'\n')
    fo.close()
    rl = np.array(rl)
    plt.plot(rl,rl**2.*xil,'k-')
    plt.xlabel(r'$s$ (Mpc/h)')
    plt.ylabel(r'$s^2\xi$')
    plt.show()
    return True

def plotxi():
	fl = dirxi+'xi024LRG0.51.15st0.dat'
	fe = dirxi+'xi024ELG0.61.45st0.dat'
	fq = dirxi+'xi024QSO0.82.25st0.dat'
	fb = dirxi+'xi024BGS0.10.45st0.dat'
	dl = np.loadtxt(fl).transpose() 
	de = np.loadtxt(fe).transpose()
	dq = np.loadtxt(fq).transpose()
	db = np.loadtxt(fb).transpose()
	plt.plot(dl[0],dl[1]*dl[0]**2.,color='r',label=r'LRGs, $0.5 < z < 1.1$')
	plt.plot(dl[0],de[1]*dl[0]**2.,color='b',label=r'ELGs, $0.6 < z < 1.4$')
	plt.plot(dl[0],dq[1]*dl[0]**2.,color='purple',label=r'quasars, $0.8 < z < 2.2$')
	plt.plot(dl[0],db[1]*dl[0]**2.,color='brown',label=r'BGS, $0.1 < z < 0.4$')
	plt.legend()
	plt.xlabel(r'$r$ ($h^{-1}$Mpc)')
	plt.ylabel(r'$\xi_0$')
	plt.title('e2e simulation')
	plt.savefig(dirxi+'xi0e2e.png')
	plt.show()

def plotxi_compgf(type,zmin,zmax):
	fr = type+str(zmin)+str(zmax)
	ff = dirxi+'xi024'+fr+'5st0.dat'
	fg = dirxi+'xi024oneperg'+fr+'5st0.dat'
	df = np.loadtxt(ff).transpose() 
	dg = np.loadtxt(fg).transpose()
	plt.plot(df[0],df[1]*df[0]**2.,'r--',label='one per cent f')
	plt.plot(df[0],dg[1]*df[0]**2.,'-',color='firebrick',label='one per cent g')
	plt.legend()
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)')
	plt.ylabel(r'$s^2\xi_0$')
	plt.title(r'e2e simulation '+type +' '+str(zmin) +'$<z<$'+str(zmax))
	plt.savefig(dirxi+'xi0gf'+type+'.png')
	plt.show()
	plt.plot(df[0],df[2]*df[0],'r--',label='one per cent f')
	plt.plot(df[0],dg[2]*df[0],'-',color='firebrick',label='one per cent g')
	plt.legend()
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)')
	plt.ylabel(r'$s\xi_2$')
	plt.title(r'e2e simulation '+type +' '+str(zmin) +'$<z<$'+str(zmax))
	plt.savefig(dirxi+'xi2gf'+type+'.png')
	plt.show()
	plt.plot(df[0],df[3]*df[0],'r--',label='one per cent f')
	plt.plot(df[0],dg[3]*df[0],'-',color='firebrick',label='one per cent g')
	plt.legend()
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)')
	plt.ylabel(r'$s\xi_4$')
	plt.title(r'e2e simulation '+type +' '+str(zmin) +'$<z<$'+str(zmax))
	plt.savefig(dirxi+'xi4gf'+type+'.png')
	plt.show()

def plotxi_compfkp(type,zmin,zmax):
	fr = 'xi024oneperg'+type
	ff = dirxi+fr+'fkp'+str(zmin)+str(zmax)+'5st0.dat'
	fg = dirxi+fr+str(zmin)+str(zmax)+'5st0.dat'
	df = np.loadtxt(ff).transpose() 
	dg = np.loadtxt(fg).transpose()
	plt.plot(df[0],df[1]*df[0]**2.,'r-',label='FKP weight')
	plt.plot(df[0],dg[1]*df[0]**2.,'--',color='firebrick',label='no FKP weight')
	plt.legend()
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)')
	plt.ylabel(r'$s^2\xi_0$')
	plt.title(r'e2e simulation '+type +' '+str(zmin) +'$<z<$'+str(zmax))
	plt.savefig(dirxi+'xi0gf'+type+'.png')
	plt.show()
	plt.plot(df[0],df[2]*df[0],'r--',label='FKP weight')
	plt.plot(df[0],dg[2]*df[0],'-',color='firebrick',label='no FKP weight')
	plt.legend()
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)')
	plt.ylabel(r'$s\xi_2$')
	plt.title(r'e2e simulation '+type +' '+str(zmin) +'$<z<$'+str(zmax))
	plt.savefig(dirxi+'xi2gf'+type+'.png')
	plt.show()
	plt.plot(df[0],df[3]*df[0],'r--',label='FKP weight')
	plt.plot(df[0],dg[3]*df[0],'-',color='firebrick',label='no FKP weight')
	plt.legend()
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)')
	plt.ylabel(r'$s\xi_4$')
	plt.title(r'e2e one per cent simulation version g '+type +' '+str(zmin) +'$<z<$'+str(zmax))
	plt.savefig(dirxi+'xi4gf'+type+'.png')
	plt.show()


def plotxi_comptrue():
	fl = dirxi+'xi024onepergLRGfkp0.51.15st0.dat'
	flt = dirxi+'xi024onepergLRGztruefkp0.51.15st0.dat'
	fet = dirxi+'xi024onepergELGztruefkp0.61.45st0.dat'
	fe = dirxi+'xi024onepergELGfkp0.61.45st0.dat'
	fq = dirxi+'xi024onepergQSOfkp0.82.25st0.dat'
	fqt = dirxi+'xi024onepergQSOztruefkp0.82.25st0.dat'
	fb = dirxi+'xi024onepergBGSfkp0.10.45st0.dat'
	fbt = dirxi+'xi024onepergBGSztruefkp0.10.45st0.dat'
	dl = np.loadtxt(fl).transpose() 
	de = np.loadtxt(fe).transpose()
	dq = np.loadtxt(fq).transpose()
	db = np.loadtxt(fb).transpose()
	dlt = np.loadtxt(flt).transpose() 
	det = np.loadtxt(fet).transpose()
	dqt = np.loadtxt(fqt).transpose()
	dbt = np.loadtxt(fbt).transpose()
	plt.plot(dl[0],dl[1]*dl[0]**2.,color='r',label=r'LRGs, $0.5 < z < 1.1$')
	plt.plot(dl[0],de[1]*dl[0]**2.,color='b',label=r'ELGs, $0.6 < z < 1.4$')
	plt.plot(dl[0],dq[1]*dl[0]**2.,color='purple',label=r'quasars, $0.8 < z < 2.2$')
	plt.plot(dl[0],db[1]*dl[0]**2.,color='brown',label=r'BGS, $0.1 < z < 0.4$')
	plt.plot(dl[0],dlt[1]*dl[0]**2.,'--r',label='no fiber assignment')
	plt.plot(dl[0],det[1]*dl[0]**2.,'--b')
	plt.plot(dl[0],dqt[1]*dl[0]**2.,'--',color='purple')
	plt.plot(dl[0],dbt[1]*dl[0]**2.,'--',color='brown')
	plt.legend()
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)')
	plt.ylabel(r'$s^2\xi_0$')
	plt.title('e2e simulation, fiber weight correction')
	plt.ylim(-100,120)
	plt.savefig(dirxi+'xi0e2egcomptrue.png')
	plt.show()

def plotxiBGS_comptrue(zmin=0.1,zmax=0.4):
	fl = 'fkp'+str(zmin)+str(zmax)+'5st0.dat'
	fb = dirxi+'xi024onepergBGS'+fl
	fbtt = dirxi+'xi024onepergtestBGS'+fl
	fbt = dirxi+'xi024onepergBGSztrue'+fl
	db = np.loadtxt(fb).transpose()
	dbt = np.loadtxt(fbt).transpose()
	dbtt = np.loadtxt(fbtt).transpose()
	plt.plot(db[0],db[1]/dbt[1],color='brown',label='fiducial, '+str(zmin)+r'$ < z < $'+str(zmax))
	plt.plot(db[0],dbt[1]/dbt[1],'--',color='brown',label='no fiber assignment')
	plt.plot(db[0],dbtt[1]/dbt[1],':',color='brown',label='test')

	plt.legend()
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)')
	plt.ylabel(r'$\xi_0/\xi_{0,{\rm no-fiber-assignment}}$')
	plt.title('BGS e2e simulation, fiber weight correction')
	plt.ylim(.8,1.2)
	plt.xlim(0,40)
	plt.savefig(dirxi+'xi0e2egBGScomptrue.png')
	plt.show()

def plotxiELG_comptrue(zmin=0.8,zmax=1.6):
	fl = 'fkp'+str(zmin)+str(zmax)
	bs = '5st0.dat'
	fb = dirxi+'xi024e2e_onepergELG'+fl+'0'+bs
	fb1 = dirxi+'xi024e2e_onepergELG'+fl+'00'+bs
	#fbtt = dirxi+'xi024onepergtestBGS'+fl
	fbt = dirxi+'xi024e2e_onepergELGztrue'+fl+'0'+bs
	db = np.loadtxt(fb).transpose()
	db1 = np.loadtxt(fb1).transpose()
	dbt = np.loadtxt(fbt).transpose()
	#dbtt = np.loadtxt(fbtt).transpose()
	plt.plot(db[0],db[1]-dbt[1],color='b',label='fiducial, '+str(zmin)+r'$ < z < $'+str(zmax))
	plt.plot(db[0],db1[1]-dbt[1],':',color='b',label='in dark area, '+str(zmin)+r'$ < z < $'+str(zmax))
	plt.plot(db[0],dbt[1]-dbt[1],'--',color='b',label='no fiber assignment')
	#plt.plot(db[0],dbtt[1]/dbt[1],':',color='b',label='test')

	plt.legend()
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)')
	plt.ylabel(r'$\xi_0-\xi_{0,{\rm no-fiber-assignment}}$')
	plt.title('ELG e2e simulation, fiber weight correction')
	plt.ylim(-0.02,0.01)
	plt.xlim(0,100)
	plt.tight_layout()
	plt.savefig(dirxi+'xi0e2egELGcomptruentile.png')
	plt.show()
	
	plt.clf()

	plt.plot(db[0],db[2]-dbt[2],color='b',label='fiducial, '+str(zmin)+r'$ < z < $'+str(zmax))
	plt.plot(db[0],db1[2]-dbt[2],':',color='b',label='in dark area, '+str(zmin)+r'$ < z < $'+str(zmax))
	plt.plot(db[0],dbt[2]-dbt[2],'--',color='b',label='no fiber assignment')
	#plt.plot(db[0],dbtt[1]/dbt[1],':',color='b',label='test')

	plt.legend()
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)')
	plt.ylabel(r'$\xi_2-\xi_{2,{\rm no-fiber-assignment}}$')
	plt.title('ELG e2e simulation, fiber weight correction')
	plt.ylim(-0.05,0.01)
	plt.xlim(0,100)
	plt.tight_layout()
	plt.savefig(dirxi+'xi2e2egELGcomptruentile.png')
	plt.show()


def plotxiLRG_comptrue(zmin=0.8,zmax=1.6):
	fl = 'fkp'+str(zmin)+str(zmax)+'5st0.dat'
	fb = dirxi+'xi024onepergLRG'+fl
	#fbtt = dirxi+'xi024onepergtestBGS'+fl
	fbt = dirxi+'xi024onepergLRGztrue'+fl
	db = np.loadtxt(fb).transpose()
	dbt = np.loadtxt(fbt).transpose()
	#dbtt = np.loadtxt(fbtt).transpose()
	plt.plot(db[0],db[1]-dbt[1],color='r',label='fiducial, '+str(zmin)+r'$ < z < $'+str(zmax))
	plt.plot(db[0],dbt[1]-dbt[1],'--',color='r',label='no fiber assignment')
	#plt.plot(db[0],dbtt[1]/dbt[1],':',color='b',label='test')

	plt.legend()
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)')
	plt.ylabel(r'$\xi_0-\xi_{0,{\rm no-fiber-assignment}}$')
	plt.title('LRG e2e simulation, fiber weight correction')
	plt.ylim(-0.02,0.01)
	plt.xlim(0,100)
	plt.tight_layout()
	plt.savefig(dirxi+'xi0e2egLRGcomptrue.png')
	plt.show()
	
	plt.clf()
	
	plt.plot(db[0],db[2]-dbt[2],color='r',label='fiducial, '+str(zmin)+r'$ < z < $'+str(zmax))
	plt.plot(db[0],dbt[2]-dbt[2],'--',color='r',label='no fiber assignment')
	#plt.plot(db[0],dbtt[1]/dbt[1],':',color='b',label='test')

	plt.legend()
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)')
	plt.ylabel(r'$\xi_2-\xi_{2,{\rm no-fiber-assignment}}$')
	plt.title('LRG e2e simulation, fiber weight correction')
	plt.ylim(-0.02,0.01)
	plt.xlim(0,100)
	plt.tight_layout()
	plt.savefig(dirxi+'xi2e2egLRGcomptrue.png')
	plt.show()
	


def plotxi2_comptrue():
	fl = dirxi+'xi024onepergLRGfkp0.51.15st0.dat'
	flt = dirxi+'xi024onepergLRGztruefkp0.51.15st0.dat'
	fet = dirxi+'xi024onepergELGztruefkp0.61.45st0.dat'
	fe = dirxi+'xi024onepergELGfkp0.61.45st0.dat'
	fq = dirxi+'xi024onepergQSOfkp0.82.25st0.dat'
	fqt = dirxi+'xi024onepergQSOztruefkp0.82.25st0.dat'
	fb = dirxi+'xi024onepergBGSfkp0.10.45st0.dat'
	fbt = dirxi+'xi024onepergBGSztruefkp0.10.45st0.dat'
	dl = np.loadtxt(fl).transpose() 
	de = np.loadtxt(fe).transpose()
	dq = np.loadtxt(fq).transpose()
	db = np.loadtxt(fb).transpose()
	dlt = np.loadtxt(flt).transpose() 
	det = np.loadtxt(fet).transpose()
	dqt = np.loadtxt(fqt).transpose()
	dbt = np.loadtxt(fbt).transpose()
	plt.plot(dl[0],dl[2]*dl[0],color='r',label=r'LRGs, $0.5 < z < 1.1$')
	plt.plot(dl[0],de[2]*dl[0],color='b',label=r'ELGs, $0.6 < z < 1.4$')
	plt.plot(dl[0],dq[2]*dl[0],color='purple',label=r'quasars, $0.8 < z < 2.2$')
	plt.plot(dl[0],db[2]*dl[0],color='brown',label=r'BGS, $0.1 < z < 0.4$')
	plt.plot(dl[0],dlt[2]*dl[0],'--r',label='no fiber assignment')
	plt.plot(dl[0],det[2]*dl[0],'--b')
	plt.plot(dl[0],dqt[2]*dl[0],'--',color='purple')
	plt.plot(dl[0],dbt[2]*dl[0],'--',color='brown')
	plt.legend()
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)')
	plt.ylabel(r'$s\xi_2$')
	plt.title('e2e simulation, fiber weight correction')
	plt.ylim(-8,4)
	plt.savefig(dirxi+'xi2e2egcomptrue.png')
	plt.show()

def plotxi4_comptrue():
	fl = dirxi+'xi024onepergLRGfkp0.51.15st0.dat'
	flt = dirxi+'xi024onepergLRGztruefkp0.51.15st0.dat'
	fet = dirxi+'xi024onepergELGztruefkp0.61.45st0.dat'
	fe = dirxi+'xi024onepergELGfkp0.61.45st0.dat'
	fq = dirxi+'xi024onepergQSOfkp0.82.25st0.dat'
	fqt = dirxi+'xi024onepergQSOztruefkp0.82.25st0.dat'
	fb = dirxi+'xi024onepergBGSfkp0.10.45st0.dat'
	fbt = dirxi+'xi024onepergBGSztruefkp0.10.45st0.dat'
	dl = np.loadtxt(fl).transpose() 
	de = np.loadtxt(fe).transpose()
	dq = np.loadtxt(fq).transpose()
	db = np.loadtxt(fb).transpose()
	dlt = np.loadtxt(flt).transpose() 
	det = np.loadtxt(fet).transpose()
	dqt = np.loadtxt(fqt).transpose()
	dbt = np.loadtxt(fbt).transpose()
	plt.plot(dl[0],dl[3]*dl[0],color='r',label=r'LRGs, $0.5 < z < 1.1$')
	plt.plot(dl[0],de[3]*dl[0],color='b',label=r'ELGs, $0.6 < z < 1.4$')
	plt.plot(dl[0],dq[3]*dl[0],color='purple',label=r'quasars, $0.8 < z < 2.2$')
	plt.plot(dl[0],db[3]*dl[0],color='brown',label=r'BGS, $0.1 < z < 0.4$')
	plt.plot(dl[0],dlt[3]*dl[0],'--r',label='no fiber assignment')
	plt.plot(dl[0],det[3]*dl[0],'--b')
	plt.plot(dl[0],dqt[3]*dl[0],'--',color='purple')
	plt.plot(dl[0],dbt[3]*dl[0],'--',color='brown')
	plt.legend()
	plt.xlabel(r'$s$ ($h^{-1}$Mpc)')
	plt.ylabel(r'$s\xi_4$')
	plt.title('e2e simulation, fiber weight correction')
	plt.ylim(-4,4)
	plt.savefig(dirxi+'xi4e2egcomptrue.png')
	plt.show()



if __name__ == '__main__':
    #export OMP_NUM_THREADS=64
    #module load gsl
    import subprocess
    truez = ''
    type = 'LRG'
    zmin=0.5
    zmax=1.1
    prep4czxi(type,0.5,1.1,truez=truez)
    subprocess.run(['chmod','+x','czpc.sh'])
    subprocess.run('./czpc.sh')
    calcxi_dataCZ(type,0.5,1.1,truez=truez)
    #plotxi_compfkp(type,0.5,1.1)
# 
#	plotxiLRG_comptrue(zmin,zmax)

    type = 'ELG'
    test = 'test'
    zmin=0.6
    zmax=1.4
	#prep4czxi(type,0.6,1.4,program='gray',truez=truez,test=test,pg=-1)
	#subprocess.run(['chmod','+x','czpc.sh'])
	#subprocess.run('./czpc.sh')
#   calcxi_dataCZ(type,0.6,1.4,truez=truez,pg=-1)
# 	prep4czxi(type,0.6,1.4,program='gray',truez=truez,mintile=1,test=test)
# 	subprocess.run(['chmod','+x','czpc.sh'])
# 	subprocess.run('./czpc.sh')
# 	calcxi_dataCZ(type,0.6,1.4,truez=truez,mintile=1)
# 	truez = 'ztrue'
# 	prep4czxi(type,0.6,1.4,program='gray',truez=truez,test=test)
# 	subprocess.run(['chmod','+x','czpc.sh'])
# 	subprocess.run('./czpc.sh')
# 	calcxi_dataCZ(type,0.6,1.4,truez=truez)

	#plotxi_compfkp(type,0.6,1.4)
	#plotxiELG_comptrue(zmin,zmax)

# 	type = 'QSO'
# 	prep4czxi(type,0.8,2.2,truez=truez)
# 	subprocess.run(['chmod','+x','czpc.sh'])
# 	subprocess.run('./czpc.sh')
# 	calcxi_dataCZ(type,0.8,2.2,truez=truez)
# 	plotxi_compfkp(type,0.8,2.2)
# 
	#type = 'BGS'
	
	#zmin =0.1
	#zmax = 0.4
	#prep4czxi(type,zmin,zmax,program='bright',truez=truez,test=test)
	#subprocess.run(['chmod','+x','czpc.sh'])
	#subprocess.run('./czpc.sh')
	#calcxi_dataCZ(type,zmin,zmax,truez=truez,test=test)
	#plotxi_compfkp(type,zmin,zmax)
	#plotxiBGS_comptrue(zmin,zmax)

	#plotxi_comptrue()
	#plotxi2_comptrue()
	#plotxi4_comptrue()


# 	ppxilcalc_LSDfjack_bs(type,tile,night,zmin=.5,zmax=1.1)
# 	ppxilcalc_LSDfjack_bs(type,tile,night,zmin=.5,zmax=1.1,bs=5)

