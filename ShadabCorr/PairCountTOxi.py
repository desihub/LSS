#Author: Shadab Alam, March 2015
#This program takes the pair count and normalization from 
#PairCount's code output and compute 2d correlation and
#monopole and quadrupole moment of corrlation

import numpy as np
import os
import sys

#import pylab as pl
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as pl

f_write=[]
f_plot=[]

def load_PairCount_output(root):
   #load sbins
   lines=open(root+'-DD.dat').readlines()
   sbins=np.array([np.float(x) for x in lines[0].split()])
   ns=sbins.size-1
   #load mu bins
   mubins=np.array([np.float(x) for x in lines[1].split()])
   nmu=mubins.size-1
   #load pair counts
   DD=np.loadtxt(root+'-DD.dat',skiprows=2)
   DR=np.loadtxt(root+'-DR.dat',skiprows=2)
   try:
      RR=np.loadtxt(root+'-RR.dat',skiprows=2)
   except:
      RR=np.loadtxt(root[:-4]+'0001-RR.dat',skiprows=2)
   #load wieghts
   lines=open(root+'-norm.dat').readlines()
   Wsum_data=np.float(lines[0].split()[-1])
   Wsum_rand=np.float(lines[1].split()[-1])

   Wrat=Wsum_data/Wsum_rand
   DR=DR*Wrat
   RR=RR*Wrat*Wrat
   #print ns,nmu,DD.shape, Wsum_data, Wsum_rand
   #print sbins, mubins

   #combined all output in a dictionary and return
   pcdict={'sbins': sbins, 'ns': ns, 'mubins':mubins, 'nmu':nmu, 
	 'DD': DD, 'DR': DR, 'RR': RR ,
	 'Wsum_data': Wsum_data, 'Wsum_rand': Wsum_rand, 'Wrat':Wrat}
   return pcdict

def load_PairCount_output_cross(root):
   #load sbins
   lines=open(root+'-D1D2.dat').readlines()
   sbins=np.array([np.float(x) for x in lines[0].split()])
   ns=sbins.size-1
   #load mu bins
   mubins=np.array([np.float(x) for x in lines[1].split()])
   nmu=mubins.size-1
   #load pair counts
   D1D2=np.loadtxt(root+'-D1D2.dat',skiprows=2)
   D1R2=np.loadtxt(root+'-D1R2.dat',skiprows=2)
   R1D2=np.loadtxt(root+'-R1D2.dat',skiprows=2)
   R1R2=np.loadtxt(root+'-R1R2.dat',skiprows=2)

   #load wieghts
   lines=open(root+'-norm.dat').readlines()
   Wsum_data1=np.float(lines[0].split()[-1])
   Wsum_rand1=np.float(lines[1].split()[-1])
   Wsum_data2=np.float(lines[2].split()[-1])
   Wsum_rand2=np.float(lines[3].split()[-1])

  
   Wrat1=Wsum_data1/Wsum_rand1
   Wrat2=Wsum_data2/Wsum_rand2

   #Normalize the pair counts
   D1R2=D1R2*Wrat2
   R1D2=R1D2*Wrat1
   R1R2=R1R2*Wrat1*Wrat2

   #combined all output in a dictionary and return
   pcdict={'sbins': sbins, 'ns': ns, 'mubins':mubins, 'nmu':nmu, 
	 'D1D2': D1D2, 'D1R2': D1R2, 'R1D2': R1D2,'R1R2': R1R2 }

   return pcdict





def combine_count(root1,root2):
  #load the count for each of the root
  sbins1, ns1, mubins1, nmu1, DD1,DR1,RR1, Wsum_data1,Wsum_rand1=load_PairCount_output(root1)
  sbins2, ns2, mubins2, nmu2, DD2,DR2,RR2, Wsum_data2,Wsum_rand2=load_PairCount_output(root2)

  Wrat=(Wsum_data1+Wsum_data2)/(Wsum_rand1+Wsum_rand2)
  #giving normalized pair counts
  DD=DD1+DD2
  DR=(DR1+DR2)*Wrat
  RR=(RR1+ RR2)*Wrat*Wrat

  if(np.sum(sbins1==sbins2)==sbins1.size and ns1==ns2 
     and np.sum(mubins1==mubins2)==mubins1.size and nmu1==nmu2):
     return sbins1, ns1, mubins1, nmu1, DD,DR,RR
  else:
     return 'error'

def compute_xi02(root,pcdict,xitype='auto', nscomb=1,write=0):
   if(write):
      outfile=root+'-xi2D.dat'
      f_write.append(outfile)
      fout=open(outfile,'w')

   if(xitype=='auto'):
      pclist=['DD', 'DR','RR']
   elif(xitype=='cross'):
      pclist=['D1D2', 'D1R2','R1D2','R1R2']

   #for rebinning
   nsrebins=np.int(np.ceil(pcdict['ns']/np.float(nscomb)))
   srebins=np.zeros(nsrebins)
   nmurebins=pcdict['nmu']
   murebins=np.zeros(nmurebins)

   xi2drebin=np.zeros(nsrebins*nmurebins).reshape(nsrebins,nmurebins)

   rebin=-1
   for ii in range(0,pcdict['ns'],nscomb):
      ss=np.mean(pcdict['sbins'][ii:ii+nscomb+1])
      rebin=rebin+1
      srebins[rebin]=ss
      for jj in range(0,nmurebins):
         mu=np.mean(pcdict['mubins'][jj:jj+2])
         if(rebin==0):
            murebins[jj]=mu

	 pcbin={} 
	 for pctype in pclist:
	    pcbin[pctype]=np.sum(pcdict[pctype][ii:ii+nscomb,jj])

         if(xitype=='auto'):
            xi2d=(pcbin['DD']-2*pcbin['DR']+pcbin['RR'])/pcbin['RR']
	 elif(xitype=='cross'):
            xi2d=(pcbin['D1D2']-pcbin['D1R2']-pcbin['R1D2']+pcbin['R1R2'])/pcbin['R1R2']

         xi2drebin[rebin,jj]=xi2d
         if(write):
            fout.write('%12.8lf %12.8lf %12.8lf' %(ss,mu,xi2d))
            for pctype in pclist:
	       fout.write('%12.8lf '%pcbin[pctype])
	    fout.write('\n')

   if(write):
      fout.close()

   #print srebins, murebins
   return srebins, murebins, xi2drebin

#Assume sampling in rp-rpi computes xi2D and wp
def compute_xi2D_rppi(xi2droot,wproot,pcdict,xitype='auto',nscomb=1,write=0):
   if(write):
      outfile=xi2droot+'-xi2D.dat'
      f_write.append(outfile)
      fout=open(outfile,'w')

      outfile=wproot+'-wp.dat'
      f_write.append(outfile)
      fwp=open(outfile,'w')
      

   if(xitype=='auto'):
      pclist=['DD', 'DR','RR']
   elif(xitype=='cross'):
      pclist=['D1D2', 'D1R2','R1D2','R1R2']

   wp=np.zeros(pcdict['nper'])
   #for rebinning
   nper_rebins=pcdict['nper']
   npar_rebins=np.int(np.ceil(pcdict['npar']/np.float(nscomb)))
   rper_rebins=np.zeros(nper_rebins)
   rpar_rebins=np.zeros(npar_rebins)

   xi2drebin=np.zeros(nper_rebins*npar_rebins).reshape(nper_rebins,npar_rebins)

   for ii in range(0,pcdict['nper'],1):
      ss=np.mean(pcdict['rper'][ii:ii+2])
      rper_rebins[ii]=ss
      rebin=-1
      for jj in range(0,pcdict['npar'],nscomb):
         mu=np.mean(pcdict['rpar'][jj:jj+nscomb+1])
         rebin=rebin+1
         rpar_rebins[rebin]=mu

	 pcbin={} 
	 for pctype in pclist:
	    pcbin[pctype]=np.sum(pcdict[pctype][ii,jj:jj+nscomb])

         if(xitype=='auto'):
            xi2d=(pcbin['DD']-2*pcbin['DR']+pcbin['RR'])/pcbin['RR']
	 elif(xitype=='cross'):
            xi2d=(pcbin['D1D2']-pcbin['D1R2']-pcbin['R1D2']+pcbin['R1R2'])/pcbin['R1R2']

         xi2drebin[ii,rebin]=xi2d
         if(write):
            fout.write('%12.8lf %12.8lf %12.8lf' %(ss,mu,xi2d))
            for pctype in pclist:
	       fout.write('%12.8lf '%pcbin[pctype])
	    fout.write('\n')

      dr_par=np.mean(rpar_rebins[1:]-rpar_rebins[:-1])
      wp[ii]=np.sum(xi2drebin[ii,:])*dr_par
      if(write):
         fwp.write('%12.8lf %12.8lf\n'%(ss,wp[ii]))

   if(write):
      fout.close()
      fwp.close()

   return rper_rebins, rpar_rebins, xi2drebin,wp
#end of compute wp and xi2d with rprpi sampling

#a function to call for pair count to xi2d and xi02 in polar co-ordinate
def xi2d_xi02_froot(root,xi2droot,xi02root,xitype='auto',nscomb=1,samp='rmu',write=0):
   #compute the correlation functions
   if(xitype=='auto'):
      pcdict=load_PairCount_output(root)
   elif(xitype=='cross'):
      pcdict=load_PairCount_output_cross(root)

   if(samp=='logr-theta'):
      pcdict['sbins']=np.power(10,pcdict['sbins'])

   srebins, murebins, xi2drebin=compute_xi02(xi2droot,pcdict,xitype=xitype,nscomb=nscomb,write=write)
   xi02=Xi_Legendre(xi02root,srebins, murebins, xi2drebin,samp=samp,write=write)

   return srebins,murebins,xi2drebin,xi02

#a function to call for pair count to xi2d and wp in cartesian co-ordinate
def xi2d_wp_froot(root,xi2droot,wproot,xitype='auto',nscomb=1,samp='log',write=0):
   #compute the correlation functions
   #rper,nper,rpar,npar, DD,DR,RR, Wsum_data,Wsum_rand=load_PairCount_output(root)
   if(xitype=='auto'):
      pcdict=load_PairCount_output(root)
   elif(xitype=='cross'):
      pcdict=load_PairCount_output_cross(root)

   #change dictionary key for cartesian coordinate
   keyreplace={'rper': 'sbins','nper': 'ns','rpar':'mubins','npar':'nmu'}
   for key in keyreplace.keys():
      pcdict[key]=pcdict.pop(keyreplace[key])

   if(samp=='logrp-pi'):
      pcdict['rper']=np.power(10,pcdict['rper'])

   rper_rebins,rpar_rebins,xi2drebin,wp=compute_xi2D_rppi(xi2droot,wproot,
                   pcdict,xitype=xitype,nscomb=nscomb,write=write)
   
   return rper_rebins, rpar_rebins, xi2drebin,wp

def xi2d_wp_xi02_froot_jn(inroot,xi2droot,outroot,xitype='auto',nscomb=1,samp='logrp-pi',compgZ=0,NJN=1):
   polar_coord=['rmu','rtheta','logr-theta']
   cart_coord=['rp-pi','logrp-pi']
   #files to write the compiled xi2d and wp functions from jackife
   xi2dfile=xi2droot+'-'+samp+'-NJN-%d.txt'%NJN

   for ii in range(0,NJN+1):
      if(NJN==0):
         root=inroot+'-'+samp
      elif(ii<NJN): #individual regions
         root=inroot+'-'+samp+'-JNdir/jk-%d'%ii
      else: #All data
         root=inroot+'-All-'+samp
      if(samp in polar_coord):
         srebins,murebins,xi2drebin,xi02=xi2d_xi02_froot(root,
	                    '','',xitype=xitype,samp=samp,nscomb=nscomb,write=0)
	 if(compgZ==1):
	    gZ=compute_gZ('',srebins, murebins, xi2drebin,samp=samp,write=0)

      elif(samp in cart_coord):
         rper_rebins, rpar_rebins, xi2drebin,wp=xi2d_wp_froot(root,
                           '','',xitype=xitype,samp=samp,nscomb=nscomb,write=0)
      else:
	 print 'Invalid sampling: %s'%samp
	 sys.exit()

      if(ii==0):
         n2d=xi2drebin.size
         xi2dmat=np.zeros(n2d*(NJN+1)).reshape(n2d,NJN+1)
         if(samp in polar_coord):
	    nxi02=xi02.shape[0]
	    rxi02=xi02[:,0]
	    xi0mat=np.zeros(nxi02*(NJN+1)).reshape(nxi02,NJN+1)
	    xi2mat=np.zeros(nxi02*(NJN+1)).reshape(nxi02,NJN+1)
	    xi1mat=np.zeros(nxi02*(NJN+1)).reshape(nxi02,NJN+1) #dipole usefule only for rtheta
	    if(compgZ==1):
	       gZmat=np.zeros(nxi02*(NJN+1)).reshape(nxi02,NJN+1)
	 else:
            nwp=rper_rebins.size
            wpmat = np.zeros(nwp*(NJN+1)).reshape(nwp,NJN+1)

      xi2dmat[:,ii]=xi2drebin.reshape(n2d)
      if(samp in polar_coord):
	 xi0mat[:,ii]=xi02[:,1]
	 xi2mat[:,ii]=xi02[:,2]
	 xi1mat[:,ii]=xi02[:,3]
	 if(compgZ==1):
	    gZmat[:,ii]=gZ[:,1]
      else:
         wpmat[:,ii]=wp


   if(samp in polar_coord):
      xi2daxes=[srebins,murebins]
      col1=rxi02
      outtags=['xi0','xi2']
      outmats=[xi0mat,xi2mat]
      if(samp=='rmu'):
         xi2dtag=['r','mu']
      elif(samp=='rtheta' or samp=='logr-theta'):
         xi2dtag=['r','theta']
         outmats.append(xi1mat)
	 outtags.append('xi1')
      else:
	 print 'Invalid sampling: %s'%samp

      if(compgZ==1):
	 outtags.append('gZ')
         outmats.append(gZmat)

   else:
      xi2daxes=[rper_rebins,rpar_rebins]
      xi2dtag=['rper','rpar']
      outmats=[wpmat]
      col1=rper_rebins
      outtags=['wp']

   #compute the mean and std err for xi2d and write to file
   meanxi2d=np.mean(xi2dmat[:,:NJN],axis=1)
   errxi2d=np.sqrt(NJN-1)*np.std(xi2dmat[:,:NJN],axis=1)
   xi2dwrite=np.column_stack([meanxi2d,errxi2d,
                                   xi2dmat[:,NJN],xi2dmat[:,:NJN]])
   with file(xi2dfile,'w') as fxi2d:
      fxi2d.write('#xi2d with jacknife, root:%s NJN=%d\n'%(root,NJN))
      for ww,arr in enumerate(xi2daxes):
         w_str='#%s: '%xi2dtag[ww]
         for rr in arr:
            w_str="%s %10.5f"%(w_str,rr)
         fxi2d.write(w_str+'\n')

      if(NJN>0):
         fxi2d.write('#meanxi2d sigmaxi2d Allxi2d jacknifecolumns\n')
         np.savetxt(fxi2d,xi2dwrite,fmt='% 18.8lf')
      elif(NJN==0):
         fxi2d.write('#xi2d\n')
	 np.savetxt(fxi2d,xi2dmat[:,NJN],fmt='% 18.8lf')

   print 'written: ',xi2dfile


   for pole, outmat in enumerate(outmats):      
     #compute the mean and std err for wp and write to file
      meanout=np.mean(outmat[:,:NJN],axis=1)
      errout=np.sqrt(NJN-1)*np.std(outmat[:,:NJN],axis=1)
      outwrite=np.column_stack([col1,meanout,errout,outmat[:,NJN],outmat[:,:NJN]])
      outfile=outroot+'-'+outtags[pole]+'-'+samp+'-NJN-%d.txt'%NJN
      with file(outfile,'w') as fout:
         fout.write('#%s with jacknife, root:%s NJN=%d\n'%(outtags[pole],root,NJN))
	 if(NJN>0):
            fout.write('#r(Mpc/h) mean sigma All jacknifecolumns\n')
            np.savetxt(fout,outwrite,fmt='% 18.8lf')
	 elif(NJN==0):
            fout.write('#r(Mpc/h) val\n')
	    np.savetxt(fout,np.column_stack([col1,outmat[:,NJN]]),fmt='% 18.8lf')
      print 'written:',outfile


   return 0



def compute_xi02_cross(root,sbins, ns, mubins, nmu, 
        D1D2,D1R2,R1D2,R1R2, nscomb=1,write=0):
   if(write):
      outfile=root+'cross-xi2D.dat'
      f_write.append(outfile)
      fout=open(outfile,'w')

   #for rebinning
   nsrebins=np.int(np.ceil(ns/np.float(nscomb)))
   srebins=np.zeros(nsrebins)
   nmurebins=nmu
   murebins=np.zeros(nmu)

   xi2drebin=np.zeros(nsrebins*nmu).reshape(nsrebins,nmu)

   rebin=-1
   for ii in range(0,ns,nscomb):
      ss=np.mean(sbins[ii:ii+nscomb+1])
      rebin=rebin+1
      srebins[rebin]=ss
      for jj in range(0,nmu):
         mu=np.mean(mubins[jj:jj+2])
         if(rebin==0):
            murebins[jj]=mu
         D1D2bin=np.sum(D1D2[ii:ii+nscomb,jj])
         D1R2bin=np.sum(D1R2[ii:ii+nscomb,jj])
         R1D2bin=np.sum(R1D2[ii:ii+nscomb,jj])
         R1R2bin=np.sum(R1R2[ii:ii+nscomb,jj])

         xi2d=(D1D2bin-D1R2bin-R1D2bin+R1R2bin)/R1R2bin
         xi2drebin[rebin,jj]=xi2d
         if(write):
            fout.write('%12.8lf %12.8lf %12.8lf\n' %(ss,mu,xi2d))


   if(write):
      fout.close()

   #print srebins, murebins
   return srebins, murebins, xi2drebin
#End of compute xi02 cross correlation

def Xi_Legendre(root,srebins, murebins, xi2drebin,samp='rmu',write=1):
   if(samp=='rmu'):
      dmu=np.average(murebins[1:]-murebins[:-1])
      diffmu=1 #to account for the fact that integral is 0 to 1
   elif(samp=='rtheta' or samp=='logr-theta'): #in this case murebin is actually theta and not mu
      dtheta=np.average(murebins[1:]-murebins[:-1])
      dmu=np.sin(murebins)*dtheta
      murebins=np.cos(murebins) #this converts theta to mu for integration
      diffmu=2 #to account for the fact that itegral is -1 to 1

   mu2=np.power(murebins,2)
   P0 =np.ones(mu2.size)/2.0    #(2l+1)Pl(mu)*sqrt(1-mu^2) monopole term
   P2=2.5*(3*mu2-1)/2.0         #p2=(3mu^2-1)/2  quadrupole term
   P1=3*murebins/2.0
   

   ns=srebins.size 
   xi02=np.zeros(ns*4).reshape(ns,4)

   if(write):
      outfile=root+'-xi02.dat'
      f_write.append(outfile)
      fout=open(outfile,'w')

   for ii in range(0,ns):
      xi02[ii,0]=srebins[ii]
      xi02[ii,1]=np.sum(xi2drebin[ii,:]*P0*dmu)*2.0/diffmu
      xi02[ii,2]=np.sum(xi2drebin[ii,:]*P2*dmu)*2.0/diffmu
      xi02[ii,3]=np.sum(xi2drebin[ii,:]*P1*dmu)*2.0/diffmu
      if(write==1 and smap=='rmu'):
         fout.write('%12.8lf %12.8lf %12.8lf\n' %(xi02[ii,0],xi02[ii,1],xi02[ii,2]))
      if(write==1 and smap=='rtheta'):
         fout.write('%12.8lf %12.8lf %12.8lf %12.8lf\n' %(xi02[ii,0],xi02[ii,1],xi02[ii,2],xi02[ii,3]))

   
   #plot_xi02(xi02)
   return xi02


def compute_gZ(root,srebins, murebins, xi2drebin,samp='rtheta',write=1):
   '''This function computes the gravitational redshift from 2d correlation'''

   H=100 #hubble constant

   if(samp=='rtheta' or samp=='logr-theta'): #in this case murebin is actually theta and not mu
      dtheta=np.average(murebins[1:]-murebins[:-1])
      dmu=np.sin(murebins)*dtheta
      murebins=-np.cos(murebins) #this converts theta to mu for integration
      diffmu=2 #to account for the fact that itegral is -1 to 1
   else:
      print 'invalid sampling for gravitaional redshift: %s'%samp
      return 0

   ns=srebins.size 
   gZ=np.zeros(ns*2).reshape(ns,2)
   numer=np.zeros(ns)
   denom=np.zeros(ns)

   if(write):
      outfile=root+'-gz.dat'
      f_write.append(outfile)
      fout=open(outfile,'w')

   #print np.min(murebins),np.max(murebins)
   for ii in range(0,ns):
      gZ[ii,0]=srebins[ii]
      #applying a mu limit
      mulim=1.0 #1.0
      ind1=murebins<mulim; ind2=murebins>-mulim; ind12=ind1*ind2
      #remove infinities and Nan
      indNaN=np.isnan(xi2drebin[ii,:])
      indinf=np.isinf(xi2drebin[ii,:])
      indcheck=indNaN+indinf
      #remove the infinity and its symmetric points
      indsym=np.copy(indcheck)
      
      for muinf in murebins[indcheck]:
         indsym=indsym+(np.abs(murebins+muinf)<1e-4)
      indcheck=indsym
      #print indcheck
      indall=~indcheck*ind12
      if(np.sum(indcheck)>0):
	 print 'infinities:(r,mu,xi)',srebins[ii],murebins[indcheck],xi2drebin[ii,indcheck]
	 #print ~indcheck*ind12
	 #import pylab as pl
	 #print indcheck
	 #pl.plot(murebins[~indcheck],xi2drebin[ii,~indcheck],'s-')
	 #pl.show()

      #summing
      #denom[ii]=np.sum((1+xi2drebin[ii,indall])*gZ[ii,0]*dtheta)
      #numer[ii]=H*np.sum((1+xi2drebin[ii,indall])*gZ[ii,0]*gZ[ii,0]*murebins[indall]*dtheta)
      #simpson rule for integration
      #Interpolation integration library
      import scipy.integrate as ssI
      from scipy import interpolate
      denom[ii]=ssI.simps(1+xi2drebin[ii,indall],murebins[indall])
      numer[ii]=H*ssI.simps((1+xi2drebin[ii,indall])*gZ[ii,0]*murebins[indall],murebins[indall])

      gZ[ii,1]=numer[ii]/denom[ii]

      if(write==1):
         fout.write('%12.8lf %12.8lf\n' %(gZ[ii,0],gZ[ii,1]))

   if(write==1):
      fout.close()

   if(0):#To integrate along r
      Rmin=1.0
      nR=15
      Rmax=60.0
      dLR=(np.log(Rmax)-np.log(Rmin))/nR
      #interpolate numer and denomr
      Idenomr = interpolate.splrep(srebins,denom, s=0,k=1)
      Inumerr = interpolate.splrep(srebins,numer, s=0,k=1)
      rr2d_tmp=np.linspace(0,Rmax,2000)
      denom_rtmp=interpolate.splev(rr2d_tmp,Idenomr, der=0)
      numer_rtmp=interpolate.splev(rr2d_tmp,Inumerr, der=0)

      rr=np.zeros(nR)
      gZnew=np.zeros(ns*2).reshape(ns,2)
      for ii in range(0,nR):
         r1=np.exp(np.log(Rmin)+ii*dLR)
         r2=np.exp(np.log(Rmin)+(ii+1)*dLR)
         rr[ii]=0.5*(r1+r2)
         ind1=rr2d_tmp>r1
         ind2=rr2d_tmp<r2
         denominator =ssI.simps(denom_rtmp[ind1*ind2],rr2d_tmp[ind1*ind2])
         numerator   =ssI.simps(numer_rtmp[ind1*ind2],rr2d_tmp[ind1*ind2])

         gZnew[ii,0]=rr[ii]
         gZnew[ii,1]=numerator/denominator
         if(0):
            import pylab as pl
            pl.figure(10)
            pl.plot(gZnew[ii,0],denominator,'b*')
            pl.plot(rr2d_tmp[ind1*ind2],denom_rtmp[ind1*ind2],'b-')
            pl.figure(20)
            pl.plot(gZnew[ii,0],numerator,'r*')
            pl.plot(rr2d_tmp[ind1*ind2],numer_rtmp[ind1*ind2],'r-')
            pl.figure(30)
            pl.plot(gZnew[ii,0],gZnew[ii,1],'k*')
            pl.plot(rr2d_tmp[ind1*ind2],numer_rtmp[ind1*ind2]/denom_rtmp[ind1*ind2],'k-')
            pl.show()
      #return np.column_stack([gZ,gZnew])
      #pl.plot(gZnew[:,0],gZnew[:,1])
      #pl.show()
      return gZnew

   return gZ
#end compute_gZ

def plot_xi02(xi02,mark='o-',lab='',lw=1):
   r=xi02[:,0]
   r2=r*r
   pl.plot(r,r2*xi02[:,1],mark,linewidth=lw,label=lab)
   pl.plot(r,r2*xi02[:,2],mark,linewidth=lw)
   pl.xlabel('r (Mpc/h)',fontsize=32)
   pl.ylabel(r'$r^2 \xi_l$',fontsize=32)
   if(lab!=''):
      pl.legend()

   return 0

def mean_cov(froot='',fold='',sky='NS',Njack=28,Sjack=9):
   xiarr=0
   reg=0
   if(sky=='NS'):
      jack=Njack+Sjack
   elif(sky=='N'):
      jack=Njack
   else:
      jack=Sjack

   title=''
   pl.figure()
   for s in sky:
      if s=='N':
         ext='-xi02.dat'
         nloop=Njack
         mark='r-'
         title=title+' N=red '
      elif s=='S':
         ext='-xi02.dat'
         nloop=Sjack
         mark='b-'
         title=title+' S=blue '
      for ii in range(0,nloop):
         if(fold=='../XI02/'):
            fname=froot+str(ii)+'-'+s+ext
         else:
            fname=froot+str(reg)+'-'+s+str(ii)+ext
         print fname
         data=np.loadtxt(fold+fname)
         if(ii==0 and xiarr==0):
            ns=data.shape[0]
            rr=data[:,0]
            rr2=rr*rr
            xi02=np.zeros(ns*2*jack).reshape(2*ns,jack)
            xiarr=1
         xi02[:ns,reg]=data[:,1]
         xi02[ns:,reg]=data[:,2]
         pl.plot(rr,rr2*data[:,1],mark)
         pl.plot(rr,rr2*data[:,2],mark)
         reg=reg+1

   #mean and error bar      
   mean_xi=np.mean(xi02,axis=1)
   cov=np.cov(xi02)*jack
   xierr=np.sqrt(np.diagonal(cov))
   pl.errorbar(rr,rr2*mean_xi[:ns],yerr=rr2*xierr[:ns],fmt='o',color='k',markersize=10,label='mean')
   pl.errorbar(rr,rr2*mean_xi[ns:],yerr=rr2*xierr[ns:],fmt='o',color='k',markersize=10)

   #plot full sample
   data=np.loadtxt(fold+'All-'+sky+'-xi02.dat')
   pl.plot(rr,rr2*data[:,1],'m*-',linewidth=4,markersize=10,label='full sample')
   pl.plot(rr,rr2*data[:,2],'m*-',linewidth=4,markersize=10)
   pl.xlabel(r'$r Mpc/h $',fontsize=32)
   pl.ylabel(r'$r^2 \xi_l$',fontsize=32)
   pl.title(title,fontsize=32)
   pl.legend()
   pl.tight_layout()
   outplot='plots/'+froot+'-'+sky+'.png'
   f_plot.append(outplot)
   pl.savefig(outplot)

   pl.figure()
   pl.pcolor(cov)
   pl.colorbar()

   corr=np.copy(cov)
   for ii in range(0,cov.shape[0]):
      for jj in range(0,cov.shape[1]):
         corr[ii,jj]=cov[ii,jj]/np.sqrt(cov[ii,ii]*cov[jj,jj])


   pl.figure()
   pl.pcolor(corr)
   pl.colorbar()
   outplot='plots/'+froot+'-corr.png'
   f_plot.append(outplot)
   pl.savefig(outplot)
   
 
   return 0

def Mocks_xi_comb(rootN='N',rootS='S',nmock1=0,nmock2=1000,outfold='xi2d',plots_dir='plots',nscomb=1,write=0):

   #compute north only and south only and north+south correlations
   pl.figure(1)  #for xi02 each region
   pl.figure(2)  #for xi02 combined for each region
   for ii in range(nmock1,nmock2):
      if(ii%50==0):
         print ii
      #compute the correlation functions for north
      root1=rootN+'-'+str(ii).zfill(4)
      outroot=outfold+'ngc/XI2D/'+root1.split('/')[-1]
      sbins, ns, mubins, nmu, DD,DR,RR, Wsum_data,Wsum_rand=load_PairCount_output(root1)
      Wrat=Wsum_data/Wsum_rand
      srebins, murebins, xi2drebin=compute_xi02(outroot,sbins, ns, mubins, nmu,
                                  DD,DR,RR, Wrat,nscomb=nscomb,write=write)
      xi02=Xi_Legendre(root1,srebins, murebins, xi2drebin,outfold=outfold+'ngc/XI02/',write=write)
      pl.figure(1)
      plot_xi02(xi02,mark='r-')

      #compute the correlation function for south
      root2=rootS+'-'+str(ii).zfill(4)
      outroot=outfold+'sgc/XI2D/'+root2.split('/')[-1]
      sbins, ns, mubins, nmu, DD,DR,RR, Wsum_data,Wsum_rand=load_PairCount_output(root2)
      Wrat=Wsum_data/Wsum_rand
      srebins, murebins, xi2drebin=compute_xi02(outroot,sbins, ns, mubins, nmu,
                                  DD,DR,RR, Wrat,nscomb=nscomb,write=write)
      xi02=Xi_Legendre(root2,srebins, murebins, xi2drebin,outfold=outfold+'sgc/XI02/',write=write)
      pl.figure(2)
      plot_xi02(xi02,mark='b-')

      #compute the combined correlation function
      sbins, ns, mubins, nmu, DD,DR,RR=combine_count(root1,root2)
      outroot=outfold+'combgc/XI2D/'+root1.split('/')[-1]
      srebins, murebins, xi2drebin=compute_xi02(outroot,sbins, ns, mubins, nmu,
                                  DD,DR,RR, 1.0,nscomb=nscomb,write=write)
      xi02=Xi_Legendre(root1,srebins, murebins, xi2drebin,outfold=outfold+'combgc/XI02/',write=write)
      pl.figure(3)
      plot_xi02(xi02,mark='k-')


   pl.figure(1)
   pl.title(rootN.split('/')[-1])
   pl.savefig(plots_dir+rootN.split('/')[-1]+'-ngc.png')

   pl.figure(2)
   pl.title(rootS.split('/')[-1])
   pl.savefig(plots_dir+rootS.split('/')[-1]+'-sgc.png')

   pl.figure(3)
   pl.title(rootN.split('/')[-1])
   pl.savefig(plots_dir+rootN.split('/')[-1]+'-comb.png')

   return 0
#End of Mocks


def AllJN_realize(nscomb=8):
   rootN='../CMASS-N-JN/'
   rootS='../CMASS-S-JN/'
   outfold1='../XI02/'
   outfoldcomb='../XI02-NS/'

   njack_N=28
   njack_S=9
   reg=0

   #compute north only and north + souht correlations
   root2=rootS+'All-S'
   pl.figure(1)  #for xi02 each region
   pl.figure(2)  #for xi02 combined for each region
   for ii in range(0,njack_N):
      root1=rootN+'JN'+str(ii)+'-N'
      #compute the correlation functions
      sbins, ns, mubins, nmu, DD,DR,RR, Wsum_data,Wsum_rand=load_PairCount_output(root1)
      Wrat=Wsum_data/Wsum_rand
      srebins, murebins, xi2drebin=compute_xi02(root1,sbins, ns, mubins, nmu,
                                  DD,DR,RR, Wrat,nscomb=nscomb,write=1)
      xi02=Xi_Legendre(root1,srebins, murebins, xi2drebin,outfold=outfold1,write=1)
      pl.figure(3)
      plot_xi02(xi02,mark='ro-',lab='N-'+str(ii))
      pl.figure(1)
      plot_xi02(xi02,mark='r-')

      #compute the combined correlation function
      sbins, ns, mubins, nmu, DD,DR,RR=combine_count(root1,root2)
      root=rootN+'COMB'+str(reg)+'-N'+str(ii)
      srebins, murebins, xi2drebin=compute_xi02(root,sbins, ns, mubins, nmu,
                                  DD,DR,RR, 1.0,nscomb=nscomb,write=1)
      xi02=Xi_Legendre(root,srebins, murebins, xi2drebin,outfold=outfoldcomb,write=1)
      pl.figure(3)
      plot_xi02(xi02,mark='bo-',lab='COMB')
      outplot='plots/Reg-N'+str(ii)+'.png'
      f_plot.append(outplot)
      pl.savefig(outplot)
      pl.close(3)
      pl.figure(2)
      plot_xi02(xi02,mark='r-')

      reg=reg+1

   #compute south only and north + south correlations
   root2=rootN+'All-N'
   for ii in range(0,njack_S):
      root1=rootS+'JN'+str(ii)+'-S'
      #compute the correlation functions
      sbins, ns, mubins, nmu, DD,DR,RR, Wsum_data,Wsum_rand=load_PairCount_output(root1)
      Wrat=Wsum_data/Wsum_rand
      srebins, murebins, xi2drebin=compute_xi02(root1,sbins, ns, mubins, nmu,
                                  DD,DR,RR, Wrat,nscomb=nscomb,write=1)
      xi02=Xi_Legendre(root1,srebins, murebins, xi2drebin,outfold=outfold1,write=1)
      pl.figure(3)
      plot_xi02(xi02,mark='ro-',lab='S-'+str(ii))
      pl.figure(1)
      plot_xi02(xi02,mark='b-')

      #compute the combined correlation function
      sbins, ns, mubins, nmu, DD,DR,RR=combine_count(root1,root2)
      root=rootN+'COMB'+str(reg)+'-S'+str(ii)
      srebins, murebins, xi2drebin=compute_xi02(root,sbins, ns, mubins, nmu,
                                  DD,DR,RR, 1.0,nscomb=nscomb,write=1)
      xi02=Xi_Legendre(root,srebins, murebins, xi2drebin,outfold=outfoldcomb,write=1)
      pl.figure(3)
      plot_xi02(xi02,mark='bo-',lab='COMB')
      outplot='plots/Reg-S'+str(ii)+'.png'
      f_plot.append(outplot)
      pl.savefig(outplot)
      pl.close(3)
      pl.figure(2)
      plot_xi02(xi02,mark='b-')

      reg=reg+1

   root1=rootN+'All-N'
   sbins, ns, mubins, nmu, DD,DR,RR, Wsum_data,Wsum_rand=load_PairCount_output(root1)
   Wrat=Wsum_data/Wsum_rand
   srebins, murebins, xi2drebin=compute_xi02(root1,sbins, ns, mubins, nmu,
                                  DD,DR,RR, Wrat,nscomb=nscomb,write=1)
   xi02=Xi_Legendre(root1,srebins, murebins, xi2drebin,outfold=outfold1,write=1)
   pl.figure(3)
   plot_xi02(xi02,mark='ro-',lab=root1)
   pl.figure(1)
   plot_xi02(xi02,mark='r--',lw=1,lab=root1)

   root2=rootS+'All-S'
   sbins, ns, mubins, nmu, DD,DR,RR, Wsum_data,Wsum_rand=load_PairCount_output(root1)
   Wrat=Wsum_data/Wsum_rand
   srebins, murebins, xi2drebin=compute_xi02(root1,sbins, ns, mubins, nmu,
                                  DD,DR,RR, Wrat,nscomb=nscomb,write=1)
   xi02=Xi_Legendre(root1,srebins, murebins, xi2drebin,outfold=outfold1,write=1)
   pl.figure(3)
   plot_xi02(xi02,mark='bo-',lab=root2)
   pl.figure(1)
   plot_xi02(xi02,mark='b--',lw=1,lab=root2)

   sbins, ns, mubins, nmu, DD,DR,RR=combine_count(root1,root2)
   srebins, murebins, xi2drebin=compute_xi02('All-NS',sbins, ns, mubins, nmu,
                               DD,DR,RR, 1.0,nscomb=nscomb,write=1)
   xi02=Xi_Legendre('All-NS',srebins, murebins, xi2drebin,outfold=outfoldcomb,write=1)
   plot_xi02(xi02,mark='k--',lab=root2)
   outplot='plots/All-NS.png'
   f_plot.append(outplot)
   pl.savefig(outplot)
   pl.close(3)

   pl.figure(1)
   pl.title('Individual region, red=N, blue=S',fontsize=32)
   outplot='plots/JN-All-xi02.png'
   f_plot.append(outplot)
   pl.savefig(outplot)

   pl.figure(2)
   plot_xi02(xi02,mark='k--',lab='All-NS')
   pl.title('COMB xi JN region, red=N, blue=S',fontsize=32)
   outplot='plots/JN-COMB-All-xi02.png'
   f_plot.append(outplot)
   pl.savefig(outplot)

def plot_CMASS():
   file='/global/u1/s/shadaba/Projects/Cosmomc/dr11_Comb_8mpc.txt'
   covfile='/global/u1/s/shadaba/Projects/Cosmomc/sdss_DR11_cov_xi02_8Mpch-0-24.mat'
   data=np.loadtxt(file)
   cov=np.loadtxt(covfile)
   diag_err=np.sqrt(np.diagonal(cov))

   rr=data[:-1,0]
   rr2=rr*rr
 
   pl.errorbar(rr,rr2*data[:-1,1],yerr=rr2*diag_err[:24],fmt='o',color='k',markersize=5,label='CMASS')
   pl.errorbar(rr,rr2*data[:-1,2],yerr=rr2*diag_err[24:],fmt='o',color='k',markersize=5)

   #pl.plot(rr,rr2*data[:,1],'ko-',linewidth=2,markersize=4,label='CMASS')
   #pl.plot(rr,rr2*data[:,2],'ko-',linewidth=2,markersize=4)
   pl.legend()


if __name__ == "__main__":
  #To compute xi02 and xi2d in rmu sampling
 if(1):
   print sys.argv
   if(len(sys.argv)!=4):
      xidir='example/' ; froot='out'
      root=xidir+'PairCount/'+froot
      xi2droot=xidir+'XI2D/'+froot
      xi02root=xidir+'XI02/'+froot
   else:
      root=sys.argv[1]
      xi2droot=sys.argv[2]
      xi02root=sys.argv[3]

   xi2d_wp_xi02_froot_jn(root,xi2droot,xi02root,nscomb=1,samp='rmu',NJN=0)


 print 'List of files written:',
 for ff in f_write:
     print ff

 print 'List of plots created:',
 for pp in f_plot:
     print pp


