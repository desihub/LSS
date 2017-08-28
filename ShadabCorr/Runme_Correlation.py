#!/usr/bin/env python

"""
This program computes 2d correlation function.
"""

from __future__ import print_function,division
import numpy as np
import f
import time

import sys
import os
import argparse

#MULTI PROCESSING
import multiprocessing as mp
#import mycosmo as mcosm
#import General_function as GF

#Try to imports pandas if available
try:
   import pandas as pd
   pd_flag=1
except:
   pd_flag=0
   print('''pandas not found: you can install pandas from 
	    https://pypi.python.org/pypi/pandas/0.18.0/#downloads 
	    The code will work without pandas but pandas will speed up loading large text files
	    ''')


__author__ = "Hongyu Zhu, Shadab Alam"
__version__ = "1.0"
__email__  = "shadaba@andrew.cmu.edu"

parser = argparse.ArgumentParser(description='Calculate 2d correlation function:')
#these should be input
parser.add_argument('-plots',type=int,default=1)
parser.add_argument('-sampmode' ,type=int,default=0,
                  help='''select the sampling mode 0=rmu , 1=rpara-rperp, 2=rtheta, 3=log(rpar)-rperp (The range of r_perp in samplim should be entered in the log space), 4=log(r)-theta (The range of r in samplim should be entered in the log space)''') 
parser.add_argument('-data' ,default='example/data.txt',help='''data file with columns being  X,Y,Z,weight,JN,Vx,Vy,Vz
The minimum columns should be X,Y,Z,weight. X,Y,Z should be in Mpc/h and Vx,Vy,Vz is in simulation unit used incase of RSD''')
parser.add_argument('-rand' ,default='example/randoms.txt',help='same as data files. No velocity is used for randoms.')
parser.add_argument('-data2' ,default='',help='second data file for cross-correlation, same structure as first data file')
parser.add_argument('-rand2' ,default='',help='second random file for cross-correlation file')
parser.add_argument('-nbins'  ,nargs='+',type=int,default=[32,100],help='number of bins in xi sampling')
parser.add_argument('-samplim',nargs='+',type=float, default=[0,160,0,1],help='The range for 2d correlation function calculation depending on the modes it will take lower and higher limit for the two axis')
parser.add_argument('-pbc'  ,type=int,default=0, help='This is to say whether the PBC should be applied or not') 
parser.add_argument('-los'  ,type=int,default=0, help='''This is to choose a definition for line of sight. 0=los is defined as mid point of the pair, 1=los is defined along the Z axis,\n(only for cross correlation 2,3) \n2=los is define along the galaxies in the first data file, 3= second data file''') 
parser.add_argument('-njn'  ,type=int,default=0, help='0 no jacknife column, >0 number of jacknife regions') 
parser.add_argument('-outfile' ,default='example/PairCount/out',help="set the output file name")

parser.add_argument('-interactive'  ,type=int, default=1, help="set 1 to print some information while running") # set zero to not print out things

parser.add_argument('-filetype'  ,default='txt',help="The script can read few different file type: txt, polartxt, fhandle,fits") 
parser.add_argument('-z1z2'  ,nargs='+',type=float,default=[0,100], help='To restrict the redshift range if RA,DEC,Redshift is provided.') 
parser.add_argument('-randfactor'  ,type=int,default=0, help='''This is used incase random file is not provided for periodic box to generate uniform random in the box with a fixed seed. The number of randoms is set as the randomfactor times the number of data''')

parser.add_argument('-H0'  ,type=float,default=67.6, help='Hubble constant when cosmology is needed') 
parser.add_argument('-omM'  ,type=float,default=0.315, help='matter density  when cosmology is needed') 

parser.add_argument('-nproc'  ,type=int,default=2, help='number of multiprocess to run') 
parser.add_argument('-nnode'  ,type=int,default=1, help='number of nodes to run') 
parser.add_argument('-nodeid'  ,type=int,default=0, help='nodeid for this job')

args = parser.parse_args()
#print(args)

#basic input checks
assert len(args.samplim)==4
assert len(args.nbins) ==2
if(args.pbc==1 and args.los !=1):
   print("*** INVALID INPUT combination ***")
   print("periodic boundary condition can be applied only with los=1(along z axis)")
   sys.exit()

if(args.los==2 or args.los==3):
   if(args.data2=='' or args.rand2==''):
      print("*** INVALID INPUT combination ***")
      print("The los=2,3 can be used only for corss correlation")
      print("Please provide second data and random file for cor-correlation")
      sys.exit()



#sampmode to sampcode relation
sampcodes=['rmu','rp-pi','rtheta','logrp-pi','logr-theta']

rlim  = np.array(args.samplim,dtype='double')
nbins = np.array(args.nbins,dtype='int')
nhocells = 200

t0 = time.time() #start time

def getminmax(data,rand=''):
   #determin blen
   blen  = np.array([0,0,0],dtype='double')
   POS_min = np.array([0,0,0],dtype='double')
   POS_max = np.array([0,0,0],dtype='double')
   for ii in range(0,3):
      #compute the minimum
      d1=np.min(data[:,ii])
      try:
         d2=np.min(rand[:,ii])
         POS_min[ii]=min(d1,d2)
      except:
	 POS_min[ii]=d1
      #compute the maximum
      d1=np.max(data[:,ii])
      try:
         d2=np.max(rand[:,ii])
         POS_max[ii]=max(d1,d2)
      except:
	 POS_max[ii]=d1
      #compute the blen
      blen[ii]=POS_max[ii]-POS_min[ii]

   return POS_min,POS_max, blen

#This uses pandas if available otherwise numpy loadtxt to read files
def load_txt(fname):
   if(pd_flag==1):
      tb=pd.read_table(fname,delim_whitespace=True,comment='#',header=None)
      tb=tb.as_matrix()
   else:
      tb=np.loadtxt(fname)
   return tb

#load the data file
def prep_data_rand(datafile,randfile,RSD=0):
   if(args.interactive>1):
      print("Loading files:\n "+datafile+" \n "+randfile)

   if(args.filetype=='txt' or args.filetype=='polartxt'):
      data=load_txt(datafile)
      if(args.filetype=='polartxt'):
         import mycosmo as mcosm
	 XYZ,interp=mcosm.RDZ2XYZ(data[:,0],data[:,1],data[:,2],args.H0,args.omM,1-args.omM,interp='')
         data=np.column_stack([XYZ,np.ones(XYZ.shape[0])])
      if(os.path.isfile(randfile)):
         rand=load_txt(randfile)
         if(args.filetype=='polartxt'):
            import mycosmo as mcosm
	    XYZ,interp=mcosm.RDZ2XYZ(rand[:,0],rand[:,1],rand[:,2],args.H0,args.omM,1-args.omM,interp='')
            rand=np.column_stack([XYZ,np.ones(XYZ.shape[0])])

	 #selecting subsample of random
	 if(args.randfactor>0):
	    nrand=args.randfactor*data.shape[0]
            idx=np.random.choice(np.arange(rand.shape[0]),nrand,replace=False)
	    print('nrand selected: ',rand.shape[0],nrand)
	    rand=rand[idx,:]
   elif(args.filetype=='fhandle'):
      if (args.pbc==1):
         data=readBOXwithifhandle(datafile)
      else:
         data,rand=readSKYwithifhandle(datafile,randfile)
   elif(args.filetype=='fits'):
      data,rand=readSKYwithifits(datafile,randfile)
   else:
      print('unknown file type')
      sys.exit()

   #generate randoms with fixed seed for periodic box if random file 
   #is not provides
   if(args.pbc==1 and randfile==''):
      print("random file not given, Generating internal random")
      nrand=np.int(args.randfactor*data.shape[0])
      rand=randoms_fixedseed(data,nrand)

   #add jacknife region if its pbc, los=zaxis and data have only x,y,z,wt
   if(args.pbc==1 and args.los==1 and data.shape[1]==4):
      data,rand=add_pbc_jncol(data,rand)


   ngal =np.int32(data.shape[0])
   nrand=np.int32(rand.shape[0])
   if(args.interactive>=0):
      print("%10d galaxies from %s"%(ngal,datafile))
      print("%10d randoms from %s"%(nrand,randfile))
      print("\nTime:%d sec"%np.int(time.time() - t0))

   POS_min,POS_max, blen=getminmax(data,rand=rand)

   if(args.interactive>1):
      print("Extent of the data/random:")
      print("MIN: ",', '.join('%8.2f'%x for x in POS_min) )
      print("MAX: ",', '.join('%8.2f'%x for x in POS_max) )
      print("Length: ",', '.join('%8.2f'%x for x in blen),'\n' )

   #compute sum of weights and write norm file
   if(args.njn==0):
      sumwtdata=np.sum(data[:,3])
      sumwtrand=np.sum(rand[:,3])
   else:
      sumwtdata=np.zeros(args.njn+1)
      sumwtrand=np.zeros(args.njn+1)
      #all weight
      sumwtdata[args.njn]=np.sum(data[:,3])
      sumwtrand[args.njn]=np.sum(rand[:,3])
      #jn weights
      for ii in range(0,args.njn):
	 ind=data[:,4]==ii
	 sumwtdata[ii]= sumwtdata[args.njn]-np.sum(data[ind,3])
	 ind=rand[:,4]==ii
	 sumwtrand[ii]= sumwtrand[args.njn]-np.sum(rand[ind,3])

   #convert arrays to contiguous array
   if(args.njn==0):
      data_c = np.ascontiguousarray(data[:,[0,1,2,3]], dtype=np.float64)
      rand_c = np.ascontiguousarray(rand[:,[0,1,2,3]], dtype=np.float64)
   else:
      data_c = np.ascontiguousarray(data[:,[0,1,2,3,4]], dtype=np.float64)
      rand_c = np.ascontiguousarray(rand[:,[0,1,2,3,4]], dtype=np.float64)
      if(np.max(data[:,4]) != args.njn-1 or np.max(rand[:,4]) != args.njn-1):
	 print("number of jacknife given and in the file do not match")
	 print("file njn data,random:  %d %d given njn: %d"
	       %(np.max(data[:,4]),np.max(rand[:,4]),args.njn))
	 sys.exit()

   return data_c, rand_c, blen, POS_min, POS_max, sumwtdata, sumwtrand

def readBOXwithifhandle(datafile,Lbox=3000):
   '''This function reads file with Martins filehandle, assumes cubic box and read only XYZ co-ordinate,
   set weight =1 and in case of jacknife creates equalsize regions, also generate uniform randoms in the box size,
   Also assumes co-ordinate to be in between zeros and 1 which is scaled by Lbox in Mpc/h'''

   #import Martins file handler
   import ndfilehandler as FH
   data = FH.read_file(datafile)['pos']*Lbox
   data=np.column_stack([data,np.ones(data.shape[0])])
   return data

def readSKYwithifhandle(datafile,randfile):
   '''This function reads file with Martins filehandle, 
   assumes SKY projected mocks read only RA,DEC and Z,
   set weight =1 and in do not support jacknife at the moment,
   Assumes cosmology '''
   H0=68; omM=0.30; omL=1-omM

   #import Martins file handler
   import ndfilehandler as FH
   #read RA,DEC and redshift
   RA  = FH.read_file(datafile)['RA']
   DEC = FH.read_file(datafile)['DEC']
   #Red = FH.read_file(datafile)['ZR']
   Red = FH.read_file(datafile)['Z']

   XYZ,interp=mcosm.RDZ2XYZ(RA,DEC,Red,H0,omM,omL,interp='')
   data=np.column_stack([XYZ,np.ones(XYZ.shape[0])])

   #Working on randoms
   RA  = FH.read_file(randfile)['RA']
   DEC = FH.read_file(randfile)['DEC']
   #Red = FH.read_file(randfile)['ZR']
   Red = FH.read_file(randfile)['Z']
   #Setup subsample random
   nrand=np.int(args.randfactor*data.shape[0])
   idx = np.random.choice(np.arange(RA.size), nrand, replace=False)
   print('nrand:',RA.size,nrand)
   #subsampled RA,DEC and Red
   RA=RA[idx]; DEC=DEC[idx];Red=Red[idx] 
   print('RED:',np.min(Red) , np.max(Red), np.mean(Red))

   XYZ,interp=mcosm.RDZ2XYZ(RA,DEC,Red,H0,omM,omL,interp='')
   rand=np.column_stack([XYZ,np.ones(XYZ.shape[0])])

   return data,rand


def readSKYwithifits(datafile,randfile):
   '''This function reads file with fits format, 
   assumes SKY projected mocks read only RA,DEC and Z,
   set weight =1 and in do not support jacknife at the moment,
   Assumes cosmology '''
   H0=68; omM=0.30; omL=1-omM; #z=3.0; dz=0.25

   #fixed JN files for now
   dir='/home/shadaba/Projects/DESI/Jacknife/'
   Nfile='North-v1.0_qso_10rand-zlim-0.0-4.0-NorthBorder_RA-DEC_JN-11-4_v1.txt'
   Sfile='South-v1.0_qso_10rand-zlim-0.0-4.0-SouthBorder_RA-DEC_JN-5-4_v1.txt'
   JNfile=[dir+Nfile,dir+Sfile]

   #import fits file handler
   import fitsio as F
   findat=F.FITS(datafile,'r')
   #Zsel="ZR >= "+str(args.z1z2[0])+" && ZR <"+str(args.z1z2[1])
   Zsel="Z >= "+str(args.z1z2[0])+" && Z <"+str(args.z1z2[1])
   indsel=findat[1].where(Zsel)

   #read RA,DEC and redshift
   RA  = findat[1]['RA'][indsel]
   DEC = findat[1]['DEC'][indsel]
   #Red = findat[1]['ZR'][indsel]
   Red = findat[1]['Z'][indsel]

   XYZ,interp=mcosm.RDZ2XYZ(RA,DEC,Red,H0,omM,omL,interp='')
   data=np.column_stack([XYZ,np.ones(XYZ.shape[0])])
   if(args.njn>0):
      JN_reg=GF.compute_JN_data(RA,DEC,JNfile)
      data=np.column_stack([data,JN_reg])

   #Working on randoms
   finran=F.FITS(randfile,'r')
   Zsel="Z >= "+str(args.z1z2[0])+" && Z <"+str(args.z1z2[1])
   indsel=finran[1].where(Zsel)

   #read RA,DEC and redshift
   RA  = finran[1]['RA'][indsel]
   DEC = finran[1]['DEC'][indsel]
   Red = finran[1]['Z'][indsel]

   #Setup subsample random
   nrand=np.int(args.randfactor*data.shape[0])
   print('nrand:',RA.size,nrand)
   if(RA.size>nrand):
      idx = np.random.choice(np.arange(RA.size), nrand, replace=False)
      #subsampled RA,DEC and Red
      RA=RA[idx]; DEC=DEC[idx];Red=Red[idx] 
   print('RED:',np.min(Red) , np.max(Red), np.mean(Red))

   XYZ,interp=mcosm.RDZ2XYZ(RA,DEC,Red,H0,omM,omL,interp='')
   rand=np.column_stack([XYZ,np.ones(XYZ.shape[0])])
   if(args.njn>0):
      JN_reg=GF.compute_JN_data(RA,DEC,JNfile)
      rand=np.column_stack([rand,JN_reg])


   return data,rand

##
def randoms_fixedseed(data,nrand):
   '''box size is decide by the data'''

   np.random.seed(600)
   POS_min,POS_max, blen=getminmax(data)
   #check if all the dimension have same lengtj
   assert (np.abs(blen[0]-blen[1])<1e-4*blen[0] and 
	 np.abs(blen[1]-blen[2])<1e-3*blen[0]), "Not a cubic box all dimensions are not same: %12.8lf %12.8lf %12.8lf "%(blen[0],blen[1],blen[2])

   Lbox=blen[1]
   rand=np.random.random(nrand*3).reshape(nrand,3)*Lbox
   rand=np.column_stack([rand,np.ones(nrand)])

   return rand

def add_pbc_jncol(data,rand):
   '''If the input is a periodic box and los is along z axis then jacknife region is simply equal area region in the x-y space which can be done in using this function and not needed to be supplied with data file make sure that njn is a perfect square'''

   #adding jacknife regions
   if(args.njn>0 and args.los==1):
      POS_min,POS_max, blen=getminmax(data,rand=rand)
      NJNx=np.int(np.sqrt(args.njn))
      NJNy=np.int(args.njn/NJNx)
      for ii in (0,2):
	 if(ii==0): mat=data
	 else: mat=rand
	 
         #get the x and y indx as integers
	 indx=np.int_(NJNx*(mat[:,0]-POS_min[0])/blen[0])
	 indy=np.int_(NJNy*(mat[:,1]-POS_min[1])/blen[1])
         #apply modulo operation on x an y index
	 indx=np.mod(indx,NJNx)
	 indy=np.mod(indy,NJNy)

	 #convert index to integers
	 #indx.astype(np.int64); indy.astype(np.int64);
	 jnreg=NJNy*indx+indy
	 mat=np.column_stack([mat,jnreg])

         if(ii==0): data=mat
	 else: rand=mat

      return data,rand
   else:
      print('not appropriate input to add jacknife internally')
      sys.exit()
      return 0


def write_pair_count(pair_count,outfile):
   xsamp=np.linspace(rlim[0],rlim[1],nbins[0]+1)
   ysamp=np.linspace(rlim[2],rlim[3],nbins[1]+1)
   with file(outfile, 'w') as fwrite:
      xsamp_str=''
      for xx in xsamp:
         xsamp_str="%s %10.5f"%(xsamp_str,xx)
      fwrite.write(xsamp_str+'\n')

      ysamp_str=''
      for yy in ysamp:
         ysamp_str="%s %10.5f"%(ysamp_str,yy)
      fwrite.write(ysamp_str+'\n')

      np.savetxt(fwrite,pair_count,fmt='% 25.15e')
   return 0


def collect_pcnodes(type='DD',cleandir=0):
   for ii in range(0,args.nnode):
      if(args.njn==0):
	 tag=''
      else:
	 tag='-All'

      outroot=args.outfile+'_nodeid'+str(ii)
      JNdir=outroot+'-'+sampcodes[args.sampmode]+'-JNdir/'
      pcfile=outroot+tag+'-'+sampcodes[args.sampmode]+'-'+type+'.dat'
      data=np.loadtxt(pcfile,skiprows=2)
      if(ii==0):
	 pcsum=data
      else:
         pcsum=pcsum+data
      #remove the file
      os.system('rm -f '+pcfile)
     
      if(args.njn>0):
	 if(ii==0):
	    pcsumJN=np.zeros(data.shape[0]*data.shape[1]*(args.njn+1)).reshape(
		                 data.shape[0],data.shape[1],(args.njn+1))
	 for jj in range(0,args.njn):
            pcfile=JNdir+'jk-'+str(jj)+'-'+type+'.dat'
            data=np.loadtxt(pcfile,skiprows=2)
	    pcsumJN[:,:,jj]=pcsumJN[:,:,jj]+data

	 #remove the directoty of node and JN for auto correlation
	 if(cleandir==1):
            os.system('rm -r '+JNdir)

   #end of node loop
   if(args.njn>0):
      pcsumJN[:,:,args.njn]=pcsum;
      pcsum=pcsumJN
      
   write_output(pcsum,type=type,combnode=1)
   return 0

def combine_clean(corrtype='auto'):
   #check if all the nodes are done
   t0=time.time()
   check=0
   while(check!=args.nnode-1):
      check=0;
      for ii in range(1,args.nnode):
         donefile=args.outfile+'-'+sampcodes[args.sampmode]
         donefile=donefile+'_nodeid'+str(ii)+'.done'
	 if(os.path.isfile(donefile)):
            check=check+1
      if(args.interactive>1):
         print('checking node %d nodes finished out of %d: %d sec'%(
	    check,args.nnode,time.time()-t0))

      if(check!=args.nnode-1):
	 time.sleep(2)
      else:
	 print('This is node 0, Waiting time for other nodes:%d sec'%(time.time()-t0))
  
   #Collect the pair counts
   if(corrtype=='auto'):
      collect_pcnodes(type='DD',cleandir=0)
      collect_pcnodes(type='DR',cleandir=0)
      collect_pcnodes(type='RR',cleandir=1) #RR must be called in the end after DD and DR call
   elif(corrtype=='cross'):
      collect_pcnodes(type='D1D2',cleandir=0)
      collect_pcnodes(type='D1R2',cleandir=0)
      collect_pcnodes(type='R1D2',cleandir=0)
      collect_pcnodes(type='R1R2',cleandir=1) #RR must be called in the end after DD and DR call
   else:
      print('Invalid corrtype = %s'%corrtype)
      sys.exit()

   #clean the donefiles
   if(args.interactive>1):
      print('clenaing the files for individual nodes')
   for ii in range(1,args.nnode):
      donefile=args.outfile+'-'+sampcodes[args.sampmode]
      donefile=donefile+'_nodeid'+str(ii)+'.done'
      os.system('rm -f '+donefile)
      

def write_output(xc,type='DD',SDwt=np.zeros(1),SRwt=np.zeros(1),SDwt2=np.zeros(1),SRwt2=np.zeros(1),combnode=0):
   if(args.nnode==1 or combnode==1):
      outroot=args.outfile
   else:
      outroot=args.outfile+'_nodeid'+str(args.nodeid)

   if(args.njn==0):
      outfile=outroot+'-'+sampcodes[args.sampmode]+'-'+type+'.dat'
      write_pair_count(xc,outfile)
      #write the norm file
      if(SDwt!=0 and SRwt !=0 and args.nodeid==0): #only for first node
         outfile=args.outfile+'-'+sampcodes[args.sampmode]+'-norm.dat'
         with file(outfile, 'w') as fwrite:
            fwrite.write(args.data+': %15.10e\n'%(SDwt))
            fwrite.write(args.rand+': %15.10e'%(SRwt))
            if(SDwt2!=0 and SRwt2 !=0): #writing the weights of second file
               fwrite.write('\n'+args.data2+': %15.10e\n'%(SDwt2))
               fwrite.write(args.rand2+': %15.10e'%(SRwt2))
   else:
      #first write the all pair count
      outfile=outroot+'-All-'+sampcodes[args.sampmode]+'-'+type+'.dat'
      xcAll=xc[:,:,args.njn]
      write_pair_count(xcAll,outfile)
      #write the norm file
      if(SDwt.size==args.njn+1 and args.nodeid==0): #only for first node
         outfile=args.outfile+'-All-'+sampcodes[args.sampmode]+'-norm.dat'
         with file(outfile, 'w') as fwrite:
            fwrite.write(args.data+': %15.10e\n'%(SDwt[args.njn]))
            fwrite.write(args.rand+': %15.10e'%(SRwt[args.njn]))
            if(SDwt2.size==args.njn+1): #only for first node
               fwrite.write('\n'+args.data2+': %15.10e\n'%(SDwt2[args.njn]))
               fwrite.write(args.rand2+': %15.10e'%(SRwt2[args.njn]))

      #now write the each jacknife in the folder
      JNdir=outroot+'-'+sampcodes[args.sampmode]+'-JNdir/'
      if(not os.path.isdir(JNdir)):
         os.system('mkdir '+JNdir)
      if(args.nodeid==0):
         JNdirnorm=args.outfile+'-'+sampcodes[args.sampmode]+'-JNdir/'
         if(not os.path.isdir(JNdirnorm)):
            os.system('mkdir '+JNdirnorm)

      for ii in range(0,args.njn):
	 outfile=JNdir+'jk-'+str(ii)+'-'+type+'.dat'
         xcJN=xc[:,:,ii]
	 if(combnode==0):
	    xcJN_All=xcAll-xcJN
            write_pair_count(xcJN_All,outfile)
	 else:
            write_pair_count(xcJN,outfile)

         #write the norm file
         if(SDwt.size==args.njn+1 and args.nodeid==0):
            outfile=JNdirnorm+'jk-'+str(ii)+'-norm.dat'
            with file(outfile, 'w') as fwrite:
               fwrite.write(args.data+': %15.10e\n'%(SDwt[ii]))
               fwrite.write(args.rand+': %15.10e'%(SRwt[ii]))
               if(SDwt2.size==args.njn+1):
                  fwrite.write('\n'+args.data2+': %15.10e\n'%(SDwt2[ii]))
                  fwrite.write(args.rand2+': %15.10e'%(SRwt2[ii]))

   return 0


def mp_pair_count(data1_c, data2_c, rlim_c, nbins, nhocells, blen_c, pos_min_c):
   def corr2d_warp(data1_c, data2_c, rlim_c, nbins, nhocells, blen_c, pos_min_c,argsin,out_q):
      #pc={}
      #pc[mp.current_process()] 
      pc= f.corr2dpy(data1_c, data2_c, rlim_c, nbins[0],nbins[1], nhocells, 
                 blen_c, pos_min_c, argsin.sampmode, argsin.njn, 
		 argsin.pbc,argsin.los,argsin.interactive)
      if(argsin.interactive>1):
         print('using pid:',os.getpid())
      #print(pc)
      #outdict={}
      #outdict[mp.current_process()]=pc
      out_q.put(pc)
   
   ndata=data1_c.shape[0]
   if(args.interactive>1):
      print('creating queue')
   # Each process will get 'chunksize' nums and a queue to put his out
   # dict into
   out_q = mp.Queue()
   chunknode = int(np.ceil(ndata / float(args.nnode)))
   chunksize = int(np.ceil(chunknode / float(args.nproc)))
   procs = []
 
   nodebeg=args.nodeid*chunknode
   if(args.interactive>1):
      print('submitting processes:')
   for ii in range(args.nproc):
      ind1=chunksize*ii+nodebeg
      ind2=chunksize*(ii+1)+nodebeg
      if(ii==args.nproc-1):
	 if(args.nodeid==args.nnode-1):
            ind2=ndata
	 else:
	    ind2=chunknode+nodebeg

      if(args.interactive>2):
	 print("submit: node= %d ,proc= %d : %d %d %d"%(
	    args.nodeid,ii,ind1,ind2, ind2-ind1))
      
      p = mp.Process(
                target=corr2d_warp,
                args=(data1_c[ind1:ind2,:],data2_c, rlim_c, nbins, 
                nhocells, blen_c, pos_min_c,args,out_q))
      procs.append(p)
      p.start()

   # Collect all results into a single result dict. We know how many dicts
   # with results to expect.
   resultdict = {}
   for ii in range(args.nproc):
      #resultdict[ii]=out_q.get(timeout=3600.)
      resultdict[ii]=out_q.get()
      if(args.interactive>2):
         print('collect proc:',ii)

      if(ii==0):
	 pcAll=resultdict[ii]
      else:
	 pcAll=pcAll+resultdict[ii]

   if(args.interactive>1):
      print('joining results')
   # Wait for all worker processes to finish
   for p in procs:
      p.join()

   return pcAll

#to combine the min and max extent of data for two sample
def combine_survey(POS1_min,POS1_max,blen1,POS2_min,POS2_max,blen2):
   diffmin=np.zeros(3); diffmax=np.zeros(3)
   for ii in range(0,3):
      diffmin[ii]=np.abs(POS1_min[ii]-POS2_min[ii])
      diffmax[ii]=np.abs(POS1_max[ii]-POS2_max[ii])

      POS1_min[ii]=min(POS1_min[ii],POS2_min[ii])
      POS1_max[ii]=max(POS1_max[ii],POS2_max[ii])

      #compute the blen
      blen1[ii]=POS1_max[ii]-POS1_min[ii]

   if(args.interactive>1):
      print("Extent of the data/random combined:")
      print("MIN: ",', '.join('%8.2f'%x for x in POS1_min) )
      print("MAX: ",', '.join('%8.2f'%x for x in POS1_max) )
      print("Length: ",', '.join('%8.2f'%x for x in blen1),'\n' )

   return POS1_min,POS1_max,blen1

def compute_auto(data_c, rand_c, sumwtdata,sumwtrand,rlim_c, nbins, blen_c, pos_min_c):
   #compute DD pair count and write to file
   print('\nWorkind on DD (. every 100k points):')
   t1 = time.time() #pair count start time
   xcDD=mp_pair_count(data_c, data_c, rlim_c, nbins, nhocells, blen_c, pos_min_c)
   write_output(xcDD,type='DD',SDwt=sumwtdata, SRwt=sumwtrand)
   if(args.interactive>0):
      print("\nFinished DD in %d sec with %d process and %d %d particles, max scale %10.4f"%(
	 np.int(time.time() - t1),args.nproc,data_c.shape[0],data_c.shape[0],rlim_c[1]))

   #compute DR pair count and write to file
   print('\nWorkind on DR (. every 100k points):')
   xcDR=mp_pair_count(data_c, rand_c, rlim_c, nbins, nhocells, blen_c, pos_min_c)
   write_output(xcDR,type='DR')
   if(args.interactive>0):
      print("\nFinished DR in %d sec with %d process and %d %d particles, max scale %10.4f"%(
	 np.int(time.time() - t1),args.nproc,data_c.shape[0],rand_c.shape[0],rlim_c[1]))

   #compute RR pair count and write to file
   print('\nWorkind on RR (. every 100k points):')
   xcRR=mp_pair_count(rand_c, rand_c, rlim_c, nbins, nhocells, blen_c, pos_min_c)
   write_output(xcRR,type='RR')
   if(args.interactive>0):
      print("\nFinished RR in %d sec with %d process and %d %d particles, max scale %10.4f"%(
	 np.int(time.time() - t1),args.nproc,rand_c.shape[0],rand_c.shape[0],rlim_c[1]))

   #if(args.interactive>0):
   print("\nFinished job in %d sec with %d process and nrand %d ,max scale %10.4f"%(np.int(time.time() - t0),args.nproc,rand_c.shape[0],rlim_c[1]))

   if(args.nodeid!=0):
      #write a file to indicate that it has finished
      donefile=args.outfile+'-'+sampcodes[args.sampmode]
      donefile=donefile+'_nodeid'+str(args.nodeid)+'.done'
      os.system('touch '+donefile)
   elif(args.nnode>1): 
      #wait until all other nodes are done and then combine 
      #results from all nodes and clean the disk
      combine_clean(corrtype='auto')

   return 0
#end of compute_auto

def compute_cross(data_c, rand_c, sumwtdata,sumwtrand,
      data2_c, rand2_c, sumwtdata2,sumwtrand2, rlim_c, nbins, blen_c, pos_min_c):

   #compute DD pair count and write to file
   print('\nWorkind on D1D2 (. every 100k points):')
   t1 = time.time() #pair count start time
   xcDD=mp_pair_count(data_c, data2_c, rlim_c, nbins, nhocells, blen_c, pos_min_c)
   write_output(xcDD,type='D1D2',SDwt=sumwtdata, SRwt=sumwtrand,SDwt2=sumwtdata2, SRwt2=sumwtrand2)
   if(args.interactive>0):
      print("\nFinished D1D2 in %d sec with %d process and %d %d particles, max scale %10.4f"%(
	 np.int(time.time() - t1),args.nproc,data_c.shape[0],data2_c.shape[0],rlim_c[1]))

   #compute D1R2 pair count and write to file
   print('\nWorkind on D1R2 (. every 100k points):')
   xcDR=mp_pair_count(data_c, rand2_c, rlim_c, nbins, nhocells, blen_c, pos_min_c)
   write_output(xcDR,type='D1R2')
   if(args.interactive>0):
      print("\nFinished D1R2 in %d sec with %d process and %d %d particles, max scale %10.4f"%(
	 np.int(time.time() - t1),args.nproc,data_c.shape[0],rand2_c.shape[0],rlim_c[1]))

   #compute D1R2 pair count and write to file
   print('\nWorkind on R1D2 (. every 100k points):')
   xcDR=mp_pair_count(rand_c, data2_c, rlim_c, nbins, nhocells, blen_c, pos_min_c)
   write_output(xcDR,type='R1D2')
   if(args.interactive>0):
      print("\nFinished R1D2 in %d sec with %d process and %d %d particles, max scale %10.4f"%(
	 np.int(time.time() - t1),args.nproc,rand_c.shape[0],data2_c.shape[0],rlim_c[1]))

   #compute RR pair count and write to file
   print('\nWorkind on R1R2 (. every 100k points):')
   xcRR=mp_pair_count(rand_c, rand2_c, rlim_c, nbins, nhocells, blen_c, pos_min_c)
   write_output(xcRR,type='R1R2')
   if(args.interactive>0):
      print("\nFinished RR in %d sec with %d process and %d %d particles, max scale %10.4f"%(
	 np.int(time.time() - t1),args.nproc,rand_c.shape[0],rand2_c.shape[0],rlim_c[1]))

   #if(args.interactive>0):
   print("\nFinished job in %d sec with %d process and nrand %d ,max scale %10.4f"%(np.int(time.time() - t0),args.nproc,rand_c.shape[0],rlim_c[1]))

   if(args.nodeid!=0):
      #write a file to indicate that it has finished
      donefile=args.outfile+'-'+sampcodes[args.sampmode]
      donefile=donefile+'_nodeid'+str(args.nodeid)+'.done'
      os.system('touch '+donefile)
   elif(args.nnode>1): 
      #wait until all other nodes are done and then combine 
      #results from all nodes and clean the disk
      combine_clean(corrtype='cross')

   return 0
#end of compute_cross


if __name__=="__main__":  
   #load and prepare the data 
   data_c, rand_c, blen, POS_min, POS_max, sumwtdata,sumwtrand =prep_data_rand(args.data,args.rand,RSD=0)

   if(args.data2!='' and args.rand2!=''):
      data2_c, rand2_c, blen2, POS2_min, POS2_max, sumwtdata2,sumwtrand2 =prep_data_rand(args.data2,args.rand2,RSD=0)
      POS_min,POS_max,blen=combine_survey(POS_min,POS_max,blen,POS2_min,POS2_max,blen2)

   rlim_c =np.ascontiguousarray(rlim)
   blen_c =np.ascontiguousarray(blen)
   pos_min_c=np.ascontiguousarray(POS_min)

   if(args.interactive>=0):
      print("\nPrepared data and random:%d sec"%np.int(time.time() - t0))

   if(args.data2==''):
      compute_auto(data_c, rand_c, sumwtdata,sumwtrand,rlim_c, nbins, blen_c, pos_min_c)
   else:
      compute_cross(data_c, rand_c, sumwtdata,sumwtrand,data2_c, rand2_c,
	    sumwtdata2,sumwtrand2,rlim_c, nbins, blen_c, pos_min_c)

