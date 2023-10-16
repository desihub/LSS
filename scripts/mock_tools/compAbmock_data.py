import numpy as np
from matplotlib import pyplot as plt
import sys
import fitsio
import os

lssdir = '/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/'
version = 'v0.6/blinded'

def getmeanmockxi(tp='LRG_ffa',zr='0.4_0.6',nmock=25,wt='default_FKP',mockmin=0):
    xil = []
    for i in range(mockmin,mockmin+nmock):
        dirxi = '/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/mock'+str(i)+'/'
        xi = np.loadtxt(dirxi+'/xi/smu/xipoles_'+tp+'_GCcomb_'+zr+'_'+wt+'_lin4_njack0_nran4_split20.txt').transpose()
        xil.append(xi)
    xi = sum(xil)/nmock
    err0 = np.zeros(len(xi))
    err2 = np.zeros(len(xi))
    err4 = np.zeros(len(xi))
    for i in range(mockmin,mockmin+nmock):
        dirxi = '/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/mock'+str(i)+'/'
        xii = np.loadtxt(dirxi+'/xi/smu/xipoles_'+tp+'_GCcomb_'+zr+'_'+wt+'_lin4_njack0_nran4_split20.txt').transpose()
        err0 += (xi[2]-xii[2])**2
        err2 += (xi[3]-xii[3])**2
        err4 += (xi[4]-xii[4])**2
    err0 = np.sqrt(err0/nmock)
    err2 = np.sqrt(err2/nmock)
    err4 = np.sqrt(err4/nmock)
    return xi,err0,err2,err4
    
def getmeanmockxi_diff(tp='LRG_ffa',zr='0.4_0.6',nmock=25,wt='default_FKP',wt2='default_FKP_addRF',mockmin=0):
    xil = []
    for i in range(mockmin,mockmin+nmock):
        dirxi = '/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/mock'+str(i)+'/'
        xi = np.loadtxt(dirxi+'/xi/smu/xipoles_'+tp+'_GCcomb_'+zr+'_'+wt+'_lin4_njack0_nran4_split20.txt').transpose()
        xi2 = np.loadtxt(dirxi+'/xi/smu/xipoles_'+tp+'_GCcomb_'+zr+'_'+wt2+'_lin4_njack0_nran4_split20.txt').transpose()
        xil.append(xi-xi2)
    xi = sum(xil)/nmock
    err0 = np.zeros(len(xi[0]))
    err2 = np.zeros(len(xi[0]))
    err4 = np.zeros(len(xi[0]))
    for i in range(mockmin,mockmin+nmock):
        dirxi = '/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/mock'+str(i)+'/'
        xii = np.loadtxt(dirxi+'/xi/smu/xipoles_'+tp+'_GCcomb_'+zr+'_'+wt+'_lin4_njack0_nran4_split20.txt').transpose()
        xi2 = np.loadtxt(dirxi+'/xi/smu/xipoles_'+tp+'_GCcomb_'+zr+'_'+wt2+'_lin4_njack0_nran4_split20.txt').transpose()
        err0 += (xi[2]-(xii[2]-xi2[2]))**2
        err2 += (xi[3]-(xii[3]-xi2[3]))**2
        err4 += (xi[4]-(xii[4]-xi2[4]))**2
    err0 = np.sqrt(err0/nmock)
    err2 = np.sqrt(err2/nmock)
    err4 = np.sqrt(err4/nmock)
   
    return xi,err0,err2,err4
    
def comp_datamockimsysdiff(tp='LRG',zr='0.4_0.6',wts=['noweight','RF']):
    plt.clf()
    mockwts = []
    datawts = []
    for wt in wts:
        if wt == 'noweight':
            mockwts.append('default_FKP')
            datawts.append('default_removeSN_FKP')
        if wt == 'RF':
            mockwts.append('default_FKP_addRF')
            datawts.append('default_swapinRF_FKP')
        if wt == 'SN':
            mockwts.append('default_FKP_addSN')
            datawts.append('default_FKP')
    mockmean = getmeanmockxi_diff(tp+'_ffa',zr=zr,wt=mockwts[0],wt2=mockwts[1])
    xi1 = np.loadtxt(lssdir+version+'/xi/smu/xipoles_'+tp+'_GCcomb_'+zr+'_'+datawts[0]+'_lin4_njack0_nran4_split20.txt').transpose()
    xi2 = np.loadtxt(lssdir+version+'/xi/smu/xipoles_'+tp+'_GCcomb_'+zr+'_'+datawts[1]+'_lin4_njack0_nran4_split20.txt').transpose()
    xidiff = xi1-xi2
    outf = lssdir+version+'/xi/smu//mockcomp/xipoles_'+tp+'_GCcomb_'+zr+wts[0]+wts[1]+'.png'
    plt.errorbar(xi1[0],xi1[0]*xidiff[2],xi1[0]*mockmean[1],fmt='ko')
    plt.plot(xi1[0],xi1[0]*mockmean[0][2],'k-')
    plt.errorbar(xi1[0],xi1[0]*xidiff[3],xi1[0]*mockmean[2],fmt='k^',alpha=0.5)
    plt.plot(xi1[0],xi1[0]*mockmean[0][3],'k--',alpha=0.5)
    plt.xlabel('s (Mpc/h)')
    plt.ylabel(r's ($\xi$'+wts[0]+r'-$\xi$'+wts[1]+')')
    plt.title(tp+' '+zr)
    plt.savefig(outf)
    plt.show()
    plt.clf()
    return True
   
