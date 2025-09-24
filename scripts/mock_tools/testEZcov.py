import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys
import fitsio
import os




def cat_ells(xis,indrange=None):
    if indrange is None:
        indrange = [0,len(xis[0])]
    xil = [xis[0][indrange[0]:indrange[1]]]
    for i in range(1,len(xis)):
        xil.append(xis[i][indrange[0]:indrange[1]])
    xic = np.concatenate(xil)
    return xic

def get_xi_cov_desipipe_baseline_txt(smin=20,smax=200,zr='0.4-0.6',tp='LRG',rec='recon_recsym',ells=[0,2,4],Nmock=1000,flavor='ffa',mockversion='v1',thetacut='',mocktype='EZmock'):
    from pycorr import TwoPointCorrelationFunction
    dirm = '/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/'+mocktype+'/desipipe/'+mockversion+'/'+flavor+'/baseline_2pt/mock'
    xil = []
    sepl = []
    start = 1    
    
    fnm = dirm +'1/'+rec+'/xi/smu/allcounts_'+tp+'_GCcomb_z'+zr+thetacut+'_d4_poles.txt' #loading first to get binning setup
    
    rbs = 4
    indmin = smin//rbs
    indmax = smax//rbs                      
    rebinned = np.loadtxt(fnm).transpose()
    
    xis = []
    for ell in ells:
        xis.append(rebinned[2+ell//2])
    xin0 = cat_ells(xis,indrange=[indmin,indmax])    
        
    nbin = len(xin0)
    print(nbin)
    xiave = np.zeros((nbin))
    cov = np.zeros((nbin,nbin))

    Ntot = 0
    fac = 1.
    for i in range(start,start+Nmock):
        nr = str(i)
        xinpy = dirm +str(i)+'/'+rec+'/xi/smu/allcounts_'+tp+'_GCcomb_z'+zr+thetacut+'_d4_poles.txt' 
        if os.path.isfile(xinpy):
            rebinned = np.loadtxt(xinpy).transpose()
            xis = []
            for ell in ells:
                xis.append(rebinned[2+ell//2])

            xic = cat_ells(xis,indrange=[indmin,indmax])            
            xiave += xic
            Ntot += 1.
    print( Ntot)        
    xiave = xiave/float(Ntot)
    for i in range(1,Nmock+1):
        nr = str(i)
        xinpy = dirm +str(i)+'/'+rec+'/xi/smu/allcounts_'+tp+'_GCcomb_z'+zr+thetacut+'_d4_poles.txt'
        if os.path.isfile(xinpy):
            rebinned = np.loadtxt(xinpy).transpose()
            xis = []
            for ell in ells:
                xis.append(rebinned[2+ell//2])
            
            xic = cat_ells(xis,indrange=[indmin,indmax])
            for j in range(0,nbin):
                xij = xic[j]
                for k in range(0,nbin):
                    xik = xic[k]
                    cov[j][k] += (xij-xiave[j])*(xik-xiave[k])

    cov = cov/float(Ntot)                   
        
    return xiave,cov

def get_xiave_desipipe_ab_baseline(zr='0.4-0.6',tp='LRG',rec='recon_recsym',nmock=25,flavor='complete',mockversion='v3',reg='GCcomb',thetacut=''):
    from pycorr import TwoPointCorrelationFunction
    dirr = '/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/desipipe/'+mockversion+'/'+flavor+'/baseline_2pt/mock'
    xil = []
    sepl = []

    for i in range(0,nmock):
        fn = dirr + str(i)+'/'+rec+'/xi/smu/allcounts_'+tp+'_'+reg+'_z'+zr+'_d4_poles.txt'
        result = np.loadtxt(fn).transpose()
        sep, xis = result[0],result[2:5]
        xil.append(xis)
        sepl.append(sep)
    xi = sum(xil)/nmock
    sep = sum(sepl)/nmock
    return sep,xi

def get_xi_desipipe_ab_baseline(mockn,zr='0.4-0.6',tp='LRG',rec='recon_recsym',nmock=25,flavor='complete',mockversion='v3',reg='GCcomb',thetacut=''):
    from pycorr import TwoPointCorrelationFunction
    dirr = '/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/desipipe/'+mockversion+'/'+flavor+'/baseline_2pt/mock'
    fn = dirr + str(mockn)+'/'+rec+'/xi/smu/allcounts_'+tp+'_'+reg+'_z'+zr+'_d4_poles.txt'
    result = np.loadtxt(fn).transpose()
    return result[0],result[2:5]

def compchi2stats(tracer,zmin,zmax,smin,smax,rec='',ells=[0,2,4],thetacut=''):
    indmin = smin//4
    indmax = smax//4
    zr = str(zmin)+'-'+str(zmax)
    twa = ''
    if tracer == 'ELG_LOP':
        twa = 'notqso'
    xiave,cov = get_xi_cov_desipipe_baseline_txt(zr=zr,tp=tracer,smin=smin,smax=smax,rec=rec,ells=ells,thetacut=thetacut)
    sep,xiave_abamtl = get_xiave_desipipe_ab_baseline(zr=zr,tp=tracer+twa,rec=rec,flavor='altmtl',mockversion='v3_1',thetacut=thetacut)
    _,xiave_abffa = get_xiave_desipipe_ab_baseline(zr=zr,tp=tracer,rec=rec,flavor='ffa',mockversion='v3',thetacut=thetacut)
    xiave_abamtl_cut =  cat_ells(xiave_abamtl[:len(ells)],indrange=[indmin,indmax])
    xiave_abffa_cut =  cat_ells(xiave_abffa[:len(ells)],indrange=[indmin,indmax])
    icov = np.linalg.inv(cov)
    chi2la = []
    chi2lf = []
    for i in range(0,25):
        _,xi_amtl = get_xi_desipipe_ab_baseline(i,zr=zr,tp=tracer+twa,rec=rec,flavor='altmtl',mockversion='v3_1',thetacut=thetacut)
        xi_amtl_cut = cat_ells(xi_amtl[:len(ells)],indrange=[indmin,indmax])
        damtl = xi_amtl_cut-xiave_abamtl_cut
        chi2_amtl = np.dot(damtl,np.dot(damtl,icov))
        chi2la.append(chi2_amtl)
        _,xi_ffa = get_xi_desipipe_ab_baseline(i,zr=zr,tp=tracer,rec=rec,flavor='ffa',mockversion='v3',thetacut=thetacut)
        xi_ffa_cut = cat_ells(xi_ffa[:len(ells)],indrange=[indmin,indmax])
        dffa = xi_ffa_cut-xiave_abffa_cut
        chi2_ffa = np.dot(dffa,np.dot(dffa,icov))
        chi2lf.append(chi2_ffa)
        #print(i,xi_amtl)
        #print(chi2_ffa,chi2_amtl)
    meana = np.mean(chi2la)
    meanf = np.mean(chi2lf)
    fig = plt.figure()
    a = plt.hist(chi2lf,histtype='step',label=r'ffa,$\bar{\chi}^2=$'+str(round(meanf,3)),lw=3,color='b')
    b = plt.hist(chi2la,histtype='step',label=r'altml,$\bar{\chi}^2=$'+str(round(meana,3)),lw=3,color='r')
    plt.plot([meana,meana],[0,max(max(a[0]),max(b[0]))],'r:')
    
    plt.plot([meanf,meanf],[0,max(max(a[0]),max(b[0]))],'b:')

    titl = tracer+' '+str(zmin)+'<z<'+str(zmax)+' for '+str(smin)+'<s<'+str(smax)+' multipoles '+str(ells)
    if rec != '':
        titl += ' recsym'
    if thetacut != '':
        titl += ' thetacut'
    plt.title(titl)
    plt.ylabel('number of mocks')
    plt.xlabel(r'$\chi^2$')
    plt.plot([len(xiave),len(xiave)],[0,max(max(a[0]),max(b[0]))],'k:')
    plt.legend()
    #plt.show()
    return fig#chi2la,chi2lf

figs = []

smin=50
smax=150
recl = ['recon_recsym','']
ellsl = [[0,2,4],[0,2],[0]]

tp = 'QSO'
zrl = [(0.8,2.1)]
for rec in recl:
    for zr in zrl:
        for ells in ellsl:
            fig = compchi2stats(tp,zr[0],zr[1],smin,smax,rec=rec,ells=ells)
            figs.append(fig)

outdir = '/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/'
with PdfPages(outdir+'testEZmockcov_'+tp+'_smin'+str(smin)+'smax'+str(smax)+'.pdf') as pdf:
    for fig in figs:
        pdf.savefig(fig)
        plt.close()

figs = []
tp = 'ELG_LOP'
zrl = [(0.8,1.1),(1.1,1.6)]
for rec in recl:
    for zr in zrl:
        for ells in ellsl:
            fig = compchi2stats(tp,zr[0],zr[1],smin,smax,rec=rec,ells=ells)
            figs.append(fig)
with PdfPages(outdir+'testEZmockcov_'+tp+'_smin'+str(smin)+'smax'+str(smax)+'.pdf') as pdf:
    for fig in figs:
        pdf.savefig(fig)
        plt.close()
            
figs = []
tp = 'LRG'
zrl = [(0.4,0.6),(0.6,0.8),(0.8,1.1)]
for rec in recl:
    for zr in zrl:
        for ells in ellsl:
            fig = compchi2stats(tp,zr[0],zr[1],smin,smax,rec=rec,ells=ells)
            figs.append(fig)

with PdfPages(outdir+'testEZmockcov_'+tp+'_smin'+str(smin)+'smax'+str(smax)+'.pdf') as pdf:
    for fig in figs:
        pdf.savefig(fig)
        plt.close()

smin=20
smax=200

recl = ['']
thetacut='_thetacut0.05'
ellsl = [[0,2,4],[0,2],[0]]

figs = []
tp = 'QSO'
zrl = [(0.8,2.1)]
for rec in recl:
    for zr in zrl:
        for ells in ellsl:
            fig = compchi2stats(tp,zr[0],zr[1],smin,smax,rec=rec,ells=ells,thetacut=thetacut)
            figs.append(fig)

outdir = '/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/'
with PdfPages(outdir+'testEZmockcov_'+tp+thetacut+'_smin'+str(smin)+'smax'+str(smax)+'.pdf') as pdf:
    for fig in figs:
        pdf.savefig(fig)
        plt.close()

figs = []
tp = 'ELG_LOP'
zrl = [(0.8,1.1),(1.1,1.6)]
for rec in recl:
    for zr in zrl:
        for ells in ellsl:
            fig = compchi2stats(tp,zr[0],zr[1],smin,smax,rec=rec,ells=ells,thetacut=thetacut)
            figs.append(fig)
with PdfPages(outdir+'testEZmockcov_'+tp+thetacut+'_smin'+str(smin)+'smax'+str(smax)+'.pdf') as pdf:
    for fig in figs:
        pdf.savefig(fig)
        plt.close()
            
figs = []
tp = 'LRG'
zrl = [(0.4,0.6),(0.6,0.8),(0.8,1.1)]
for rec in recl:
    for zr in zrl:
        for ells in ellsl:
            fig = compchi2stats(tp,zr[0],zr[1],smin,smax,rec=rec,ells=ells,thetacut=thetacut)
            figs.append(fig)

with PdfPages(outdir+'testEZmockcov_'+thetacut+tp+'_smin'+str(smin)+'smax'+str(smax)+'.pdf') as pdf:
    for fig in figs:
        pdf.savefig(fig)
        plt.close()

thetacut=''

recl = ['recon_recsym','']
ellsl = [[0,2,4],[0,2],[0]]

tp = 'QSO'
zrl = [(0.8,2.1)]
for rec in recl:
    for zr in zrl:
        for ells in ellsl:
            fig = compchi2stats(tp,zr[0],zr[1],smin,smax,rec=rec,ells=ells)
            figs.append(fig)

outdir = '/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/'
with PdfPages(outdir+'testEZmockcov_'+tp+'_smin'+str(smin)+'smax'+str(smax)+'.pdf') as pdf:
    for fig in figs:
        pdf.savefig(fig)
        plt.close()

figs = []
tp = 'ELG_LOP'
zrl = [(0.8,1.1),(1.1,1.6)]
for rec in recl:
    for zr in zrl:
        for ells in ellsl:
            fig = compchi2stats(tp,zr[0],zr[1],smin,smax,rec=rec,ells=ells)
            figs.append(fig)
with PdfPages(outdir+'testEZmockcov_'+tp+'_smin'+str(smin)+'smax'+str(smax)+'.pdf') as pdf:
    for fig in figs:
        pdf.savefig(fig)
        plt.close()
            
figs = []
tp = 'LRG'
zrl = [(0.4,0.6),(0.6,0.8),(0.8,1.1)]
for rec in recl:
    for zr in zrl:
        for ells in ellsl:
            fig = compchi2stats(tp,zr[0],zr[1],smin,smax,rec=rec,ells=ells)
            figs.append(fig)

with PdfPages(outdir+'testEZmockcov_'+tp+'_smin'+str(smin)+'smax'+str(smax)+'.pdf') as pdf:
    for fig in figs:
        pdf.savefig(fig)
        plt.close()

