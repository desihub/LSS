#sys.path.append(os.getenv('HOME')+'/Y1KPplots/py/')
#Need to get pack in path, via, e.g., PYTHONPATH=$PYTHONPATH:$HOME/Y1KPplots/py/
import numpy as np
from matplotlib import pyplot as plt
import sys
import fitsio
import os
import KP3style


outdir = '/global/cfs/cdirs/desi/survey/catalogs/Y1/KP3plots/KeyPaper/'

axlabsize = 14

zrl = ['0.4_0.6','0.6_0.8','0.8_1.1']
for zr in zrl:
    KP3style.plot_poles_points_werr('LRG',zr)
    mockxi = []
    for i in range(0,20):
        dirmockxi = '/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/mock'+str(i)+'/xi/smu/'
        d = np.loadtxt(dirmockxi+'xipoles_LRG_complete_gtlimaging_GCcomb_'+zr+'_default_lin4_njack0_nran4_split20.txt').transpose()
        mockxi.append(d)
    xiave = mockxi[0]
    for i in range(1,20):
        xiave += mockxi[i]
    xiave /= 20.
    clr = KP3style.colors['LRG'+zr]
    ells = [0,2,4]
    for ell in ells:
        lt = KP3style.ltype[ell]
        av = KP3style.aval[ell]
        plt.plot(xiave[0],xiave[0]**2.*xiave[2+ell//2],lt,color=clr,alpha=av)

#KP3style.plot_poles_points_werr('LRG','0.6_0.8')
#KP3style.plot_poles_points_werr('LRG','0.8_1.1')
plt.ylim(-140,100)
plt.xlabel(r'$s$ ($h^{-1}$Mpc) ',size=axlabsize )
plt.ylabel(r'$s^2\xi_{\ell}$ LRG',size=axlabsize,labelpad=-2)
#plt.title('LRG')
plt.grid()
plt.savefig(outdir+'LRGxi.png')
plt.clf()

rpcut = '_rpcut2.5'
binfac = 5
nmock = 10
for zr in zrl:
    KP3style.plot_pk_poles_points('LRG',zr,rpcut=rpcut,binfac=binfac)
    mockpk = []
    for i in range(0,nmock):
        dirmockpk = '/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/mock'+str(i)+'/pk/'
        d = np.loadtxt(dirmockpk+'pkpoles_LRG_complete_gtlimaging_GCcomb_'+zr+'_default_lin'+str(binfac)+rpcut+'.txt',dtype=complex).transpose()
        mockpk.append(d)
    pkave = mockpk[0]
    for i in range(1,nmock):
        pkave += mockpk[i]
    pkave /= nmock
    clr = KP3style.colors['LRG'+zr]
    ells = [0,2,4]
    for ell in ells:
        lt = KP3style.ltype[ell]
        av = KP3style.aval[ell]
        plt.plot(pkave[1],pkave[1]*pkave[3+ell//2].real,lt,color=clr,alpha=av)

#KP3style.plot_pk_poles_points('LRG','0.6_0.8')
#KP3style.plot_pk_poles_points('LRG','0.8_1.1')
plt.xlabel(r'$k$ ($h$Mpc$^{-1}$) ',size=axlabsize )
plt.ylabel(r'$kP$ LRG',size=axlabsize,labelpad=-2)
plt.grid()
plt.savefig(outdir+'LRGpk'+rpcut+'.png')
plt.clf()

zrl = ['0.8_1.1','1.1_1.6']
for zr in zrl:
    KP3style.plot_poles_points_werr('ELG_LOPnotqso',zr)
    mockxi = []
    for i in range(0,20):
        dirmockxi = '/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/mock'+str(i)+'/xi/smu/'
        d = np.loadtxt(dirmockxi+'xipoles_ELG_LOP_complete_gtlimaging_GCcomb_'+zr+'_default_lin4_njack0_nran4_split20.txt').transpose()
        mockxi.append(d)
    xiave = mockxi[0]
    for i in range(1,20):
        xiave += mockxi[i]
    xiave /= 20.
    clr = KP3style.colors['ELG_LOPnotqso'+zr]
    ells = [0,2,4]
    for ell in ells:
        lt = KP3style.ltype[ell]
        av = KP3style.aval[ell]
        plt.plot(xiave[0],xiave[0]**2.*xiave[2+ell//2],lt,color=clr,alpha=av)

#KP3style.plot_poles_points_werr('ELG_LOPnotqso','0.8_1.1',xi_version='v0.6/blinded/',wt='default_FKP')
#KP3style.plot_poles_points_werr('ELG_LOPnotqso','1.1_1.6',xi_version='v0.6/blinded/',wt='default_FKP')
plt.ylim(-60,40)
plt.xlabel(r'$s$ ($h^{-1}$Mpc) ',size=axlabsize)
plt.ylabel(r'$s^2\xi_{\ell}$ ELG',size=axlabsize)
plt.grid()
plt.savefig(outdir+'ELGxi.png')
plt.clf()

for zr in zrl:
    KP3style.plot_pk_poles_points('ELG_LOPnotqso',zr,rpcut=rpcut,binfac=binfac)
    mockpk = []
    for i in range(0,nmock):
        dirmockpk = '/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/mock'+str(i)+'/pk/'
        d = np.loadtxt(dirmockpk+'pkpoles_ELG_LOP_complete_gtlimaging_GCcomb_'+zr+'_default_lin'+str(binfac)+rpcut+'.txt',dtype=complex).transpose()
        mockpk.append(d)
    pkave = mockpk[0]
    for i in range(1,nmock):
        pkave += mockpk[i]
    pkave /= nmock
    clr = KP3style.colors['ELG_LOPnotqso'+zr]
    ells = [0,2,4]
    for ell in ells:
        lt = KP3style.ltype[ell]
        av = KP3style.aval[ell]
        plt.plot(pkave[1],pkave[1]*pkave[3+ell//2].real,lt,color=clr,alpha=av)


#KP3style.plot_pk_poles_points('ELG_LOPnotqso','0.8_1.1')
#KP3style.plot_pk_poles_points('ELG_LOPnotqso','1.1_1.6')
plt.xlabel(r'$k$ ($h$Mpc$^{-1}$) ',size=axlabsize )
plt.ylabel(r'$kP$ ELG',size=axlabsize,labelpad=-2)
plt.grid()
plt.savefig(outdir+'ELGpk'+rpcut+'.png')
plt.clf()


KP3style.plot_poles_points_werr('QSO','0.8_2.1')
zr = '0.8_2.1'
KP3style.plot_poles_points_werr('QSO',zr)
mockxi = []
for i in range(0,20):
    dirmockxi = '/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/mock'+str(i)+'/xi/smu/'
    d = np.loadtxt(dirmockxi+'xipoles_QSO_complete_gtlimaging_GCcomb_'+zr+'_default_lin4_njack0_nran4_split20.txt').transpose()
    mockxi.append(d)
xiave = mockxi[0]
for i in range(1,20):
    xiave += mockxi[i]
xiave /= 20.
clr = KP3style.colors['QSO'+zr]
ells = [0,2,4]
for ell in ells:
    lt = KP3style.ltype[ell]
    av = KP3style.aval[ell]
    plt.plot(xiave[0],xiave[0]**2.*xiave[2+ell//2],lt,color=clr,alpha=av)

plt.ylim(-70,70)
plt.xlabel(r'$s$ ($h^{-1}$Mpc) ',size=axlabsize)
plt.ylabel(r'$s^2\xi_{\ell}$ QSO',size=axlabsize)
plt.grid()
plt.savefig(outdir+'QSOxi.png')
plt.clf()

KP3style.plot_pk_poles_points('QSO','0.8_2.1',rpcut=rpcut,binfac=binfac)
mockpk = []

for i in range(0,nmock):
    dirmockpk = '/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/mock'+str(i)+'/pk/'
    d = np.loadtxt(dirmockpk+'pkpoles_QSO_complete_gtlimaging_GCcomb_'+zr+'_default_lin'+str(binfac)+rpcut+'.txt',dtype=complex).transpose()
    mockpk.append(d)
pkave = mockpk[0]
for i in range(1,nmock):
    pkave += mockpk[i]
pkave /= nmock
clr = KP3style.colors['QSO'+zr]
ells = [0,2,4]
for ell in ells:
    lt = KP3style.ltype[ell]
    av = KP3style.aval[ell]
    plt.plot(pkave[1],pkave[1]*pkave[3+ell//2].real,lt,color=clr,alpha=av)

plt.xlabel(r'$k$ ($h$Mpc$^{-1}$) ',size=axlabsize )
plt.ylabel(r'$kP$ QSO',size=axlabsize,labelpad=-2)
plt.grid()
plt.savefig(outdir+'QSOpk'+rpcut+'.png')
plt.clf()


KP3style.plot_poles_points_werr('BGS_BRIGHT-21.5','0.1_0.4')
plt.ylim(-100,100)
plt.xlabel(r'$s$ ($h^{-1}$Mpc) ',size=axlabsize)
plt.ylabel(r'$s^2\xi_{\ell}$ BGS',size=axlabsize,labelpad=-2)
plt.grid()
plt.savefig(outdir+'BGSxi.png')
plt.clf()

KP3style.plot_pk_poles_points('BGS_BRIGHT-21.5','0.1_0.4')
plt.xlabel(r'$k$ ($h$Mpc$^{-1}$) ',size=axlabsize )
plt.ylabel(r'$kP$ BGS',size=axlabsize,labelpad=-2)
plt.grid()
plt.savefig(outdir+'BGSpk.png')
plt.clf()


#make legend
tps = ['LRG0.4_0.6','LRG0.6_0.8','LRG0.8_1.1','ELG_LOPnotqso0.8_1.1','ELG_LOPnotqso1.1_1.6','QSO0.8_2.1','BGS_BRIGHT-21.50.1_0.4']
xl = [-2,-1]
yl = [-2,-1]
a = np.ones(10)*-1
plt.figure(figsize=(6,1))
for tp in tps:
    if tp == 'LRG0.4_0.6':
        lab = 'LRG, 0.4<z<0.6'
    if tp == 'LRG0.6_0.8':
        lab = 'LRG, 0.6<z<0.8'
    if tp == 'LRG0.8_1.1':
        lab = 'LRG, 0.8<z<1.1'
    if tp == 'ELG_LOPnotqso0.8_1.1':
        lab = 'ELG, 0.8<z<1.1'
    if tp == 'ELG_LOPnotqso1.1_1.6':
        lab = 'ELG, 1.1<z<1.6'
    if tp == 'QSO0.8_2.1':
        lab = 'QSO, 0.8<z<2.1'
    if tp == 'BGS_BRIGHT-21.50.1_0.4':
        lab = 'BGS, 0.1<z<0.4'
    plt.hist(a,color=KP3style.colors[tp],label=lab)
    #plt.plot(xl,yl,'*',color=KP3style.colors[tp],label=lab)
poles=[0,2,4]
for pole in poles:
    plt.plot(xl,yl,'k'+KP3style.ptype[pole],alpha=KP3style.aval[pole],label=r'$\ell=$'+str(pole))
plt.legend(ncol=3,loc='center')
plt.xlim(0,1)
plt.ylim(0,1)
plt.axis('off')
plt.savefig(outdir+'ell_legend.png')



