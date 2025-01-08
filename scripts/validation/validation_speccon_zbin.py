import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import os
import sys
import argparse

import fitsio
from astropy.table import join,Table

import LSS.common_tools as common


parser = argparse.ArgumentParser()
parser.add_argument("--version", help="catalog version",default='test')
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA",default='Y1')
parser.add_argument("--tracers", help="all runs all for given survey",default='all')
parser.add_argument("--verspec",help="version for redshifts",default='iron')
parser.add_argument("--data",help="LSS or mock directory",default='LSS')
parser.add_argument("--hpmapcut",help="string to add for hpmapcut or not",default='_HPmapcut')
args = parser.parse_args()


indir = '/dvs_ro/cfs/cdirs/desi/survey/catalogs/'+args.survey+'/'+args.data+'/'+args.verspec+'/LSScats/'+args.version+'/'
specdir =  '/dvs_ro/cfs/cdirs/desi/survey/catalogs/'+args.survey+'/'+args.data+'/'+args.verspec+'/'
outdir = indir.replace('dvs_ro','global')+'plots/speccon/'

if args.data == 'LSS':
    if not os.path.exists(outdir):
        os.mkdir(outdir)
        print('made '+outdir)


zcol = 'Z'
nran = 18

tps = [args.tracers]
if args.tracers == 'all':
    tps = ['QSO','LRG','BGS_BRIGHT','ELG_LOPnotqso']

zdw = ''#'zdone'

regl = ['_N','_S']

if args.survey == 'SV3' and args.tracers == 'all':
    tps = ['QSO','LRG','BGS_ANY','BGS_BRIGHT','ELG','ELG_HIP','ELG_HIPnotqso','ELGnotqso']
    zdw = ''
    if args.data != 'LSS':
        tps = ['QSO','LRG','ELG']


def normed_plot(data,seltot,sel,per_range=(.5,99.5),col='TSNR2_ELG',wcol='WEIGHT_ZFAIL',cl='k',ps='o',lab=''):
    selnz = data[col] > 0 #all quantities have true values > 0
    if col == 'TEMPAIRmPMIRROR':
        selnz = data[col] != -99
    datastnz = data[seltot&selnz]
    range = (np.percentile(datastnz[col],per_range[0]),np.percentile(datastnz[col],per_range[1]))
    a,bins = np.histogram(data[seltot][col],range=range)
    datas = data[sel]
    b,_ = np.histogram(datas[col],bins=bins)
    c,_ = np.histogram(datas[col],bins=bins,weights=datas[wcol])
    sp = bins[1]-bins[0]
    err = np.sqrt(b*(1-b/a))/a  # binomial error
    normc = np.sum(c)/np.sum(a)
    cchi = np.sum((c/a/normc-1.)**2./(err/normc)**2.)
    plt.errorbar(bins[:-1]+sp/2.,c/a/normc,err/normc,fmt=cl+ps,label=lab+' weighted '+str(round(cchi,3)))
    normb = np.sum(b)/np.sum(a)
    bchi = np.sum((b/a/normb-1.)**2./(err/normb)**2.)
    plt.plot(bins[:-1]+sp/2.,b/a/normb,cl+'--',label=lab+' unweighted '+str(round(bchi,3)))



for tp in tps:
    prog='dark'
    if tp[:3] == 'QSO':
        zmin = 0.8
        zmax = 3.5
        dz = 0.3
        cols = columns=['TARGETID','TILEID','LOCATION','Z_not4clus','WEIGHT_ZFAIL','PHOTSYS','PRIORITY','ZWARN']
    if tp[:3] == 'ELG':
        zmin = 0.6
        zmax = 1.6
        dz = 0.1
        cols = columns=['TARGETID','TILEID','LOCATION','Z_not4clus','WEIGHT_ZFAIL','PHOTSYS','ZWARN','o2c','DELTACHI2']
    if tp[:3] == 'LRG':
        zmin = 0.4
        zmax = 1.1
        dz = 0.1
        cols = columns=['TARGETID','TILEID','LOCATION','Z_not4clus','WEIGHT_ZFAIL','PHOTSYS','ZWARN','DELTACHI2']
    if tp[:3] == 'BGS':
        dz = 0.1
        zmin = 0.1
        zmax = 0.5
        prog='bright'
        cols = columns=['TARGETID','TILEID','LOCATION','Z_not4clus','WEIGHT_ZFAIL','PHOTSYS','ZWARN','DELTACHI2']
    df = fitsio.read(indir+tp+zdw+'_full'+args.hpmapcut+'.dat.fits',columns=cols)
    speccon = Table(fitsio.read(specdir+'specobscon_'+prog+'.fits'))
    speccon['TILEID'] = speccon['TILEID'].astype(int)
    speccon['TEMPAIRmPMIRROR'] = -99.
    speccon['SPEED'] = -99.
    sel = speccon['PMIRTEMP'] > 0
    sel &= speccon['TAIRTEMP'] > 0
    speccon['TEMPAIRmPMIRROR'][sel] = speccon['TAIRTEMP'][sel] - speccon['PMIRTEMP'][sel]
    sel = speccon['EXPTIME'] > 0
    sel &= speccon['EFFTIME_SPEC'] > 0
    speccon['SPEED'][sel] = speccon['EFFTIME_SPEC'][sel]/speccon['EXPTIME'][sel]

    df = join(df,speccon,keys=['TARGETID','TILEID','LOCATION'],join_type='left')

    #dcl = []
    #for reg in regl:
    #    dtc = fitsio.read(indir+tp+zdw+reg+'_clustering.dat.fits')
    #    dcl.append(dtc)
    #dc = np.concatenate(dcl)
    #df = join(dtf,dc,keys=['TARGETID'],join_type='left')
    selgz = common.goodz_infull(tp[:3],df)#df['Z'].mask == False
    selo = df['ZWARN'] != 999999
    selo &= df['ZWARN']*0 == 0
    if tp[:3] == 'QSO':
        selo &= df['PRIORITY'] == 3400 #repeats throw things off

    mean_gz = sum(df[selgz]['WEIGHT_ZFAIL'])/len(df[selo])
    print('number with good z, sum of weight_zfail,  number with good obs')
    print(len(df[selgz]),sum(df[selgz]['WEIGHT_ZFAIL']),len(df[selo]))
    fl = ['SPEED','TEMPAIRmPMIRROR','WINDDIR','SEEING_ETC','EBV','SEEING_GFA','SKY_MAG_AB_GFA',\
    'SKY_MAG_G_SPEC','SKY_MAG_R_SPEC','SKY_MAG_Z_SPEC','ETCTRANS','ETCSKY','ETCTHRUB','ZD','TURBRMS','SLEWANGL','AIRMASS','MOON_ILLUM',\
    'TRANSPARENCY_GFA','WINDSPD','HUMIDITY','ACQFWHM','PMIRTEMP','TAIRTEMP','PARALLAC','ROTOFFST']
    for feat in fl:
        df[feat] = df[feat].filled(-99)
        print('doing '+feat)
        zm = zmin
        figs = []
        selz = df['Z_not4clus'] > zm
        selz &= df['Z_not4clus'] < zmax
        seln = df['PHOTSYS'] == 'N'
        fig = plt.figure()
        normed_plot(df,selo&seln,selo&selgz&selz&seln,cl='r',ps='d',lab='N',col=feat)
        normed_plot(df,selo&~seln,selo&selgz&selz&~seln,cl='b',ps='o',lab='S',col=feat)
        plt.legend()
        plt.grid()
        plt.ylabel(tp+' relative z success')
        plt.xlabel(feat)
        plt.title(str(round(zm,3))+'<z<'+str(round(zmax,3)))
        figs.append(fig)
        print(zm)
        while zm < zmax:
            selz = df['Z_not4clus'] > zm
            selz &= df['Z_not4clus'] < zm+dz
            seln = df['PHOTSYS'] == 'N'
            fig = plt.figure()
            normed_plot(df,selo&seln,selo&selgz&selz&seln,cl='r',ps='d',lab='N',col=feat)
            normed_plot(df,selo&~seln,selo&selgz&selz&~seln,cl='b',ps='o',lab='S',col=feat)
            plt.legend()
            plt.grid()
            plt.ylabel(tp+' relative z success')
            plt.xlabel(feat)
            plt.title(str(round(zm,3))+'<z<'+str(round(zm+dz,3)))
            #plt.savefig(outdir+tp+'_'+str(round(zm,3))+'ltzlt'+str(round(zm+dz,3))+'_relsuccess_tnsr.png')
            #plt.clf()
            figs.append(fig)
            zm += round(dz,1)
            print(zm)
        with PdfPages(outdir+tp+'_relsuccess_'+feat+'_zbins.pdf') as pdf:
            for fig in figs:
                pdf.savefig(fig)
                plt.close()
        print('done with '+feat)
    print('done with '+tp)

