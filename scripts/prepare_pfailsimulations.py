import numpy as np
import scipy
import sys
from scipy import stats
from scipy.stats import norm
from scipy.stats import binom
import fitsio
import glob
import os
import matplotlib.pyplot as plt
import statistics
import argparse
import astropy
from astropy.table import Table
from astropy.table import join
from astropy.table import Column
from astropy.time import Time
from astropy.io import fits
from mpi4py import MPI
rank = MPI.COMM_WORLD.Get_rank()
import warnings
warnings.filterwarnings('ignore')

sys.path.insert(0,'/global/cfs/cdirs/desi/users/akrolew/LSS/py/LSS')
from ssr_tools_new import model_ssr, model_ssr_zfac
import LSS.common_tools as common

fibers_per_task = int(sys.argv[6])

# createbins 
def create_bins(lower_bound, width, upper_bound):
    bins = []
    for low in range(lower_bound, upper_bound, width):
        bins.append((low, low+width))
    return bins

def find_bin(value, bins):
    for i in range(0, len(bins)):
        if bins[i][0] <= value < bins[i][1]:
            return i
    return -1


#function to preprocess data from fitsfile

def ash_code(tp):  # Make changes for Y1/daily will run
    from desitarget import targetmask
    if tp =='LRG':
        lb=lower[0]
    elif tp =='BGS_BRIGHT':
        lb = lower[1]
    elif tp =='BGS_FAINT':
        lb=lower[1]
    elif tp =='ELG_LOPnotqso':
        lb = lower[2]
    elif tp =='ELG_VLOnotqso':
        lb = lower[2]
    elif tp =='QSO':
        lb = lower[3]
    
    
    
    if (survey == 'main') or (survey == 'Y1') or (survey == 'DA02') or (survey == 'DA2'):
        if ((specver == 'jura-v1')  or (specver == 'kibo-v1')) and (tp[:3] != 'BGS'):
            zf = basedir+'/'+survey+'/LSS/'+specver+'/datcomb_dark_spec_zdone.fits'
        elif ((specver == 'jura-v1') or (specver == 'kibo-v1')) and (tp[:3] == 'BGS'):
            zf = basedir+'/'+survey+'/LSS/'+specver+'/datcomb_bright_spec_zdone.fits'
        else:
            if tp == 'BGS_FAINT':
                zf = basedir+'/'+survey+'/LSS/'+specver+('/LSScats/%s/' % ver)+'BGS_ANY'+'_full_noveto.dat.fits'
            elif tp == 'ELG_VLOnotqso':
                zf = basedir+'/'+survey+'/LSS/'+specver+('/LSScats/%s/' % ver)+'ELG'+'_full_noveto.dat.fits'
            else:
                zf = basedir+'/'+survey+'/LSS/'+specver+('/LSScats/%s/' % ver)+tp+'_full_noveto.dat.fits'
        dz = Table(fitsio.read(zf))
        if (specver == 'jura-v1') or (specver == 'kibo-v1'):
            if tp[:3] == 'BGS':
                targs = Table(fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/main/LSS/BGS_ANYtargetsDR9v1.1.1.fits'))
            else:
                targs = Table(fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/main/LSS/%stargetsDR9v1.1.1.fits' % tp[:3]))
            dz = join(dz, targs, keys='TARGETID',join_type='left',uniq_col_name='{col_name}{table_name}',table_names=['','_TARG'])
            if tp == 'LRG':
                wtype = ((dz['DESI_TARGET'] & 2**0) != 0)
            elif tp == 'ELG_LOPnotqso':
                wtype = ((dz['DESI_TARGET'] & 2**5) != 0)
                wtype &= ((dz['DESI_TARGET'] & 2**2) == 0)
            elif tp == 'QSO':
                wtype = ((dz['DESI_TARGET'] & 2**2) != 0)
            elif tp == 'BGS_BRIGHT':
            	wtype = ((dz['BGS_TARGET'] & 2**1) != 0)
            dz = dz[wtype]
        if tp == 'ELG_LOPnotqso':
            wtype = ((dz['DESI_TARGET'] & 4) == 0) #remove QSO
            dz = dz[wtype]
        elif tp == 'ELG_VLOnotqso':
            wtype = ((dz['DESI_TARGET'] & 4) == 0) #remove QSO
            dz = dz[wtype]
            wtype = ((dz['DESI_TARGET'] & targetmask.desi_mask['ELG_VLO']) != 0) #keep VLO
            dz = dz[wtype] 
            
        elif tp == 'BGS_FAINT':
            wtype = ((dz['BGS_TARGET'] & targetmask.bgs_mask['BGS_BRIGHT']) == 0) #remove BGS_BRIGHT
            dz = dz[wtype]
        dz = common.cut_specdat(dz)
        from LSS.globals import main
        #pars = main(tp,specver)
        #print('len of pars.elgzf',len(pars.elgzf))
        

    elif survey == 'SV3':
        sys.exit('not written for SV3 yet')
        zf = basedir+'/'+survey+'/LSS/'+specver+'/datcomb_dark_tarspecwdup_Alltiles.fits'
        dz = Table(fitsio.read(zf))
        desitarg = 'SV3_DESI_TARGET'
        bit = 1 #for selecting LRG
        wtype = ((dz[desitarg] & bit) > 0)
        print(len(dz[wtype]))
        #dz = dz[wtype&wg]
        dz = dz[wtype]
        wz = dz['ZWARN'] != 999999 #this is what the null column becomes
        wz &= dz['ZWARN']*0 == 0 #just in case of nans
        wz &= dz['COADD_FIBERSTATUS'] == 0
        ff = dz[wz]

        zf = basedir+'/'+survey+'/LSS/'+specver+'/datcomb_bright_tarspecwdup_Alltiles.fits'
        dz = Table(fitsio.read(zf))
        desitarg = 'SV3_BGS_TARGET'
        wtype = dz[desitarg] > 0#((dz[desitarg] & bit) > 0)
        print(len(dz[wtype]))
        #dz = dz[wtype&wg]
        dz = dz[wtype]
        wz = dz['ZWARN'] != 999999 #this is what the null column becomes
        wz &= dz['ZWARN']*0 == 0 #just in case of nans
        wz &= dz['COADD_FIBERSTATUS'] == 0

        ff2 = dz[wz]

    if (tp == 'ELG_LOPnotqso') or (tp == 'ELG_VLOnotqso'):
        if (specver == 'jura-v1') or (specver == 'kibo-v1'):
            emlin = Table(fitsio.read(basedir+'/'+survey+'/LSS/'+specver+'/emlin_catalog.fits'))
            emlin['TILEID'] = emlin['TILEID'].astype('int')
            dz = join(dz, emlin, keys=['TARGETID','TILEID','LOCATION'],join_type='left',uniq_col_name='{col_name}{table_name}',table_names=['','_EMLIN'])

    if tp == 'QSO':
        if (specver == 'jura-v1') or (specver == 'kibo-v1'):
            qso = Table(fitsio.read(basedir+'/'+survey+'/QSO/%s/QSO_cat_%s_cumulative_v1.fits' % (specver.split('-')[0],specver.split('-')[0])))
            dz = join(dz,qso,keys=['TARGETID','TILEID','LOCATION'],join_type='left',uniq_col_name='{col_name}{table_name}',table_names=['','_QF'])
            dz['Z'].name = 'Z_orig' #rename the original redrock redshifts
            dz['Z_QF'].name = 'Z' #the redshifts from the quasar file should be used instead


    z_tot = dz['ZWARN'] != 999999
    z_tot &= dz['ZWARN']*0 == 0
    z_tot &= ((dz['COADD_FIBERSTATUS'] == 0) | (dz['COADD_FIBERSTATUS'] == 8))
    if tp[:3] == 'BGS':
        z_tot &= dz['TSNR2_BGS'] > 1000
    else:
        z_tot &= dz['TSNR2_ELG'] > 80
    # TSNR cuts
    print('L148',np.unique(dz[z_tot]['COADD_FIBERSTATUS']))
    


    if tp == 'LRG':
        z_suc= dz['ZWARN']==0
        z_suc &= dz['DELTACHI2']>15
        z_suc &= dz['Z']<1.5
       
    if (tp == 'ELG_LOPnotqso') or (tp == 'ELG_VLOnotqso'):
        #else:
        z_suc = np.log10(dz['OII_FLUX'] * np.sqrt(dz['OII_FLUX_IVAR']))+0.2*np.log10(dz['DELTACHI2'])> 0.9

    
    if tp == 'QSO':
        #print(5/0)
        #qsozf = pars.qsozf
        #if specver == 'guadalupe':
        #    #qsozf = '/global/cfs/cdirs/desi/users/edmondc/QSO_catalog/guadalupe/QSO_cat_guadalupe_cumulative.fits'
        #arz = Table(fitsio.read(qsozf))
        #arz.keep_columns(['TARGETID','LOCATION','TILEID','Z','Z_QN'])
        #arz['TILEID'] = arz['TILEID'].astype(int)


        #dz = join(dz,arz,keys=['TARGETID','TILEID','LOCATION'],join_type='left',uniq_col_name='{col_name}{table_name}',table_names=['','_QF'])
        #dz['Z'].name = 'Z_RR' #rename the original redrock redshifts
        #dz['Z_QF'].name = 'Z' #the redshifts from the quasar file should be used instead
        if (specver == 'jura-v1') or (specver == 'kibo-v1'):
            z_suc = (dz['Z'] != 999999.0) | (dz['Z'].mask == False)
        else:
            z_suc = dz['Z'] != 999999.0
    print('L191',np.unique(dz[z_tot]['COADD_FIBERSTATUS']))
            
    ind = np.where(dz['TARGETID'] == 39627380788040729)[0]
    print('dz ind 192',dz[ind])
    #print('z_suc ind',z_suc[ind])
    #print('z_tot ind',z_tot[ind])
    print('zwarn == 0 ind',(dz['ZWARN'] == 0)[ind])
    print('deltachi2 > 40',(dz['DELTACHI2'] > 40)[ind])
    print('deltachi2',dz['DELTACHI2'][ind])

    if tp == 'BGS_BRIGHT':
        z_suc = dz['ZWARN']==0
        z_suc &= dz['DELTACHI2']>40

    ind = np.where(dz['TARGETID'] == 39627380788040729)[0]
    print('dz ind 195',dz[ind])
    print('z_suc ind',z_suc[ind])
    print('z_tot ind',z_tot[ind])
    print('zwarn == 0 ind',(dz['ZWARN'] == 0)[ind])
    print('deltachi2 > 40',(dz['DELTACHI2'] > 40)[ind])

        
    if tp == 'BGS_FAINT':
        z_suc = dz['ZWARN']==0
        z_suc &= dz['DELTACHI2']>40
    
    # QSO: 0.8 < z < 3.5

    # ELG: 0.8 < z < 1.5

    # LRG: 0.4 < z < 1.1

    # BGS 0.1 < z < 0.5
    
    if (specver == 'jura-v1') or (specver == 'kibo-v1'):
        cat = Table(dz)
        #cat['Z_not4clus'] = cat['Z']
    
        selobs = cat['ZWARN']*0 == 0
        selobs &= cat['ZWARN'] != 999999
        if tp[:3] != 'BGS':
            selobs &= cat['TSNR2_ELG'] > 80
        else:
            selobs &= cat['TSNR2_BGS'] > 1000
    
        if tp[:3] == 'LRG':
            band = 'Z'
            mintsnr=500/12.15
            #maxtsnr =2000/12.15
            maxtsnr =1700/12.15
        elif tp[:3] == 'QSO':
            band = 'R'
            mintsnr=450/(8.60/0.255)
            maxtsnr=1800/(8.60/0.255)
        elif tp[:3] == 'BGS':
            band = 'R'
            mintsnr=120/(12.15/89.8)
            maxtsnr =300/(12.15/89.8)
        elif tp[:3] == 'ELG':
            band = 'G'
            mintsnr = 80
            maxtsnr = 200
        if tp[:3] != 'ELG':
            modelN = model_ssr(cat[selobs],tracer=tp[:3],reg='N',band=band,tsnr_min=mintsnr,tsnr_max=maxtsnr,readpars=True,outdir=basedir+'/'+survey+'/LSS/'+specver+'/LSScats/%s/' % ver,outfn_root=tp,overwrite_pars_ssrmaxflux=False)
            modelS = model_ssr(cat[selobs],tracer=tp[:3],reg='S',band=band,tsnr_min=mintsnr,tsnr_max=maxtsnr,readpars=True,outdir=basedir+'/'+survey+'/LSS/'+specver+'/LSScats/%s/' % ver,outfn_root=tp,overwrite_pars_ssrmaxflux=False)
        else:
            cat.add_column(np.log10(cat['OII_FLUX'] * np.sqrt(cat['OII_FLUX_IVAR']))+0.2*np.log10(cat['DELTACHI2']), name='o2c')
            modelN = model_ssr_zfac(cat[selobs],reg='N',outdir=basedir+'/'+survey+'/LSS/'+specver+'/LSScats/test/',outfn_root=tp[:3])
            modelS = model_ssr_zfac(cat[selobs],reg='S',outdir=basedir+'/'+survey+'/LSS/'+specver+'/LSScats/test/',outfn_root=tp[:3])
    
        if tp[:3] == 'LRG':
            modelN.fluxfittype='piecewise'
            modelN.flux_break = 3
            modelS.fluxfittype = 'piecewise'
            modelS.flux_break = 3
        elif tp[:3] != 'ELG':
            modelN.fluxfittype='linear'
            modelS.fluxfittype='linear'
    
        if tp[:3] != 'ELG':
            parsmaxflux = np.loadtxt(basedir+'/'+survey+'/LSS/'+specver+'/LSScats/test/%sNpars_ssrmaxflux.txt' % tp)
            modelN.pars_ferf = parsmaxflux
    
            parsmaxflux = np.loadtxt(basedir+'/'+survey+'/LSS/'+specver+'/LSScats/test/%sSpars_ssrmaxflux.txt' % tp)
            modelS.pars_ferf = parsmaxflux

    
        wtfN, modN = modelN.add_modpre(cat)
        wtfS, modS = modelS.add_modpre(cat)
    
        mod = np.zeros(len(cat))
        mod[cat['PHOTSYS'] == 'S'] = modS[cat['PHOTSYS'] == 'S']
        mod[cat['PHOTSYS'] == 'N'] = modN[cat['PHOTSYS'] == 'N']
    
        wtf = np.zeros(len(cat))
        wtf[cat['PHOTSYS'] == 'S'] = wtfS[cat['PHOTSYS'] == 'S']
        wtf[cat['PHOTSYS'] == 'N'] = wtfN[cat['PHOTSYS'] == 'N']
        
    ind = np.where(dz['TARGETID'] == 39627380788040729)[0]
    print('dz ind 277',dz[ind])
    print('z_suc ind',z_suc[ind])
    print('z_tot ind',z_tot[ind])

    print('L294',np.unique(dz[z_tot]['COADD_FIBERSTATUS']))
    return(dz,z_suc,z_tot,lb,mod)

#fiber cut function
def fiber_cut(cut,tp,dz,z_tot,z_suc): 
    if cut==0: tag='full'
    elif cut==1:
        badfiber_list = np.loadtxt("/global/homes/a/akrolew/lrg+bgs_3sig_bad.txt",dtype=int)
        tag='check'
    elif cut==2:
        mfailpvalues = np.loadtxt("/global/homes/a/akrolew/foldtest/"+survey+"/"+specver+"/"+tp+"_zsuc.txt")[:,6]
        tag = 'mfailp'
        mfail=mfailpvalues<norm.cdf(-4.0)
        badfiber_list = np.loadtxt("/global/homes/a/akrolew/foldtest/"+survey+"/"+specver+"/"+tp+"_zsuc.txt")[:,0]
        badfiber_list = badfiber_list[mfail]
    elif cut==3:
        if tp=='LRG' or 'BGS_ANY':
            badfiber_list = np.loadtxt("/global/homes/a/akrolew/foldtest/"+survey+"/"+specver+"/plotdata/newcut/"+tp+"_newcut.txt",dtype=int)
            tag= 'newcut'
        else:
            tag='uncut'
            badfiber_list = []
    else:
        print("wrong input for cut\n Execution stopped")
        exit()    
    if(cut): 
        badfib_mask=dz==dz
        for f in badfiber_list:
            badfib_mask[dz['FIBER'] == f] = False
        z_suc&=badfib_mask
    
    print("zsuccess rate for "+tp,len(dz[z_suc&z_tot])/len(dz[z_tot])) # success rate with criterion
    
    return(tag, dz, z_tot, z_suc)

#writing files 

def write_in_file(tp,dz,z_tot,z_suc,cut,tag):
    fibl,n_tot = np.unique(dz[z_tot]['FIBER'],return_counts=True)
    fiblg,n_g = np.unique(dz[z_suc&z_tot]['FIBER'],return_counts=True)
    fib_test = np.isin(fibl,fiblg)
    z_tot &= np.isin(dz['FIBER'],fibl[fib_test]) # looks like z tots change here due to isin(ask ashley)
    fibl,n_tot = np.unique(dz[z_tot]['FIBER'],return_counts=True)# n total is based on updated z_tot due to previous line
    # sumz=sum(dz
    if np.array_equal(fibl,fiblg):
        gfrac = n_g/n_tot
    else:
        sys.exit('need to put something in for mismatch fiber lists')
    fn = '/global/homes/a/akrolew/foldtest/'    
    fn=fn+survey+'/'+specver+'/'+tp+'_zsuc.txt'
    sumz=[]
    sumzsq=[]
    mean=np.sum(n_g)/np.sum(n_tot)
    print('mean',mean)
    #print(1/0)
    if tag != 'full':
        txtmfailp =np.loadtxt("/global/homes/a/akrolew/foldtest/"+survey+"/"+specver+"/"+tp+"_zsuc.txt")[:,6]# dont recalculate!!
        txtfib = np.loadtxt("/global/homes/a/akrolew/foldtest/"+survey+"/"+specver+"/"+tp+"_zsuc.txt")[:,0]
    binno=[]
    value=[]
    morefailp=[]
    bins = create_bins(lb,1,ub)
    #tsnr_mean=[]
    #modsucmean=[]
    if tag == 'full':
        fo = open(fn,'w')
        
    for ii in range(len(fibl)):
        m=dz['FIBER']==fibl[ii]
        m&= z_suc
        m&= z_tot
        n = n_tot[ii]
        s = n_g[ii]
        p = mean
        #modsucmean.append(np.mean(dz['mod_success_rate'][m]))
        
        
        
        zvalues= dz['Z'][m&z_tot]
        #if tp == 'ELGnotqso':
            #tsnr_mean.append(np.mean(dz['TSNR2_ELG'][m]))
        #else:
            #exec("tsnr_mean.append( np.mean(dz['TSNR2_"+tp+"'][m]))")
        sumz.append(np.sum(zvalues))
        sumzsq.append(np.sum(zvalues**2))

        if tag != 'full':
            morefailp.append(txtmfailp[txtfib == fibl[ii]])# load it from text files directly!!!
        else: 
            morefailp.append(binom.cdf(s-1, n, p))
            
        value.append(np.log10(morefailp[ii]))
    
        ####cdf = value = morefailp = 1/2 * (1 + erf(x/sqrt(2)))
        #value2 = np.sqrt(2) * scipy.special.erfinv(2 * morefailp[ii] - 1)
        binno.append(find_bin(value[ii],bins))
        
        if tag == 'full':
            fo.write(str(fibl[ii])+' '+str(n_g[ii]/n_tot[ii])+' '+str(n_g[ii])+' '
            +str(n_tot[ii])+'  '+str(sumz[ii])+' '+str(sumzsq[ii])+' '
            +str(morefailp[ii])+' '+str(binno[ii])+'\n')
    if tag == 'full':
        fo.close()
    
    
    dz['bin']=np.zeros(len(dz))
    for i in range(len(fibl)):
        dzfib= dz['FIBER']==fibl[i] 
        dz['bin'][dzfib]= binno[i]

    return (bins, binno, dz, z_tot, z_suc)

#function to create p_values for each fiber

def weighted_pvalues(dz,z_tot,z_suc):
    dz['p_value'] = np.zeros(len(dz))
    #dz['p_suc']= np.zeros(len(dz))
    #dz['mod_suc']=np.zeros(len(dz))
    #dz['p_suc_w']=np.zeros(len(dz))
    p_suc = len(dz[z_suc&z_tot])/len(dz[z_tot])
    #dz = dz[z_tot]
    #norm =  len(dz['WEIGHT_ZFAIL'])/np.sum(1/dz['WEIGHT_ZFAIL'])
    fibl,n_tot = np.unique(dz[z_tot]['FIBER'],return_counts=True)
    fiblg,n_g = np.unique(dz[z_suc&z_tot]['FIBER'],return_counts=True)
    lownt = fibl
    #lownt = np.loadtxt('/global/homes/a/akrolew/foldtest/Y1/himalayas/concat/fibermissing.txt')
    #lownt = np.loadtxt('/global/cfs/cdirs/desi/users/akrolew/QSOmissing.txt')
    p_valuefib = np.zeros(len(lownt))
    Varr = np.zeros(len(lownt))
    Nsig = np.zeros(len(lownt))
    Nlownt = np.zeros(len(lownt))
    Slownt = np.zeros(len(lownt))
    #ff2['obs_fail_p'] = np.zeros(len(ff2))
    #lownt = fibl[np.where(n_tot <100)]
    #print("number of low n tot:", len(lownt))
    start_time = time.time()
    n = np.zeros(len(lownt))
    s = np.zeros(len(lownt))
    mean_x = np.zeros(len(lownt))
    mean_y = np.zeros(len(lownt))
    mean_mod_suc = np.zeros(len(lownt))
    frac_succ = np.zeros(len(lownt))
    f = open('/global/cfs/cdirs/desi/users/akrolew/summary/%s/%i.txt'%(tp,int(sys.argv[5])  + rank),'w')
    
    for j in range(len(lownt)):
        fmask=dz['FIBER']==lownt[j]
 
        n[j] = n_tot[fibl==lownt[j]]
        s[j] = n_g[fiblg==lownt[j]]
        mean_x[j] = np.mean(dz['FIBERASSIGN_X'][fmask])
        mean_y[j] = np.mean(dz['FIBERASSIGN_Y'][fmask])
        mod_suc = dz['mod_success_rate'][fmask & z_tot]
        mean_mod_suc[j]  = np.mean(mod_suc)
        #frac_suc = s/n
        frac_succ[j] = s[j]/n[j] #frac_suc
    
    scale = np.mean(frac_succ)/np.mean(mean_mod_suc)  
    dz['mod_success_rate'] = dz['mod_success_rate']*scale
    dz['mod_success_rate'][dz['mod_success_rate'] > 1] = 1
    dz['mod_success_rate'][dz['mod_success_rate'] < 0] = 0
    mean_mod_suc = mean_mod_suc*scale
    #tot_mod_suc_rate = np.sum(dz['mod_success_rate'][z_tot]) / len(dz['mod_success_rate'][z_tot])
    #print(np.where(lownt == 3986)[0])
    for ii in range(len(lownt)):
        if (ii-int(sys.argv[5]))//fibers_per_task  == rank: #change number from 0 to 7 for each batch run
            #if ii == 661:
            
            fmask=dz['FIBER']==lownt[ii]
            #if fibl[ii] == 1995:
            #    print(5/0)
            #fmask&= z_suc
            fmask&= z_tot
            #n = n_tot[fibl==lownt[ii]]
            #s = n_g[fibl==lownt[ii]]

            #print(s)
            #print(len(dz['WEIGHT_ZFAIL'][fmask]))
            #obs = dz['WEIGHT_ZFAIL'][fmask]

            #testobs = 1/obs * p_suc * norm
            #TSNR2_Sum = np.sum(1/dz['WEIGHT_ZFAIL'][fmask])

            mod_suc = dz['mod_success_rate'][fmask] 
            #print(5/0)
            #frac_suc = s/n
            #frac_suc_w = np.sum(dz['WEIGHT_ZFAIL'][fmask])/n

            sim_n= int(3e6)
            np.random.seed(7) 


            p_succ = np.zeros(sim_n)
            #p_succ = np.sum(sample, axis=1)/len(mod_suc)
            for i in range(sim_n):
                sample = np.random.binomial(1, mod_suc)
                p_succ[i] = np.sum(sample)/len(sample)
            #p_succ = np.fromfile('/global/homes/a/akrolew/foldtest/Y1/himalayas/fibersims/%i.bin' % lownt[ii])
            #print(p_succ, frac_suc)
            print("--- Total time is %s seconds ---" % (time.time() - start_time))

            #s = np.sqrt(1/n*(np.sum(mod_suc-np.mean(mod_suc))**2))
            ss = np.std(mod_suc)
            var = (n[ii]*np.mean(mod_suc)*(1-np.mean(mod_suc))-n[ii]*ss**2)/n[ii]**2 #confirm with alex
            nsig = (frac_succ[ii]-np.mean(mod_suc))/np.sqrt(var)
            p_value = np.searchsorted(np.sort(p_succ), 1e-10+frac_succ[ii])/sim_n #changed it to frac_suc from 1-frac_suc
            dz['p_value'][fmask] = p_value
            p_valuefib[ii] = p_value
            Varr[ii] = var
            Nsig[ii] = nsig
            #Nlownt[ii] = n
            #Slownt[ii] = s

            #np.savetxt('/global/homes/a/akrolew/foldtest/Y1/himalayas/fibersims/%i.txt'%lownt[ii],p_succ)
            p_succ.tofile('/global/cfs/cdirs/desi/users/akrolew/fibersims/%s/%i.bin'%(tp,lownt[ii]))
            #print("2")
            #print(nsig)

            #print("3")
            #dz['p_suc'][fmask] = frac_suc
            #print("4")
            #print(ii,dz['p_suc'][ii])
            #dz['p_suc_w'][fmask] = frac_suc_w

            #dz['testobs'][ii] = np.mean(testobs)
            #print("loop no:",ii)
            f.write('%i  %.5f %.5e %.5f %i %i %.5f %.5f %.2f %.2f\n' % (fibl[ii],p_value, var, nsig, n[ii], s[ii], mean_mod_suc[ii],frac_succ[ii], mean_x[ii], mean_y[ii]))
            print('%i  %.5f %.5e %.5f %i %i %.5f %.5f %.2f %.2f\n' % (fibl[ii],p_value, var, nsig, n[ii], s[ii], mean_mod_suc[ii],frac_succ[ii], mean_x[ii], mean_y[ii]))

            f.flush()
            print("--- Total time is %s seconds ---" % (time.time() - start_time)) 
     
    f.close()
    return(dz,var,Varr,Nsig,p_valuefib,lownt,ss,n,mod_suc,s,Nlownt,Slownt)


basedir='/global/cfs/cdirs/desi/survey/catalogs' 
survey=sys.argv[1]#str(input("Enter survey name:\n")) #DA02 or main
specver=sys.argv[2]#str(input("Enter specver:\n"))
ver = sys.argv[4]

tracers=[sys.argv[3]]
cut = 0 #int(input("Remove the badfibers from the analysis?\n  0)No 1)check 2)mfailp 3) newcut\n"))
lower = [-22,-14,-13,-6]#use as appropriate for version
ub = 0
import time
start_time = time.time()



for tp in tracers:
    dz, z_suc, z_tot, lb,mod= ash_code(tp) #Selects catalogs and cuts them down to Z_suc using criteria for tracers.
    tag, dz, z_tot, z_suc = fiber_cut(cut, tp, dz, z_tot, z_suc)# Makes a bad fiber cut and updates z_suc/z_tot
    bins, binno, dz,z_tot, z_suc = write_in_file(tp, dz, z_tot, z_suc, cut, tag)# creates bins and binnos
    t = Table(np.array([dz['FIBER'],dz['FIBERASSIGN_X'],dz['FIBERASSIGN_Y'],mod]).T,names=('FIBER','FIBERASSIGN_X','FIBERASSIGN_Y','mod_success_rate'))
    t.write('pfailsimulations_%s/%s_cat_trunc.fits' % (specver,tp),overwrite=True)
    np.savetxt('pfailsimulations_%s/%s_ztot.txt' % (specver,tp),z_tot)
    if tp == 'QSO':
        fucker = np.zeros(len(z_suc))
        fucker[z_suc.mask==True] = 0
        fucker[z_suc.mask==False] = 1
        np.savetxt('pfailsimulations_%s/QSO_zsuc.txt' % specver,fucker)

    else:
        np.savetxt('pfailsimulations_%s/%s_zsuc.txt' % (specver,tp),z_suc)
    #weighted_pvalues(dz,z_tot,z_suc)
    #print(dz['p_value'])
    #print('time',time.time()-start_time)
    
    