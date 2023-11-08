#standard python
import sys
import os
import shutil
import unittest
from datetime import datetime
import json
import numpy as np
import fitsio
import glob
import argparse
from astropy.table import Table,join,unique,vstack
from matplotlib import pyplot as plt
import LSS.main.cattools as ct
import LSS.common_tools as common

from LSS.globals import main

if os.environ['NERSC_HOST'] == 'cori':
    scratch = 'CSCRATCH'
elif os.environ['NERSC_HOST'] == 'perlmutter':
    scratch = 'PSCRATCH'
else:
    print('NERSC_HOST is not cori or permutter but is '+os.environ['NERSC_HOST'])
    sys.exit('NERSC_HOST not known (code only works on NERSC), not proceeding') 


parser = argparse.ArgumentParser()
parser.add_argument("--realization",type=int)
parser.add_argument("--tracer", help="tracer type to be selected")
parser.add_argument("--base_dir", help="base directory for input/output",default='/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/')
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA",default='Y1')
parser.add_argument("--verspec",help="version for redshifts",default='iron')
parser.add_argument("--data_version",help="version for redshifts",default='v0.6')
parser.add_argument("--minr", help="minimum number for random files",default=0)
parser.add_argument("--maxr", help="maximum for random files, 18 are available (use parallel script for all)",default=18) 
parser.add_argument("--prepsysnet",help="prepare data to get sysnet weights for imaging systematics?",default='n')
parser.add_argument("--add_sysnet",help="add sysnet weights for imaging systematics to full files?",default='n')
parser.add_argument("--imsys_zbin",help="if yes, do imaging systematic regressions in z bins",default='y')
parser.add_argument("--regressis",help="RF weights for imaging systematics?",default='n')
parser.add_argument("--add_regressis",help="add RF weights for imaging systematics?",default='n')
parser.add_argument("--add_regressis_ran",help="add RF weights to randoms?",default='n')
parser.add_argument("--add_sysnet_ran",help="add sysnet weights to randoms",default='n')
parser.add_argument("--par",help="whether to run in parallel",default='n')

parser.add_argument("--add_regressis_ext",help="add RF weights for imaging systematics, calculated elsewhere",default='n')
parser.add_argument("--imsys_nside",help="healpix nside used for imaging systematic regressions",default=256,type=int)
parser.add_argument("--imsys_colname",help="column name for fiducial imaging systematics weight, if there is one (array of ones by default)",default=None)


args = parser.parse_args()
print(args)

tp = args.tracer
rm = int(args.minr)
rx = int(args.maxr)


datadir = '/global/cfs/cdirs/desi/survey/catalogs/'+args.survey+'/LSS/'+args.verspec+'/LSScats/'+args.data_version+'/'
    

if tp[:3] == 'BGS' or tp == 'bright' or tp == 'MWS_ANY':
    prog = 'BRIGHT'

else:
    prog = 'DARK'

progl = prog.lower()

mainp = main(args.tracer,args.verspec,survey=args.survey)

zmin = mainp.zmin
zmax = mainp.zmax


fit_maps = mainp.fit_maps


tpstr = tp
if tp == 'BGS_BRIGHT-21.5':
    tpstr = 'BGS_BRIGHT'
nside = 256


if tp[:3] == 'ELG':
    if args.imsys_zbin == 'y':
        zrl = [(0.8,1.1),(1.1,1.6)]
    else:
        zrl = [(0.8,1.6)]
if tp[:3] == 'QSO':
    if args.imsys_zbin == 'y':
        zrl = [(0.8,1.3),(1.3,2.1),(2.1,3.5)] 
    else:
        zrl = [(0.8,3.5)]   
if tp[:3] == 'LRG':
    if args.imsys_zbin == 'y':
        zrl = [(0.4,0.6),(0.6,0.8),(0.8,1.1)] 
    else:
        zrl = [(0.4,1.1)]  
if tp == 'BGS_BRIGHT-21.5':
    zrl = [(0.1,0.4)]
elif tp[:3] == 'BGS':
    zrl = [(0.01,0.5)]
    zmin = 0.01
    zmax = 0.5    

mockdir = args.base_dir+'mock'+str(args.realization)+'/'

if args.prepsysnet == 'y' or args.regressis == 'y':
    
    def make_hp(value, hpix, nside, fill_with=np.nan):
        """ A Function to create a HEALPix map
        """
        m_ = np.zeros(12*nside*nside)
        m_[:] = fill_with
        m_[hpix] = value
    
        return m_

    import healpy as hp

    dirmap = '/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/v0/kp3_maps/'
    nside = 256#64
    nest = False
    eclrs = ['gr','rz']
    debv = Table()
    for ec in eclrs:
        ebvn = fitsio.read(dirmap+'v0_desi_ebv_'+ec+'_'+str(nside)+'.fits')
        debv_a = ebvn['EBV_DESI_'+ec.upper()]-ebvn['EBV_SFD']
        debv_a = hp.reorder(debv_a,r2n=True)
        debv['EBV_DIFF_'+ec.upper()] = debv_a

dirout = mockdir
lssmapdirout = datadir+'/hpmaps/'
if args.prepsysnet == 'y':
    #logf.write('preparing data to run sysnet regression for '+tp+' '+str(datetime.now())+'\n')
    if not os.path.exists(dirout+'/sysnet'):
        os.mkdir(dirout+'/sysnet')
        print('made '+dirout+'/sysnet')    

    from LSS.imaging import sysnet_tools
    #_HPmapcut'
    datn = fitsio.read(os.path.join(dirout, tp+'_NGC'+'_clustering.dat.fits'))
    dats = fitsio.read(os.path.join(dirout, tp+'_SGC'+'_clustering.dat.fits'))
    dat = np.concatenate((datn,dats))
    dat = common.addNS(Table(dat))
    ranl = []
    for i in range(0,18):
        rann = fitsio.read(os.path.join(dirout, tp+'_NGC'+'_'+str(i)+'_clustering.ran.fits'), columns=['RA', 'DEC','WEIGHT']) 
        rans = fitsio.read(os.path.join(dirout, tp+'_SGC'+'_'+str(i)+'_clustering.ran.fits'), columns=['RA', 'DEC','WEIGHT']) 
        ran = np.concatenate((rann,rans))
        ran = common.addNS(Table(ran))
        ranl.append(ran)
    rands = np.concatenate(ranl)
    regl = ['N','S']
    
    for zl in zrl:
        zw = ''
        if args.imsys_zbin == 'y':
            zw = str(zl[0])+'_'+str(zl[1])
        for reg in regl:
            pwf = lssmapdirout+'QSO_mapprops_healpix_nested_nside'+str(nside)+'_'+reg+'.fits'
            sys_tab = Table.read(pwf)
            cols = list(sys_tab.dtype.names)
            for col in cols:
                if 'DEPTH' in col:
                    bnd = col.split('_')[-1]
                    sys_tab[col] *= 10**(-0.4*common.ext_coeff[bnd]*sys_tab['EBV'])
            for ec in ['GR','RZ']:
                if 'EBV_DIFF_'+ec in fit_maps: 
                    sys_tab['EBV_DIFF_'+ec] = debv['EBV_DIFF_'+ec]

            seld = dat['PHOTSYS'] == reg
            selr = rands['PHOTSYS'] == reg
        
            prep_table = sysnet_tools.prep4sysnet(dat[seld], rands[selr], sys_tab, zcolumn='Z', zmin=zl[0], zmax=zl[1], nran_exp=None,
                    nside=nside, nest=True, use_obiwan=False, columns=fit_maps,wtmd='wt',tp=tp[:3])
            fnout = dirout+'/sysnet/prep_'+tp+zw+'_'+reg+'.fits'
            common.write_LSS(prep_table,fnout)

if args.regressis == 'y':
    #logf.write('adding regressis weights to data catalogs for '+tp+' '+str(datetime.now())+'\n')
    #from regressis, must be installed
    from regressis import DESIFootprint,DR9Footprint

    from LSS.imaging import regressis_tools as rt
    dirreg = dirout+'/regressis_data'
    


    if not os.path.exists(dirreg):
        os.mkdir(dirreg)
        print('made '+dirreg)   
    #pwf = '/global/cfs/cdirs/desi/survey/catalogs/pixweight_maps_all/pixweight-1-dark.fits'   
    sgf = '/global/cfs/cdirs/desi/survey/catalogs/extra_regressis_maps/sagittarius_stream_'+str(nside)+'.npy' 
    dr9_footprint = DR9Footprint(nside, mask_lmc=False, clear_south=True, mask_around_des=False, cut_desi=False)
    desi_footprint = DESIFootprint(nside)
    suffix_tracer = ''
    suffix_regressor = ''

    param = dict()
    param['data_dir'] = dirreg
    param['output_dir'] = dirreg
    param['use_median'] = False
    param['use_new_norm'] = False
    if tp[:3] == 'QSO':
        param['regions'] = ['North', 'South', 'Des']
    else:
        param['regions'] = ['North', 'South_mid_ngc', 'South_mid_sgc']
    max_plot_cart = 300

    cut_fracarea = False
    seed = 42

    #logf.write('using fit maps '+str(fit_maps)+'\n')
    feature_names_ext=None
    pw_out_fn_root = dirout+'/regressis_data/'+tp+'feature_data_'
    regl = ['N','S']
    
    for reg in regl:
        pwf = lssmapdirout+'QSO_mapprops_healpix_nested_nside'+str(nside)+'_'+reg+'.fits'
        sys_tab = Table.read(pwf)
        cols = list(sys_tab.dtype.names)
        for col in cols:
            if 'DEPTH' in col:
                bnd = col.split('_')[-1]
                sys_tab[col] *= 10**(-0.4*common.ext_coeff[bnd]*sys_tab['EBV'])
        for ec in ['GR','RZ']:
            if 'EBV_DIFF_'+ec in fit_maps: 
                sys_tab['EBV_DIFF_'+ec] = debv['EBV_DIFF_'+ec]
        pw_out_fn = pw_out_fn_root+reg+'.fits'
    
        print(pw_out_fn)
        sys_tab.write(pw_out_fn,overwrite=True,format='fits')

        
    use_sgr=False
    if 'SGR' in fit_maps:
        use_sgr = True
        fit_maps.remove('SGR')

    for zl in zrl:    
        zw = str(zl[0])+'_'+str(zl[1])
        print('computing RF regressis weight for '+tp+zw)
        #logf.write('computing RF regressis weight for '+tracer_clus+zw+'\n')
        rt.get_desi_data_clus_compute_weight(dirout, 'main', tp, nside, dirreg, zl, param,foot=dr9_footprint,nran=18,\
        suffix_tracer=suffix_tracer, suffix_regressor=suffix_regressor, cut_fracarea=cut_fracarea, seed=seed,\
         max_plot_cart=max_plot_cart,pixweight_path=pw_out_fn_root,pixmap_external=debv,sgr_stream_path=sgf,\
         feature_names=fit_maps,use_sgr=use_sgr,feature_names_ext=feature_names_ext)
        #rt._compute_weight('main', tracer_clus+zw, dr9_footprint, suffix_tracer, suffix_regressor, cut_fracarea, seed, max_plot_cart,pixweight_path=pw_out_fn,pixmap_external=debv,sgr_stream_path=sgf,feature_names=fit_maps,use_sgr=use_sgr,feature_names_ext=feature_names_ext)

if args.add_regressis == 'y':
    from LSS.imaging import densvar
    from regressis import PhotoWeight
    fb = dirout+tp
    regl = ['NGC','SGC']
    for reg in regl:
        fcd = fb+'_'+reg+'_clustering.dat.fits'
        dd = Table.read(fcd)
        dd['WEIGHT_RF'] = np.ones(len(dd))

        for zl in zrl:    
            print(zl)
            zw = str(zl[0])+'_'+str(zl[1])

            fnreg = dirout+'/regressis_data/main_'+tp+zw+'_256/RF/main_'+tp+zw+'_imaging_weight_256.npy'
            rfpw = PhotoWeight.load(fnreg)
            #print(np.mean(rfpw))
            #dth,dphi = densvar.radec2thphi(dd['RA'],dd['DEC'])
            #dpix = densvar.hp.ang2pix(densvar.nside,dth,dphi,nest=densvar.nest)
            #drfw = rfpw[dpix]
        
            selz = dd['Z'] > zl[0]
            selz &= dd['Z'] <= zl[1]
            dd['WEIGHT_RF'][selz] = rfpw(dd['RA'][selz], dd['DEC'][selz], normalize_map=True)#drfw[selz]
            #norm = 
            print(np.mean(dd['WEIGHT_RF'][selz]))
    #logf.write('added RF regressis weight for '+tracer_clus+zw+'\n')

        common.write_LSS(dd,fcd)#,comments)

regl = ['NGC','SGC']



if args.add_sysnet == 'y':
    #logf.write('adding sysnet weights to data catalogs for '+tp+' '+str(datetime.now())+'\n')
    from LSS.imaging import densvar
    import healpy as hp
    fn_full = dirout+tracer+'_full'+args.use_map_veto+'.dat.fits'
    dd = Table.read(fn_full)
    dd['WEIGHT_SN'] = np.ones(len(dd))
    dth,dphi = densvar.radec2thphi(dd['RA'],dd['DEC'])
    dpix = hp.ang2pix(256,dth,dphi)

    regl_sysnet = ['N','S']
    for reg in regl_sysnet:
        for zl in zrl:
            zw = ''
            if args.imsys_zbin == 'y':
                zw = str(zl[0])+'_'+str(zl[1])
            sn_weights = fitsio.read(dirout+'/sysnet/'+tracer+zw+'_'+reg+'/nn-weights.fits')
            pred_counts = np.mean(sn_weights['weight'],axis=1)
            pix_weight = np.mean(pred_counts)/pred_counts
            pix_weight = np.clip(pix_weight,0.5,2.)
            sn_pix = sn_weights['hpix']
            hpmap = np.ones(12*256*256)
            for pix,wt in zip(sn_pix,pix_weight):
                hpmap[pix] = wt
        
            sel = dd['PHOTSYS'] == reg
            selz = dd['Z_not4clus'] > zl[0]
            selz &= dd['Z_not4clus'] <= zl[1]

            #print(np.sum(sel))
            dd['WEIGHT_SN'][sel&selz] = hpmap[dpix[sel&selz]]
    #print(np.min(dd['WEIGHT_SYS']),np.max(dd['WEIGHT_SYS']),np.std(dd['WEIGHT_SYS']))
    comments = []
    comments.append("Using sysnet for WEIGHT_SYS")

    common.write_LSS(dd,fn_full,comments)

if args.add_regressis_ran == 'y' or args.add_sysnet_ran == 'y':
    if args.add_regressis_ran == 'y':
        wtcol = 'WEIGHT_RF'
    else:
        wtcol = 'WEIGHT_SN'
    fb = dirout+tp
    fcdn = fitsio.read(fb+'_NGC_clustering.dat.fits',columns=['TARGETID',wtcol])
    fcds = fitsio.read(fb+'_SGC_clustering.dat.fits',columns=['TARGETID',wtcol])
    indata = Table(np.concatenate((fcdn,fcds)))
    indata.rename_column('TARGETID', 'TARGETID_DATA')

    def addrancol(rn):
        for reg in regl:
            fname = dirout+tp+'_'+reg+'_'+str(rn)+'_clustering.ran.fits'
            cd = Table(fitsio.read(fname))
            cols2rem = [wtcol,wtcol+'_1',wtcol+'_2']
            for col in cols2rem:
                if col in list(cd.dtype.names):
                    cd.remove_column(col)
            cd = join(cd,indata,keys=['TARGETID_DATA'],join_type='left')
            common.write_LSS(cd,fname)
    
    if args.par == 'n':
        for rn in range(rm,rx):
            addrancol(rn)
    if args.par == 'y':
        nproc = 9
        nran = rx-rm
        inds = np.arange(nran)
        from multiprocessing import Pool
        with Pool(processes=nproc) as pool:
            res = pool.map(addrancol, inds)


# if args.add_sysnet_ran == 'y':
#     fb = dirout+tp
#     fcdn = fitsio.read(fb+'_NGC_clustering.dat.fits',columns=['TARGETID','WEIGHT_SN'])
#     fcds = fitsio.read(fb+'_SGC_clustering.dat.fits',columns=['TARGETID','WEIGHT_SN'])
#     fcd = Table(np.concatenate((fcdn,fcds)))
#     indata.rename_column('TARGETID', 'TARGETID_DATA')
#     
#     regl = ['NGC','SGC']
#     for rn in range(rm,rx):
#         for reg in regl:
#             fname = dirout+tp+'_'+reg+'_'+str(rn)+'_clustering.ran.fits'
#             cd = fitsio.read(fname)
#             cd = join(cd,indata,keys=['TARGETID_DATA'],join_type='left')
#             common.write_LSS(cd,fname)
