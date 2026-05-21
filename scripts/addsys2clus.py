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

import logging
logname = 'mkCat'
logger = logging.getLogger(logname)
logger.setLevel(logging.INFO)

# create console handler and set level to debug
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)

# create formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# add formatter to ch
ch.setFormatter(formatter)

# add ch to logger
logger.addHandler(ch)


#logger = logging.getLogger('mkCat')
#logger.setLevel(level=logging.INFO)

#from desihub
#from desitarget import targetmask
#from regressis, must be installed
#from regressis import DR9Footprint
#sys.path.append('../py') #this requires running from LSS/bin, *something* must allow linking without this but is not present in code yet

#from this package
#try:
import LSS.main.cattools as ct
import LSS.common_tools as common

from LSS.globals import main
#except:
#    print('import of LSS.mkCat_singletile.cattools failed')
#    print('are you in LSS/bin?, if not, that is probably why the import failed')   

if os.environ['NERSC_HOST'] == 'cori':
    scratch = 'CSCRATCH'
elif os.environ['NERSC_HOST'] == 'perlmutter':
    scratch = 'PSCRATCH'
else:
    print('NERSC_HOST is not cori or permutter but is '+os.environ['NERSC_HOST'])
    sys.exit('NERSC_HOST not known (code only works on NERSC), not proceeding') 


parser = argparse.ArgumentParser()
parser.add_argument("--type", help="tracer type to be selected")
parser.add_argument("--basedir", help="base directory for output, default is SCRATCH",default=os.environ[scratch])
parser.add_argument("--version", help="catalog version; use 'test' unless you know what you are doing!",default='test')
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA",default='main')
parser.add_argument("--verspec",help="version for redshifts",default='loa-v1')
parser.add_argument("--ext",help="h5 or fits",default='fits')

parser.add_argument("--extra_clus_dir", help="an optional extra layer of directory structure for clustering catalog",default='')

parser.add_argument("--zcmb", help="whether or not to correct redshifts based on cmb dipole",default='n')
parser.add_argument("--clusran", help="make the random clustering files; these are cut to a small subset of columns",default='n')
parser.add_argument("--relax_zbounds", help="whether or not to use less restricted redshift bounds",default='y')
parser.add_argument("--minr", help="minimum number for random files",default=0,type=int)
parser.add_argument("--maxr", help="maximum for random files, 18 are available (use parallel script for all)",default=18,type=int) 

parser.add_argument("--prepsysnet",help="prepare data to get sysnet weights for imaging systematics?",default='n')
parser.add_argument("--addsysnet",help="add sysnet weights for imaging systematics to full files?",default='n')
parser.add_argument("--imsys_zbin",help="if yes, do imaging systematic regressions in z bins",default='y')

parser.add_argument("--doimlin",help="add weights for imaging systematics using eboss method?",default='n')
parser.add_argument("--imsys_clus_fb",help="perform linear weight fits in fine redshift bins",default='n')
parser.add_argument("--replace_syscol",help="whether to replace any existing weight_sys with new",default='n')
parser.add_argument("--nran4imsys",help="number of random files to using for linear regression",default=10,type=int)

parser.add_argument("--imsys_colname",help="column name for fiducial imaging systematics weight, if there is one (array of ones by default)",default=None)
parser.add_argument("--use_allsky_rands", help="if yes, use all sky randoms to get fractional area per pixel for SYSNet data preparation",default='n')
parser.add_argument("--par", help="if yes, use multiprocessing",default='y')

args = parser.parse_args()
common.printlog(str(args),logger)

type = args.type
tp = type
basedir = args.basedir
version = args.version
specrel = args.verspec
rm = int(args.minr)
rx = int(args.maxr)


common.printlog('running catalogs for tracer type '+type,logger)

    

if type[:3] == 'BGS' or type == 'bright' or type == 'MWS_ANY':
    prog = 'BRIGHT'

elif type[:3] == 'LGE':
    prog = 'DARK1B'

else:
    prog = 'DARK'

progl = prog.lower()

mainp = main(args.type,args.verspec,survey=args.survey,relax_zbounds=args.relax_zbounds)
mdir = mainp.mdir+progl+'/' #location of ledgers
tdir = mainp.tdir+progl+'/' #location of targets
#mtld = mainp.mtld
tiles = mainp.tiles
imbits = mainp.imbits #mask bits applied to targeting
ebits = mainp.ebits #extra mask bits we think should be applied

tsnrcut = mainp.tsnrcut
dchi2 = mainp.dchi2
tsnrcol = mainp.tsnrcol        
zmin = mainp.zmin
zmax = mainp.zmax

#leaving these like this for easier comparison to the other scripts
maindir = basedir +'/'+args.survey+'/LSS/'

ldirspec = maindir+specrel+'/'

dirout = ldirspec+'LSScats/'+version+'/'+args.extra_clus_dir
lssmapdirout = ldirspec+'LSScats/'+version+'/hpmaps'

writefunc = common.write_LSS_scratchcp
if args.ext == 'h5':
    writefunc = write_LSShdf5_scratchcp

def read_file(fn,columns=None):
    if '.fits' in fn:
        data = Table(fitsio.read(fn.replace('global','dvs_ro')))
        if columns is not None:
            data.keep_columns(columns)
    if '.h5' in fn:
        data = common.read_hdf5_blosc(fn.replace('global','dvs_ro'),columns=columns)
    return data
nside=256
tracer_clus = args.type
if args.doimlin == 'y' or args.prep4sysnet == 'y' or args.addsysnet=='y':
    syscol = 'WEIGHT_IMLIN'
    tpstr = args.type
    if "BGS" in tracer_clus:
        tpstr = "BGS_BRIGHT"
        
    if "LRG" in tracer_clus:
        tpstr = "LRG"
    tpmap = tpstr
    inds = np.arange(rm, rx)
    print(tracer_clus[:3])
    if tracer_clus[:3] == "ELG":
        if args.imsys_zbin == "split":
            zrl = [(0.8, 1.1), (1.1, 1.6)]
        elif args.imsys_zbin == 'one':
            zrl = [(0.8, 1.6)]
        zsysmin = 0.8
        zsysmax = 1.6
    
    
    if tracer_clus[:3] == "QSO":
        if args.imsys_zbin == "split":
            zrl = [(0.8, 1.3), (1.3, 2.1), (2.1, 3.5)]
        elif args.imsys_zbin == 'one':
            zrl = [(0.8, 3.5)]
        zsysmin = 0.8
        zsysmax = 3.5
    if tracer_clus[:3] == "LRG":
        if args.imsys_zbin == "split":
            zrl = [(0.4, 0.6), (0.6, 0.8), (0.8, 1.1)]
        elif args.imsys_zbin == 'one':
            zrl = [(0.4, 1.1)]
        zsysmin = 0.4
        zsysmax = 1.1
        #if args.relax_zbounds:
        #    zsysmax = 1.2
        #    zsysmin = 0.3
    
    if "BGS_BRIGHT-" in tracer_clus:
        zrl = [(0.1, 0.4)]
    elif tracer_clus[:3] == "BGS":
        zrl = [(0.01, 0.5)]
        zmin = 0.01
        zmax = 0.5

    
    if args.imsys_zbin == 'fine':
        dz = 0.1
        zm = zsysmin
        #zx = zm + dz
        #redshift_ranges = [(zm, zx)]
        redshift_ranges = []
        while zm < zsysmax:
            zx = zm + dz
            zx = round(zx, 1)
            redshift_ranges += [(zm, zx)]
            zm = zx
        
    else:
        redshift_ranges = zrl
    common.printlog('the redshift bins that will be fit are '+str(redshift_ranges),logger)
    fit_maps = mainp.fit_maps_allebv
    use_maps = fit_maps
    debv = common.get_debv()
    zcmb = common.mk_zcmbmap()


    from LSS.imaging.systematics_linear_regression import (
        make_fit_maps_dictionary,
        produce_imweights,
    )
        # define the paths for the input files
    fname_ngc_out = os.path.join(
        dirout, f"{tracer_clus}_NGC_clustering.dat."+args.ext
    )

    fname_sgc_out = os.path.join(
        dirout, f"{tracer_clus}_SGC_clustering.dat."+args.ext
    )

    # get paths for random catalogs
    randoms_fnames_out = [
        os.path.join(
            dirout,
            f"{tracer_clus}_NGC_{i}_clustering.ran."+args.ext,
        )
        for i in range(args.nran4imsys)
    ] + [
        os.path.join(
            dirout,
            f"{tracer_clus}_SGC_{i}_clustering.ran."+args.ext,
        )
        for i in range(args.nran4imsys)
    ]

    fname_sgc_in = fname_sgc_out
    fname_ngc_in = fname_ngc_out
    randoms_fnames_in = [
        randoms_fname for randoms_fname in randoms_fnames_out
    ]
    
    # Load the data and randoms
    # Get all columns since they will be used for writing later
    data_sgc = read_file(fname_sgc_in)
    data_ngc = read_file(fname_ngc_in)

    data_catalogs = vstack([data_sgc, data_ngc])#np.concatenate([data_sgc, data_ngc])
    #common.printlog(str(np.unique(data_catalogs['PHOTSYS'],return_counts=True)),logger)

    #randoms_catalogs = np.concatenate(
    randoms_catalogs = vstack(
        [read_file(fname) for fname in randoms_fnames_in]
    )
    #common.printlog(str(np.unique(randoms_catalogs['PHOTSYS'],return_counts=True)),logger)
    
if args.doimlin == 'y':
     # perform regression
    weights = produce_imweights(
        data_catalogs=data_catalogs,
        randoms_catalogs=randoms_catalogs,
        is_clustering_catalog=True,
        weight_scheme=None,
        tracer_type=tracer_clus,
        redshift_range=redshift_ranges,
        templates_maps_path_S=os.path.join(
            lssmapdirout, f"{tpmap}_mapprops_healpix_nested_nside{nside}_S.fits"
        ),
        templates_maps_path_N=os.path.join(
            lssmapdirout, f"{tpmap}_mapprops_healpix_nested_nside{nside}_N.fits"
        ),
        fit_maps=fit_maps,
        output_directory=dirout,
        output_catalog_path=None,  # writing to disk will be done later to handle SGC/NGC separately
        output_column_name=syscol,
        save_summary_plots=True,
        nbins=10,  # is the default
        tail=0.5,  # is the default
        logger=logger,
        loglevel="INFO",
    )
    
    # Data catalogs are already loaded, just recast them to astropy Tables
    data_sgc = Table(data_sgc)
    data_ngc = Table(data_ngc)
    # Catalogs are just concatenated in the order SGC, NGC
    # so this is enough to assign weights to the correct one
    transition_index = len(data_sgc)
    assert transition_index + len(data_ngc) == len(weights), "Shape mismatch!"
    # add custom column to catalog
    data_sgc[syscol] = weights[:transition_index]
    data_ngc[syscol] = weights[transition_index:]
    # overwrite the WEIGHT columns
    if args.replace_syscol:
        data_sgc["WEIGHT"] /= data_sgc["WEIGHT_SYS"]
        data_sgc["WEIGHT_SYS"] = data_sgc[syscol]
        data_sgc["WEIGHT"] *= data_sgc["WEIGHT_SYS"]

        data_ngc["WEIGHT"] /= data_ngc["WEIGHT_SYS"]
        data_ngc["WEIGHT_SYS"] = data_ngc[syscol]
        data_ngc["WEIGHT"] *= data_ngc["WEIGHT_SYS"]
    # write out everything
    writefunc(
        data_sgc,
        fname_sgc_out,
        logger=logger,
    )
    writefunc(
        data_ngc,
        fname_ngc_out,
        logger=logger,
    )

    #  also write the weights in the randoms
    #if args.imsys_clus_ran:
    fname = os.path.join(
        dirout,  f"{tracer_clus}_NGC_clustering.dat."+args.ext
    )
    dat_ngc = Table(read_file(fname, columns=["TARGETID", syscol]))
    fname = os.path.join(
        dirout, f"{tracer_clus}_SGC_clustering.dat."+args.ext
    )
    dat_sgc = Table(read_file(fname, columns=["TARGETID", syscol]))
    dat = vstack([dat_sgc, dat_ngc])
    dat.rename_column("TARGETID", "TARGETID_DATA")
    regl = ["NGC", "SGC"]
    syscolr = syscol

    # if args.replace_syscol == 'y':
    #    syscolr = 'WEIGHT_SYS'
    def _add2ran(rann):
        for reg in regl:
            ran_fn = os.path.join(
                dirout,
                f"{tracer_clus}_{reg}_{rann}_clustering.ran."+args.ext,
            )
            ran = Table(read_file(ran_fn))
            if syscolr in ran.colnames:
                ran.remove_column(syscolr)
            ran = join(ran, dat, keys=["TARGETID_DATA"])
            if args.replace_syscol:
                ran["WEIGHT"] /= ran["WEIGHT_SYS"]
                ran["WEIGHT_SYS"] = ran[syscolr]
                ran["WEIGHT"] *= ran["WEIGHT_SYS"]
            writefunc(ran, ran_fn, logger=logger)

    if args.par == "y":
        from multiprocessing import Pool

        with Pool() as pool:
            res = pool.map(_add2ran, inds)
    else:
        for rn in inds:  # range(rm,rx):
            _add2ran(rn)

if args.prepsysnet == 'y':
    common.printlog('preparing data to run sysnet regression for '+tracer_clus,logger)
    if not os.path.exists(dirout+'/sysnet'):
        os.mkdir(dirout+'/sysnet')
        print('made '+dirout+'/sysnet')    

    from LSS.imaging import sysnet_tools
    
    regl = ['N','S']
    
    for zl in zrl:
        zw = ''
        zmin,zmax=zl[0],zl[1]
        #if args.imsys_zbin == 'y':
        zw = str(zmin)+'_'+str(zmax)
        for reg in regl:
            if type == 'LRG':
                if reg == 'N':
                    fitmapsbin = fit_maps
                else:
                    if zmax == 0.6:
                        fitmapsbin = mainp.fit_maps46s
                    if zmax == 0.8:
                        fitmapsbin = mainp.fit_maps68s
                    if zmax == 1.1:
                        fitmapsbin = mainp.fit_maps81s
            else:
                fitmapsbin = fit_maps
            pwf = lssmapdirout+'/'+tpmap+'_mapprops_healpix_nested_nside'+str(nside)+'_'+reg+'.fits'
            sys_tab = Table.read(pwf)
            cols = list(sys_tab.dtype.names)
            for col in cols:
                if 'DEPTH' in col:
                    bnd = col.split('_')[-1]
                    sys_tab[col] *= 10**(-0.4*common.ext_coeff[bnd]*sys_tab['EBV'])
            for ec in ['GR','RZ']:
                if 'EBV_DIFF_'+ec in fit_maps: 
                    sys_tab['EBV_DIFF_'+ec] = debv['EBV_DIFF_'+ec]
            if 'EBV_DIFF_MPF' in fit_maps:
                sys_tab['EBV_DIFF_MPF'] = sys_tab['EBV'] - sys_tab['EBV_MPF_Mean_FW15']
            if 'ZCMB' in fit_maps:
                sys_tab['ZCMB'] = zcmb
            seld = data_catalogs['PHOTSYS'] == reg
            selr = randoms_catalogs['PHOTSYS'] == reg
            #if args.use_allsky_rands == 'y':
            allsky_fn = f"/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/allsky_rpix_{reg}_nran18_nside256_ring.fits"
            allsky_rands = fitsio.read(allsky_fn)
            allrands = allsky_rands['RANDS_HPIX'] # randoms count per hp pixel
            #    selr_all = allsky_rands['PHOTSYS'] == reg
            #    allrands = allsky_rands[selr_all]
            #else:
            #    allrands = None
            common.printlog(f"{tpstr} {reg} z{zmin}-{zmax}: {fitmapsbin}",logger)
            wtmd = 'fracz'
            common.printlog('using '+tpmap +' maps and '+wtmd+' weights')
            prep_table = sysnet_tools.prep4sysnet(data_catalogs[seld], randoms_catalogs[selr], sys_tab, zcolumn='Z', allsky_rands=allrands, 
                                                  zmin=zl[0], zmax=zl[1], nran_exp=None, nside=nside, nest=True, use_obiwan=False,
                                                  columns=fitmapsbin,wtmd=wtmd)
            fnout = dirout+'/sysnet/prep_'+tracer_clus+zw+'_'+reg+'.fits'
            if not os.path.isdir(dirout+'/sysnet/'):
                os.makedirs( dirout+'/sysnet/')
            common.write_LSS_scratchcp(prep_table,fnout,logger=logger)


if args.addsysnet == 'y':
    common.printlog('adding sysnet weights to data catalogs for '+tracer_clus,logger)
    from LSS.imaging import densvar
    import healpy as hp
    #fn_full = dirout+tracer_clus+'_full'+args.use_map_veto+'.dat.fits'
    #dd = Table.read(fn_full)
    data_catalogs['WEIGHT_SN'] = np.ones(len(data_catalogs))
    dth,dphi = densvar.radec2thphi(data_catalogs['RA'],data_catalogs['DEC'])
    dpix = hp.ang2pix(256,dth,dphi)

    regl_sysnet = ['N','S']
    for reg in regl_sysnet:
        for zl in zrl:
            #zw = ''
            #if args.imsys_zbin == 'y':
            zw = str(zl[0])+'_'+str(zl[1])
            sn_weights = fitsio.read(dirout+'/sysnet/'+tracer_clus+zw+'_'+reg+'/nn-weights.fits')
            pred_counts = np.mean(sn_weights['weight'],axis=1)
            #pix_weight = np.mean(pred_counts)/pred_counts
            #pix_weight = np.clip(pix_weight,0.5,2.)
            pix_weight = 1./pred_counts
            pix_weight = pix_weight / pix_weight.mean()
            pix_weight = np.clip(pix_weight,0.5,2.)
            sn_pix = sn_weights['hpix']
            hpmap = np.ones(12*256*256)
            for pix,wt in zip(sn_pix,pix_weight):
                hpmap[pix] = wt
        
            sel = data_catalogs['PHOTSYS'] == reg
            selz = data_catalogs['Z'] > zl[0]
            selz &= data_catalogs['Z'] <= zl[1]

            #print(np.sum(sel))
            data_catalogs['WEIGHT_SN'][sel&selz] = hpmap[dpix[sel&selz]]
    # Catalogs are just concatenated in the order SGC, NGC
    # so this is enough to assign weights to the correct one
    transition_index = len(data_sgc)
    assert transition_index + len(data_ngc) == len(data_catalogs['WEIGHT_SN']), "Shape mismatch!"
    # add custom column to catalog
    data_sgc['WEIGHT_SN'] = data_catalogs['WEIGHT_SN'][:transition_index]
    data_ngc['WEIGHT_SN'] = data_catalogs['WEIGHT_SN'][transition_index:]
    syscol = 'WEIGHT_SN'
    # overwrite the WEIGHT columns
    if args.replace_syscol:
        data_sgc["WEIGHT"] /= data_sgc["WEIGHT_SYS"]
        data_sgc["WEIGHT_SYS"] = data_sgc[syscol]
        data_sgc["WEIGHT"] *= data_sgc["WEIGHT_SYS"]

        data_ngc["WEIGHT"] /= data_ngc["WEIGHT_SYS"]
        data_ngc["WEIGHT_SYS"] = data_ngc[syscol]
        data_ngc["WEIGHT"] *= data_ngc["WEIGHT_SYS"]
    # write out everything
    writefunc(
        data_sgc,
        fname_sgc_out,
        logger=logger,
    )
    writefunc(
        data_ngc,
        fname_ngc_out,
        logger=logger,
    )

    #  also write the weights in the randoms
    #if args.imsys_clus_ran:
    
    fname = os.path.join(
        dirout,  f"{tracer_clus}_NGC_clustering.dat."+args.ext
    )
    dat_ngc = Table(read_file(fname, columns=["TARGETID", syscol]))
    fname = os.path.join(
        dirout, f"{tracer_clus}_SGC_clustering.dat."+args.ext
    )
    dat_sgc = Table(read_file(fname, columns=["TARGETID", syscol]))
    dat = vstack([dat_sgc, dat_ngc])
    dat.rename_column("TARGETID", "TARGETID_DATA")
    regl = ["NGC", "SGC"]
    syscolr = syscol

    # if args.replace_syscol == 'y':
    #    syscolr = 'WEIGHT_SYS'
    def _add2ran(rann):
        for reg in regl:
            ran_fn = os.path.join(
                dirout,
                f"{tracer_clus}_{reg}_{rann}_clustering.ran."+args.ext,
            )
            ran = Table(read_file(ran_fn))
            if syscolr in ran.colnames:
                ran.remove_column(syscolr)
            ran = join(ran, dat, keys=["TARGETID_DATA"])
            if args.replace_syscol:
                ran["WEIGHT"] /= ran["WEIGHT_SYS"]
                ran["WEIGHT_SYS"] = ran[syscolr]
                ran["WEIGHT"] *= ran["WEIGHT_SYS"]
            writefunc(ran, ran_fn, logger=logger)

    if args.par == "y":
        from multiprocessing import Pool

        with Pool() as pool:
            res = pool.map(_add2ran, inds)
    else:
        for rn in inds:  # range(rm,rx):
            _add2ran(rn)
