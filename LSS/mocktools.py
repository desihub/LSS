import fitsio
import numpy as np
from astropy.table import Table,join,vstack
# system
import os
import subprocess
import sys
import tempfile
import shutil
import re
import astropy.io.fits as pf
# time
from time import time
from datetime import datetime, timedelta

# desi
import desitarget
from desitarget.io import read_targets_in_tiles 

import LSS.common_tools as common


def get_fba_mock(mockdir,mocknum,survey='DA02',prog='dark'):
    #produces script to run to get mock fiberassign files
    mock_fn = mockdir+'/forFA'+str(mocknum)+'.fits'
    if not os.path.exists(mockdir+'/'+survey):
        os.mkdir(mockdir+'/'+survey)
        print('made '+mockdir+'/'+survey)
    if not os.path.exists(mockdir+'/'+survey+'/fba'+str(mocknum)):
        os.mkdir(mockdir+'/'+survey+'/fba'+str(mocknum))
        print('made '+mockdir+'/'+survey+'/fba'+str(mocknum))

    tile_fn = '/global/cfs/cdirs/desi/survey/catalogs/'+survey+'/LSS/tiles-'+prog.upper()+'.fits'
    tiles = Table(fitsio.read(tile_fn,columns=['TILEID','RA','DEC']))
    tiles['OBSCONDITIONS'] = 1
    tiles['IN_DESI'] = 1
    tiles['PROGRAM'] = 'MAIN'
    
    ts = str(tiles['TILEID'][0]).zfill(6)
    #get info from origin fiberassign file
    fht = fitsio.read_header('/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'+ts[:3]+'/fiberassign-'+ts+'.fits.gz')
    skyf = '/global/cfs/cdirs/desi/survey/catalogs/'+survey+'/LSS/skies-'+prog.upper()+'.fits'
    outdir = mockdir+'/'+survey+'/fba'+str(mocknum)
    tile_fn =  outdir+'/tiles.fits'
    tiles.write(tile_fn,overwrite=True)
    tars = read_targets_in_tiles(mock_fn,tiles)
    tarfn = outdir+'/targs.fits'
    Table(tars).write(tarfn,format='fits',overwrite=True)

    fo = open(outdir+'/fa-'+ts+'.sh','w')
    fo.write('#!/bin/bash\n\n')
    fo.write('source /global/common/software/desi/desi_environment.sh main\n')
    fo.write("module swap fiberassign/5.0.0\n")

    fo.write("fba_run")
    fo.write(" --targets "+tarfn)
    fo.write(" --sky "+skyf)
    fo.write(" --footprint "+tile_fn)
    rundate= fht['RUNDATE']
    fo.write(" --rundate "+rundate)
    fo.write(" --fieldrot "+str(fht['FIELDROT']))
    fo.write(" --dir "+outdir)
    fo.write(" --sky_per_petal 40 --standards_per_petal 10")
    fo.write(" --sky_per_slitblock 1")
    fo.write(" --ha "+str(fht['FA_HA']))
    fo.write(" --margin-gfa 0.4 --margin-petal 0.4 --margin-pos 0.05")
    fo.close()
    print('wrote scripts for fiberassign '+outdir+'/fa-'+ts+'.sh') 
    return(outdir+'/fa-'+ts+'.sh')   

def get_fba_mock_ran(mockdir,rannum,survey='DA02',prog='dark'):
    #produces script to run to get mock fiberassign files
    from fiberassign.targets import (TargetsAvailable)
    from fiberassign.utils import option_list, GlobalTimers
    from fiberassign.hardware import load_hardware
    from fiberassign.tiles import load_tiles, Tiles
    from fiberassign.targets import (TARGET_TYPE_SCIENCE, TARGET_TYPE_SKY,
                                 TARGET_TYPE_SUPPSKY,
                                 TARGET_TYPE_STANDARD, TARGET_TYPE_SAFE,
                                 Targets, TargetsAvailable,
                                 LocationsAvailable, load_target_file)
    from fiberassign.assign import (Assignment, write_assignment_fits,
                                write_assignment_ascii, merge_results,
                                read_assignment_fits_tile)                                 
 
 
    mock_fn = mockdir+'/ran_forFA'+str(rannum)+'.fits'
    if not os.path.exists(mockdir+'/'+survey):
        os.mkdir(mockdir+'/'+survey)
        print('made '+mockdir+'/'+survey)
    dirout = mockdir+'/'+survey+'/random_fba'+str(rannum)
    if not os.path.exists(dirout):
        os.mkdir(dirout)
        print('made '+dirout)

    if survey == 'MVMY1':
        tile_fn = '/global/cfs/cdirs/desi/users/FA_EZ_1year/fiberassign_EZ_3gpc/fba001/inputs/tiles.fits'
    else:
        tile_fn = '/global/cfs/cdirs/desi/survey/catalogs/'+survey+'/LSS/tiles-'+prog.upper()+'.fits'
    tiles = Table(fitsio.read(tile_fn,columns=['TILEID','RA','DEC']))
    tiles['OBSCONDITIONS'] = 1
    tiles['IN_DESI'] = 1
    tiles['PROGRAM'] = 'MAIN'
    
    ts = str(tiles['TILEID'][0]).zfill(6)
    #get info from origin fiberassign file
    fht = fitsio.read_header('/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'+ts[:3]+'/fiberassign-'+ts+'.fits.gz')
    rundate= fht['RUNDATE']
    tile_fn =  dirout+'/tiles.fits'
    tiles.write(tile_fn,overwrite=True)
    tarfn = dirout+'/targs.fits'
    if os.path.isfile(tarfn) == False:
        tars = read_targets_in_tiles(mock_fn,tiles)
        print(len(tars)) 
        Table(tars).write(tarfn,format='fits',overwrite=True)
        print('wrote '+tarfn)

    
    from fiberassign.targets import TargetTagalong,create_tagalong
    tgs = Targets()
    tagalong = create_tagalong()#TargetTagalong([])
    load_target_file(tgs,tagalong,tarfn)
    print('loaded target file '+tarfn)
    
    hw = load_hardware(rundate=rundate)
    tiles = load_tiles(tiles_file=tile_fn)
    from fiberassign.targets import targets_in_tiles
    tile_targetids, tile_x, tile_y = targets_in_tiles(hw, tgs, tiles,tagalong)
    tgsavail = TargetsAvailable(hw, tiles, tile_targetids, tile_x, tile_y)
    favail = LocationsAvailable(tgsavail)
    asgn = Assignment(tgs, tgsavail, favail,{}) #this is needed for fiberassign 2.4 and higher(?)

    asgn.assign_unused(TARGET_TYPE_SCIENCE)
    write_assignment_fits(tiles,tagalong, asgn, out_dir=dirout, all_targets=True)
    print('wrote assignment files to '+dirout)  

def combtiles_assign_wdup_7pass(indir,outdir,tarf,addcols=['TARGETID','RSDZ','ZWARN','PRIORITY'],fba=True,tp='dark'):

    s = 0
    td = 0
    #tiles.sort('ZDATE')
    
    outf = outdir+'/datcomb_'+tp+'assignwdup.fits'
    if fba:
        pa_hdu = 'FASSIGN'
    tl = []
    for pass_num in range(0,7):
        passdir = indir+'faruns/farun-pass'+str(pass_num)+'/'
        tiles = fitsio.read(passdir+'tiles-pass'+str(pass_num)+'.fits')
        for tile in tiles['TILEID']:
            if fba:
                ffa = passdir+'/fba-'+str(tile).zfill(6)+'.fits'
            if os.path.isfile(ffa):
                fa = Table(fitsio.read(ffa,ext=pa_hdu,columns=['TARGETID','LOCATION']))
                sel = fa['TARGETID'] >= 0
                fa = fa[sel]
                td += 1
                fa['TILEID'] = int(tile)
                tl.append(fa)
                print(td,len(tiles))
            else:
                print('did not find '+ffa)
    dat_comb = vstack(tl)
    print(len(dat_comb))
    tar_in = Table(fitsio.read(tarf))#,columns=addcols))
    cols = list(tar_in.dtype.names)
    if 'ZWARN' not in cols:
        tar_in['ZWARN'] = np.zeros(len(tar_in),dtype=int)
    tar_in.keep_columns(addcols)
    dat_comb = join(dat_comb,tar_in,keys=['TARGETID'])
    print(len(dat_comb))
    
    dat_comb.write(outf,format='fits', overwrite=True)
    print('wrote '+outf)
    return dat_comb

def combtiles_pa_wdup_7pass(indir,outdir,tarf,addcols=['TARGETID','RA','DEC'],fba=True,tp='dark',ran='dat',dtar=''):
    if ran == 'dat':
        #addcols.append('PRIORITY')
        addcols.append('PRIORITY')
        addcols.append(dtar+'DESI_TARGET')
    s = 0
    td = 0
    #tiles.sort('ZDATE')
    #print(len(tiles))
    outf = outdir+'/'+ran+'comb_'+tp+'wdup.fits'
    if fba:
        pa_hdu = 'FAVAIL'
    tl = []
    for pass_num in range(0,7):
        passdir = indir+'faruns/farun-pass'+str(pass_num)+'/'
        tiles = fitsio.read(passdir+'tiles-pass'+str(pass_num)+'.fits')
        for tile in tiles['TILEID']:
            if fba:
                ffa = passdir+'/fba-'+str(tile).zfill(6)+'.fits'
            if os.path.isfile(ffa):
                fa = Table(fitsio.read(ffa,ext=pa_hdu,columns=['TARGETID','LOCATION']))
                sel = fa['TARGETID'] >= 0
                fa = fa[sel]
                td += 1
                fa['TILEID'] = int(tile)
                tl.append(fa)
                print(td,len(tiles))
            else:
                print('did not find '+ffa)
    dat_comb = vstack(tl)
    print(len(dat_comb))
    tar_in = fitsio.read(tarf,columns=addcols)
    dat_comb = join(dat_comb,tar_in,keys=['TARGETID'])
    print(len(dat_comb))
    dat_comb.rename_column('PRIORITY', 'PRIORITY_INIT') 
    if dtar != '':
        dat_comb.rename_column(dtar+'DESI_TARGET', 'DESI_TARGET') 
    dat_comb.write(outf,format='fits', overwrite=True)
    print('wrote '+outf)
    return dat_comb


def mkclusdat_allpot(fl,ztable,tp='',dchi2=9,tsnrcut=80,rcut=None,ntilecut=0,ccut=None,ebits=None,zmin=0,zmax=6):
    '''
    make data clustering for mock with everything in the full catalog
    fl is the root of the input/output file
    weighttileloc determines whether to include 1/FRACZ_TILELOCID as a completeness weight
    zmask determines whether to apply a mask at some given redshift
    tp is the target type
    dchi2 is the threshold for keeping as a good redshift
    tnsrcut determines where to mask based on the tsnr2 value (defined below per tracer)

    '''
    wzm = '_complete_'
    if ccut is not None:
        wzm = ccut+'_' #you could change this to however you want the file names to turn out

    if rcut is not None:
        wzm += 'rmin'+str(rcut[0])+'rmax'+str(rcut[1])+'_'
    if ntilecut > 0:
        wzm += 'ntileg'+str(ntilecut)+'_'
    outf = fl+wzm+'clustering.dat.fits'
    ff = Table.read(fl+'_full.dat.fits')
    cols = list(ff.dtype.names)
    print(len(ff))
    ff = join(ff,ztable,keys=['TARGETID'])
    print('after join to z',str(len(ff)))
    ff['WEIGHT'] = np.ones(len(ff))
    
    kl = ['RA','DEC','Z','WEIGHT','TARGETID']
    wn = ff['PHOTSYS'] == 'N'

    ff.keep_columns(kl)
    print('minimum,maximum weight')
    print(np.min(ff['WEIGHT']),np.max(ff['WEIGHT']))

    #comments = ["DA02 'clustering' LSS catalog for data, all regions","entries are only for data with good redshifts"]
    #common.write_LSS(ff,outf,comments)

    outfn = fl+wzm+'N_clustering.dat.fits'
    comments = ["DA02 'clustering' LSS catalog for data, BASS/MzLS region","entries are only for data with good redshifts"]
    common.write_LSS(ff[wn],outfn,comments)

    outfn = fl+wzm+'S_clustering.dat.fits'
    comments = ["DA02 'clustering' LSS catalog for data, DECaLS region","entries are only for data with good redshifts"]
    ffs = ff[~wn]
    common.write_LSS(ffs,outfn,comments)
    
def mkclusdat_tiles(fl,ztable,bit=None,zmin=0,zmax=6):
    '''
    make data clustering given some input with RA,DEC,Z,DESI_TARGET assuming it is complete (all targets in region have a redshift)
    `fl` (string) is the root of the output file name 
    `ztable` is an input astropy table with at least RA,DEC,Z columns
    `bit` is used if the input includes all tracer types and you want to select a particular one given DESI_TARGET
    `zmin` and `zmax` are floats that apply any redshift bounds to the output catalog
    '''
    wzm = '_tiles_'

    if bit is not None:
        sel = ztable['DESI_TARGET'] & bit > 0
        ff = ztable[sel]
    else:
        ff = ztable
    common.addNS(ff)
    #ff['PHOTSYS'] = 'N'
    #sel = ff['DEC'] < 32.375 #this is imperfect for the SGC (some of it is > 32.375 but still DECaLS), fix in the future
    #ff['PHOTSYS'][sel] = 'S'    
    
    ff['WEIGHT'] = np.ones(len(ff))
    
    kl = ['RA','DEC','Z','WEIGHT','TARGETID']
    wn = ff['PHOTSYS'] == 'N'

    ff.keep_columns(kl)
    print('minimum,maximum weight')
    print(np.min(ff['WEIGHT']),np.max(ff['WEIGHT']))

    #comments = ["DA02 'clustering' LSS catalog for data, all regions","entries are only for data with good redshifts"]
    #common.write_LSS(ff,outf,comments)

    outfn = fl+wzm+'N_clustering.dat.fits'
    #edit these comments at some point
    comments = ["DA02 'clustering' LSS catalog for data, BASS/MzLS region","entries are only for data with good redshifts"]
    common.write_LSS(ff[wn],outfn,comments)

    outfn = fl+wzm+'S_clustering.dat.fits'
    comments = ["DA02 'clustering' LSS catalog for data, DECaLS region","entries are only for data with good redshifts"]
    ffs = ff[~wn]
    common.write_LSS(ffs,outfn,comments)
    
def mkclusran_tiles(ffc,fl,rann,rcols=['Z','WEIGHT']):
    '''
    `ffc` is an input astropy table with at least columns RA,DEC, with RA,DEC assumed to be randomly sampling the area
    associated with the data file
    `fl` is a string that points to the data catalog file format and is used for the random catalog format
    `rann` is the number associated with the random file
    `rcols` is the list of columns to sample from the data catalog
    '''
    
    wzm = ''
    fcdn = Table.read(fl+wzm+'N_clustering.dat.fits')
    kc = ['RA','DEC','Z','WEIGHT']
    rcols = np.array(rcols)
    wc = np.isin(rcols,list(fcdn.dtype.names))
    rcols = rcols[wc]
    print('columns sampled from data are:')
    print(rcols)

    common.addNS(ffc)
    #ffc['PHOTSYS'] = 'N'
    #sel = ffc['DEC'] < 32.375 #this is imperfect for the SGC (some of it is > 32.375 but still DECaLS), fix in the future
    #ffc['PHOTSYS'][sel] = 'S'
    #wn = ffc['PHOTSYS'] == 'N'

    #ffc.keep_columns(kc)
    #outf =  fl+wzm+str(rann)+'_clustering.ran.fits'
    #comments = ["DA02 'clustering' LSS catalog for random number "+str(rann)+", all regions","entries are only for data with good redshifts"]
    #common.write_LSS(ffc,outf,comments)

    outfn =  fl+wzm+'N_'+str(rann)+'_clustering.ran.fits'
    
    ffcn = ffc[wn]
    inds = np.random.choice(len(fcdn),len(ffcn))
    dshuf = fcdn[inds]
    for col in rcols:
        ffcn[col] = dshuf[col]
        kc.append(col)
    ffcn.keep_columns(kc)
    
    comments = ["DA02 'clustering' LSS catalog for random number "+str(rann)+", BASS/MzLS region","entries are only for data with good redshifts"]
    common.write_LSS(ffcn,outfn,comments)

    outfs =  fl+wzm+'S_'+str(rann)+'_clustering.ran.fits'
    fcds = Table.read(fl+wzm+'S_clustering.dat.fits')
    ffcs = ffc[~wn]
    inds = np.random.choice(len(fcds),len(ffcs))
    dshuf = fcds[inds]
    for col in rcols:
        ffcs[col] = dshuf[col]
    ffcs.keep_columns(kc)
    #edit these comments at some point
    comments = ["DA02 'clustering' LSS catalog for random number "+str(rann)+", DECaLS region","entries are only for data with good redshifts"]
    common.write_LSS(ffcs,outfs,comments)


    
def mock_equal_data_density(mockdir, datadir, outdir, tracer, region, zmin, zmax, nran, ran_dens):
    
    # Read mock clusering catalog
    # Read data clustering catalog
    # Read random catalog(s) 
    
    tracer_dat = tracer
    if tracer == "ELG":
        tracer_dat = "ELG_LOPnotqso"

    mock_fn = mockdir + "/" + tracer + "_" + region + "_clustering.dat.fits"
    data_fn = datadir + "/" + tracer_dat + "_" + region + "_clustering.dat.fits" 
    randoms_fn = {}
    for i in range(nran):
        randoms_fn[i] = datadir + "/" + tracer_dat + "_" + region + "_" + str(i) +"_clustering.ran.fits"

    mock_randoms_fn = mockdir + "/" + tracer + "_1_full.ran.fits"

    # Load mocks and select z range
    mock_tab = Table.read(mock_fn)
    sel_z = mock_tab["Z"] > zmin
    sel_z &= mock_tab["Z"] < zmax
    mock_tab = mock_tab[sel_z]
    num_mock_z_range = len(mock_tab)

    # Get length of data clustering catalog
    # Get length of randoms (divide by nran)
    num_data = fitsio.read_header(data_fn, ext =1)['NAXIS2']
    num_ran = 0
    for i in randoms_fn:
        num_ran += fitsio.read_header(randoms_fn[i], ext =1)['NAXIS2']    
    num_ran = num_ran / nran

    mock_ran_tab = Table.read(mock_randoms_fn)
    sel_z_ran = mock_ran_tab["PHOTSYS"] == region
    mock_ran_tab = mock_ran_tab[sel_z_ran]
    num_mock_ran = len(mock_ran_tab)


    fraction = (2500 * num_data * num_mock_ran / num_ran) / (ran_dens * num_mock_z_range)

    print("Mock fraction to be used to equal data n(z):", fraction)

    fraction_len = round(fraction * num_mock_z_range)

    rep_state = False
    sel_fraction = np.random.choice(num_mock_z_range, size = fraction_len, replace = rep_state) 

    out_fn = outdir + "/" + tracer + "_" + region + "_clustering.dat.fits"  

    print("Writing region", region, "to", out_fn)

    mock_tab[sel_fraction].write(out_fn, overwrite = True)
    print("Mock fraction catalogue writing complete.")

def create_collision_from_pota(fin, fout):
    print('Creating collision file from pota file', fin)
    df = fitsio.read(fin.replace('global','dvs_ro'))
    selcoll = df['COLLISION'] == True
    df = df[selcoll]
    print('size of collisions', len(df))
    common.write_LSS(df, fout, extname='COLLISION')
    return fout



def createrancomb_wdupspec(outdir, ranfile, alltileloc, mockassign, fdataspec):
    print('reading PRIORITY from mock and save it to random zdone')
    mockspec = Table(fitsio.read(mockassign.replace('global','dvs_ro'),columns=['LOCATION','TILEID','PRIORITY']))
    dataspec = Table(fitsio.read(fdataspec.replace('global','dvs_ro'), columns=['LOCATION','TILEID','TSNR2_ELG','TSNR2_LYA','TSNR2_BGS','TSNR2_QSO','TSNR2_LRG']))

    specD = join(mockspec, dataspec, keys=['LOCATION','TILEID'], join_type = 'left')

    #mockspec['TILELOCID'] = 10000*mockspec['TILEID'] +mockspec['LOCATION']
    #mockspec.keep_columns(['TILELOCID','PRIORITY'])
    randomdata = Table(fitsio.read(ranfile.replace('global','dvs_ro'), columns=['LOCATION','FIBER','TARGETID','RA','DEC','TILEID']))
    #randomdata['TILELOCID'] = 10000*randomdata['TILEID'] +randomdata['LOCATION']
    randomdata = join(randomdata, specD, keys=['LOCATION','TILEID'], join_type='left')

    randomdata.write(os.path.join(outdir, ranfile.split('/')[-1]), overwrite=True)
    if alltileloc is not None:
        print('copying alltileloc from spec dir to mock dir')
        shutil.copy(alltileloc, os.path.join(outdir, alltileloc.split('/')[-1]))
        return os.path.join(outdir, ranfile.split('/')[-1]), os.path.join(outdir, alltileloc.split('/')[-1])
    else:
        return os.path.join(outdir, ranfile.split('/')[-1]), None

def calc_weight_nt_misspw(data_full):
    if len(np.shape(data_full['BITWEIGHTS'])) == 1:
        nbits = 64
    else:        
        nbits = 64 * np.shape(data_full['BITWEIGHTS'])[1]
    recurr_full = data_full['PROB_OBS']*nbits
    wiip_full = (nbits+1)/(recurr_full+1)
    zerop_msk_full = (data_full['PROB_OBS']==0) & (~data_full['LOCATION_ASSIGNED'])
    ntile_range = [0, 10]
    idx_nt_full = []
    idx_nt_full_zp = []
    idx_nt_full_la = []
    for n in range(ntile_range[0], ntile_range[1]+1):
        idx_nt_full.append(np.where(data_full['NTILE'] == n))
        idx_nt_full_zp.append(np.where((data_full['NTILE'] == n) & (zerop_msk_full)))
        idx_nt_full_la.append(np.where((data_full['NTILE'] == n) & (data_full['LOCATION_ASSIGNED'])))
    s1 = 0
    s2 = 0
    s3 = 0
    #f_ntzp = np.ones(ntile_range[1]+1)
    f_ntmisspw = np.ones(ntile_range[1]+1)
    for n in range(ntile_range[0], ntile_range[1]+1):
        n_nt = len(idx_nt_full[n][0])
        n_ntzp = len(idx_nt_full_zp[n][0])
        n_ntwiip = np.sum(wiip_full[idx_nt_full_la[n][0]])
        n_ntmisspw = n_nt - n_ntwiip
#        print(n, n_nt, n_ntzp, n_ntmisspw)
        s1 += n_nt
        s2 += n_ntzp
        s3 += n_ntwiip
        #if n_nt > 0: f_ntzp[n] = n_ntzp / n_nt
        if n_nt > 0: f_ntmisspw[n] = n_ntmisspw / n_nt

    #w_ntzp_full = np.ones(len(data_full['RA']))
    w_ntmisspw_full = np.ones(len(data_full['RA']))
    for n in range(ntile_range[0], ntile_range[1]+1):
        #if 1-f_ntzp[n] > 0: w_ntzp_full[idx_nt_full[n][0]] = 1/(1-f_ntzp[n])
        if 1-f_ntmisspw[n] > 0: w_ntmisspw_full[idx_nt_full[n][0]] = 1/(1-f_ntmisspw[n])

    data_full['WEIGHT_NT_MISSPW'] = w_ntmisspw_full
    return data_full, f_ntmisspw

def calc_weight_nt_misspw_ran(data_full_ran, f_ntmisspw):
    ntile_range = [0, 10]
    idx_nt_full_ran = []
    for n in range(ntile_range[0], ntile_range[1]+1):
        idx_nt_full_ran.append(np.where(data_full_ran['NTILE'] == n))

#    f_ntmisspw = np.ones(ntile_range[1]+1)
#    w_ntzp_full_ran = np.ones(len(data_full_ran['RA']))
    w_ntmisspw_full_ran = np.ones(len(data_full_ran['RA']))
    for n in range(ntile_range[0], ntile_range[1]+1):
    #    w_ntzp_full_ran[idx_nt_full_ran[n][0]] = 1-f_ntzp[n]
        w_ntmisspw_full_ran[idx_nt_full_ran[n][0]] = 1-f_ntmisspw[n]
    data_full_ran['WEIGHT_NT_MISSPW'] = w_ntmisspw_full_ran
    return data_full_ran

def do_weight_nt_misspw(fb, ranmin=0, ranmax=18, par='n', dirout=None):

    if dirout is not None:
        fb_full_destiny = os.path.join(dirout, fb.split('/')[-1])
    else:
        fb_full_destiny = fb

    fb_full = fb + '_full_HPmapcut.dat.fits'
    data_to_concat = []
    rans_sgc_to_concat = []
    rans_ngc_to_concat = []
    
    print('aqui',fb_full)
    inputdata = Table(pf.open(fb_full)[1].data)

    ngc_mask = common.splitGC(inputdata)

    for cap in ['NGC', 'SGC']:
        #selngc = common.splitGC(inputdata)
        if cap == 'SGC':
            datacap = inputdata[~ngc_mask]
        else:
            datacap = inputdata[ngc_mask]

        datacap, f_ntmisspw = calc_weight_nt_misspw(datacap)

        data_to_concat.append(datacap)

        for rn in range(ranmin, ranmax):
            fb_full = Table.read(fb + '_%d_full_HPmapcut.ran.fits' % rn)
            selngc = common.splitGC(fb_full)
            #ngc_mask = common.splitGC(fb_full)   

            if cap == 'SGC':
                datacap = calc_weight_nt_misspw_ran(fb_full[~selngc], f_ntmisspw)
                rans_sgc_to_concat.append(datacap)
            else:
            #    rans_sgc_to_concat.append(datacap)
            #else:

                datacap = calc_weight_nt_misspw_ran(fb_full[selngc], f_ntmisspw)
                rans_ngc_to_concat.append(datacap)

    concatenated_table = vstack(data_to_concat)
    common.write_LSS(concatenated_table, fb_full_destiny + '_full_HPmapcut.dat.fits')
    #concatenated_table.write(fb_full_destiny + '_full_HPmapcut.dat.fits')
    for ran_sgc, ran_ngc, rn in zip(rans_sgc_to_concat, rans_ngc_to_concat, np.arange(ranmin, ranmax)):
        concatenated_table = vstack([ran_sgc, ran_ngc])
        common.write_LSS(concatenated_table, fb_full_destiny + '_{RANNUM}_full_HPmapcut.ran.fits'.format(RANNUM=rn))
        #concatenated_table.write(fb_full_destiny + '_{RANNUM}_full_HPmapcut.ran.fits'.format(RANNUM=rn))
    return True


