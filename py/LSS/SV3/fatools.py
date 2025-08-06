##from desiutil.iers import freeze_iers
##freeze_iers()

import fitsio
import numpy as np
from astropy.table import Table,join
# system
import os
import subprocess
import sys
import tempfile
import shutil
import re
import glob

# time
from time import time
from datetime import datetime, timedelta

# desi
import desitarget
from desitarget import io 
from desitarget.mtl import inflate_ledger
from desiutil.log import get_logger

log = get_logger()


#hardcode target directories; these are fixed
#JL - These change once at the very beginning of observations.
#JL - Adding a format string to change versioning of photometry. 

skydir = '/global/cfs/cdirs/desi/target/catalogs/dr9/0.57.0/skies'
skydirMain = '/global/cfs/cdirs/desi/target/catalogs/dr9/1.1.1/skies'
#skydirMain = '/global/cfs/cdirs/desi/target/catalogs/dr9/{0}/skies'
tdir = '/global/cfs/cdirs/desi/target/catalogs/dr9/0.57.0/targets/sv3/resolve/'
tdirMain = '/global/cfs/cdirs/desi/target/catalogs/dr9/1.1.1/targets/main/resolve/'
#tdirMain = '/global/cfs/cdirs/desi/target/catalogs/dr9/{0}/targets/main/resolve/'
# AR default REF_EPOCH for PMRA=PMDEC=REF_EPOCH=0 objects
gaia_ref_epochs = {"dr2": 2015.5}


minimal_target_columns= ['RELEASE','BRICKNAME','BRICKID','BRICK_OBJID','MORPHTYPE','RA',\
'DEC','EBV','FLUX_G','FLUX_R','FLUX_Z','FLUX_W1','FLUX_W2','FLUX_IVAR_G','FLUX_IVAR_R',\
'FLUX_IVAR_Z','FLUX_IVAR_W1','FLUX_IVAR_W2','FIBERFLUX_G','FIBERFLUX_R','FIBERFLUX_Z',\
'FIBERTOTFLUX_G','FIBERTOTFLUX_R','FIBERTOTFLUX_Z','REF_EPOCH','MASKBITS','SERSIC',\
'SHAPE_R','SHAPE_E1','SHAPE_E2','REF_ID','REF_CAT','GAIA_PHOT_G_MEAN_MAG',\
'GAIA_PHOT_BP_MEAN_MAG','GAIA_PHOT_RP_MEAN_MAG','PARALLAX','PMRA','PMDEC','PHOTSYS',\
'TARGETID','SUBPRIORITY','OBSCONDITIONS','PRIORITY_INIT','NUMOBS_INIT','SV3_DESI_TARGET',\
'SV3_BGS_TARGET','SV3_MWS_TARGET','SV3_SCND_TARGET']

minimal_target_columns_main= ['RELEASE','BRICKNAME','BRICKID','BRICK_OBJID','MORPHTYPE','RA',\
'DEC','EBV','FLUX_G','FLUX_R','FLUX_Z','FLUX_W1','FLUX_W2','FLUX_IVAR_G','FLUX_IVAR_R',\
'FLUX_IVAR_Z','FLUX_IVAR_W1','FLUX_IVAR_W2','FIBERFLUX_G','FIBERFLUX_R','FIBERFLUX_Z',\
'FIBERTOTFLUX_G','FIBERTOTFLUX_R','FIBERTOTFLUX_Z','REF_EPOCH','MASKBITS','SERSIC',\
'SHAPE_R','SHAPE_E1','SHAPE_E2','REF_ID','REF_CAT','GAIA_PHOT_G_MEAN_MAG',\
'GAIA_PHOT_BP_MEAN_MAG','GAIA_PHOT_RP_MEAN_MAG','PARALLAX','PMRA','PMDEC','PHOTSYS',\
'TARGETID','SUBPRIORITY','OBSCONDITIONS','PRIORITY_INIT','NUMOBS_INIT','DESI_TARGET',\
'BGS_TARGET','MWS_TARGET','SCND_TARGET']

def comp_neworig(tileid,dirn='/global/cfs/cdirs/desi/survey/catalogs/testfiberassign/SV3rerun/orig/', verbose = False):
    """
    check that new matches the original
    
    """
    ts = str(tileid).zfill(6)
    fa = fitsio.read('/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'+ts[:3]+'/fiberassign-'+ts+'.fits.gz') 
    fn = fitsio.read(dirn+'fba-'+ts+'.fits')
    w = fn['DEVICE_TYPE'] == 'POS'
    fn = fn[w]
    wn = fn['TARGETID'] >= 0
    fn = fn[wn]
    print(len(fn))
    wa = fa['TARGETID'] >= 0
    fa = fa[wa]
    print(len(fa))  
    ws = np.isin(fn['TARGETID'],fa['TARGETID'])
    print(np.sum(ws))   
    if np.sum(ws) == len(fa) and len(fa) == len(fn):
        return True
    else:
        return False

def comp_neworig_tgt(tileid, verbose = False):
    """
    check that new matches the original, just tgt
    
    """
    ts = str(tileid).zfill(6)
    fa = fitsio.read('/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'+ts[:3]+'/fiberassign-'+ts+'.fits.gz')
    dirn =  '/global/cfs/cdirs/desi/survey/catalogs/testfiberassign/SV3rerun/orig/'
    fn = fitsio.read(dirn+'fba-'+ts+'.fits')
    wn = fn['FA_TYPE'] != 0
    wn &= fn['FA_TYPE'] != 4

    #w = fn['DEVICE_TYPE'] == 'POS'
    fn = fn[wn]
    #wn = fn['TARGETID'] >= 0
    #fn = fn[wn]
    print(len(fn))
    wa = fa['FA_TYPE'] != 0
    wa &= fa['FA_TYPE'] != 4
    fa = fa[wa]
    print(len(fa))  
    ws = np.isin(fn['TARGETID'],fa['TARGETID'])
    print(np.sum(ws))   
    if np.sum(ws) == len(fa) and len(fa) == len(fn):
        return True
    else:
        return False

def comp_neworig_fba(tileid,dirn =  '/global/cfs/cdirs/desi/survey/catalogs/testfiberassign/SV3rerun/orig/', verbose = False):
    """
    check that new matches the original, comparing fba files
    
    """
    ts = str(tileid).zfill(6)
    ts = str(tileid).zfill(6)
    #get info from origin fiberassign file
    fht = fitsio.read_header('/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'+ts[:3]+'/fiberassign-'+ts+'.fits.gz')
    indir = fht['OUTDIR']
    if fht['DESIROOT'] == '/data/datasystems':
        indir = '/global/cfs/cdirs/desi/survey/fiberassign/SV3/' +fht['PMTIME'][:10].translate({ord('-'): None})  +'/'      
        try:
            
            f = fitsio.read(indir+ts+'-targ.fits')
        except:
            
            date = int(fht['PMTIME'][:10].translate({ord('-'): None}))-1
            indir = '/global/cfs/cdirs/desi/survey/fiberassign/SV3/'+str(date)+'/'

           
    fa = fitsio.read(indir+'fba-'+ts+'.fits')
    
    fn = fitsio.read(dirn+'fba-'+ts+'.fits')
    return np.array_equal(fa['TARGETID'],fn['TARGETID'])
#     w = fn['DEVICE_TYPE'] == 'POS'
#     fn = fn[w]
#     wn = fn['TARGETID'] >= 0
#     fn = fn[wn]
#     #print(len(fn))
#     wa = fa['OBJTYPE'] == 'TGT'
#     fa = fa[wa]
#     #print(len(fa))  
#     ws = np.isin(fn['TARGETID'],fa['TARGETID'])
#     #print(np.sum(ws))   
#     if np.sum(ws) == len(fa):# and len(fa) == len(fn):
#         return True
#     else:
#         return False

 
def redo_fba_fromorig(tileid,outdir=None,faver=None, verbose = False,survey='main',faenvdir=None):
    '''
    simply try to reproduce fiberassign from the files in the fiberassign directory
    '''
    ts = str(tileid).zfill(6)
    #get info from origin fiberassign file
    fht = fitsio.read_header('/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'+ts[:3]+'/fiberassign-'+ts+'.fits.gz')
    if survey == 'SV3':
        indir = fht['OUTDIR']
        if fht['DESIROOT'] == '/data/datasystems':
            indir = '/global/cfs/cdirs/desi/survey/fiberassign/SV3/' +fht['PMTIME'][:10].translate({ord('-'): None})  +'/'      
            try:
                f = fitsio.read(indir+ts+'-targ.fits')
            except:
        
                date = int(fht['PMTIME'][:10].translate({ord('-'): None}))-1
                indir = '/global/cfs/cdirs/desi/survey/fiberassign/SV3/'+str(date)+'/'
    else:
        indir = '/global/cfs/cdirs/desi/survey/fiberassign/'+survey+'/'+ts[:3]+'/'
    tarf = indir+ts+'-targ.fits'
    try:
        fitsio.read(tarf)
    except:
        return('Error! target file does not appear to exist for tile '+ts+' '+tilef)    
    tilef = indir+ts+'-tiles.fits'
    try:
        fitsio.read(tilef)
    except:
        return('Error! tile file does not appear to exist for tile '+ts+' '+tilef)    
    skyf = indir+ts+'-sky.fits'
    try:
        fitsio.read(skyf)
    except:
        log.critical('Error! sky file does not appear to exist')    
    scndf = indir+ts+'-scnd.fits'
    scnd = True 
    try:
        fitsio.read(scndf)
    except:
        log.info(' secondary file does not appear to exist')
        scnd = False 
           
    gfaf = indir+ts+'-gfa.fits'
    try:
        fitsio.read(gfaf)
    except:
        log.info('Error! gfa file does not appear to exist')    
    toof = indir+ts+'-too.fits'
    too = os.path.isfile(toof)
    if too:
        print('will be using too file '+toof)
    if outdir is None:
        outdir = '/global/cfs/cdirs/desi/survey/catalogs/testfiberassign/'+survey+'rerun/orig/'
      
    prog = fht['FAPRGRM'].lower()
    gaiadr = None
    if np.isin('gaiadr2',fht['FAARGS'].split()):
        gaiadr = 'dr2'
    if np.isin('gaiaedr3',fht['FAARGS'].split()):
        gaiadr = 'edr3'
    
    fo = open(outdir+'fa-'+ts+'.sh','w')
    fo.write('#!/bin/bash\n\n')
    fo.write('source /global/common/software/desi/desi_environment.sh main\n')
    if faver is None:
        faver = float(fht['FA_VER'][:3])
        if 'main' in indir:
            assert(faver > 3.0)
        if faver == 2.4:
            fo.write('export SKYBRICKS_DIR=${DESI_ROOT}/target/skybricks/v2\n')

        if faver < 2.4:
            if int(indir[-7:-1]) > 210413:
                fo.write("module swap fiberassign/2.3.0\n")
            else:
                fo.write("module swap fiberassign/"+fht['FA_VER'][:3]+'.0'+"\n")
        elif faver >= 5.0:
            fo.write("module swap fiberassign/main")
        else:
            assert(faver < 5.0)
            fo.write("module swap fiberassign/"+fht['FA_VER']+"\n")
    elif faver == 'current':
        faver = 5
        if faenvdir is not None:
            fo.write('PATH='+faenvdir+'fiberassign/bin:$PATH\n')
            fo.write('PYTHONPATH='+faenvdir+'/fiberassign/py:$PYTHONPATH\n')
            
    else:
        faver = float(faver[:3])
        if 'main' in indir:
            assert(faver > 3.0)
        if faver >= 5.0:
            fo.write("module swap fiberassign/main")
        else:
            assert(faver < 5.0)
            fo.write("module swap fiberassign/"+str(faver)+"\n")
       
    fo.write("fba_run")
    fo.write(" --targets "+tarf)
    if scnd:
        fo.write(" "+scndf)
    if too:
        fo.write(" "+toof)
    fo.write(" --sky "+skyf)
    fo.write(" --footprint "+tilef)
    rundate= fht['RUNDATE']
    if rundate == '2021-04-10T21:28:37':
        rundate = '2021-04-10T20:00:00'
    fo.write(" --rundate "+rundate)
    fo.write(" --fieldrot "+np.format_float_positional(fht['FIELDROT']))
    fo.write(" --dir "+outdir)
    #if indir != '/global/cfs/cdirs/desi/survey/fiberassign/SV3/20210416/' and indir != '/global/cfs/cdirs/desi/survey/fiberassign/SV3/20210418/':
    fo.write(" --sky_per_petal 40 --standards_per_petal 10")
    #fo.write(" --by_tile true")
    if faver >= 2.4:
        fo.write(" --sky_per_slitblock 1")
    if faver >= 3:
        fo.write(" --ha "+str(fht['FA_HA']))
        fo.write(" --margin-gfa 0.4 --margin-petal 0.4 --margin-pos 0.05")
    fo.close()    
 
        
def get_fba_fromnewmtl(tileid,mtldir=None,getosubp=False,outdir=None,faver=None, overwriteFA = False,newdir=None, verbose = False, mock = False, targver = '1.1.1', reproducing = False):
    ts = str(tileid).zfill(6)
    #get info from origin fiberassign file
    fa_fn = '/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'+ts[:3]+'/fiberassign-'+ts+'.fits.gz'
    fht = fitsio.read_header(fa_fn)
    indir = fht['OUTDIR']
    if (fht['DESIROOT'] == '/data/datasystems') and not ( ('holding' in indir.lower()) or ('main' in indir.lower())):
        indir = '/global/cfs/cdirs/desi/survey/fiberassign/SV3/' +fht['PMTIME'][:10].translate({ord('-'): None})  +'/'      
        try:
            f = fitsio.read(indir+ts+'-targ.fits')
        except:
            date = int(fht['PMTIME'][:10].translate({ord('-'): None}))-1
            indir = '/global/cfs/cdirs/desi/survey/fiberassign/SV3/'+str(date)+'/'
            
    elif ( ('holding' in indir.lower()) or ('main' in indir.lower())):
        indir = '/global/cfs/cdirs/desi/survey/fiberassign/main/' + ts[0:3] +'/'
    log.info('In directory: {}'.format(indir))

    tilef = indir+ts+'-tiles.fits'
    try:
        fitsio.read(tilef)
    except:
        '''
        try:
            if 'sv3' in indir.lower():
                date = int(fht['PMTIME'][:10].translate({ord('-'): None}))-1
                indir = '/global/cfs/cdirs/desi/survey/fiberassign/SV3/'+str(date)+'/'
            elif ('main' in indir.lower()) or ('holding' in indir.lower()):
                ###WHY NOT DATE SET HERE???
                log.info('WHY NOT DATE SET HERE FOR MAIN')
                indir = '/global/cfs/cdirs/desi/survey/fiberassign/main/' + ts[0:3] +'/'
            else:
                raise ValueError('survey not sv3 or main, will have checks for SV2/1/CMX in future.')
            tilef = indir+ts+'-tiles.fits'
            fitsio.read(tilef)

        except:
        '''
        log.critical('failed to read tile file')
        log.critical('Error! tile file does not appear to exist for tile '+ts+' '+tilef)
        log.critical('indir')
        log.critical(indir)
        return('Error! tile file does not appear to exist for tile '+ts+' '+tilef)
    skyf = indir+ts+'-sky.fits'
    try:
        fitsio.read(skyf)
    except:
        log.critical('Error! sky file does not appear to exist')    
    scndf = indir+ts+'-scnd.fits'
    scnd = True 
    try:
        fitsio.read(scndf)
    except:
        log.info(' secondary file does not appear to exist')
        scnd = False 
    gfaf = indir+ts+'-gfa.fits'
    try:
        fitsio.read(gfaf)
    except:
        log.critical('Error! gfa file does not appear to exist')   
    toof = indir+ts+'-too.fits'
    too = os.path.isfile(toof)
    if too:
        log.info('will be using too file '+toof)
    if outdir is None:
        outdir = '/global/cfs/cdirs/desi/survey/catalogs/testfiberassign/SV3rerun/'
    if getosubp == True or (mtldir == None and newdir == None):
        outdir += 'orig/'
    if newdir == None:
        if mtldir == None:
            tarfn = indir+ts+'-targ.fits' 
        else:
            tarfn = outdir+ts+'-targ.fits'   
    else:
        tarfn = newdir+ts+'-targ.fits' 
    prog = fht['FAPRGRM'].lower()
    gaiadr = None
    if np.isin('gaiadr2',fht['FAARGS'].split()):
        gaiadr = 'dr2'
    if np.isin('gaiaedr3',fht['FAARGS'].split()):
        gaiadr = 'edr3'

    if mtldir is not None:
        if 'sv3' in indir.lower():
            if verbose:
                log.debug('sv3 survey')
            altcreate_mtl(tilef,
            mtldir+prog,        
            gaiadr,
            fht['PMCORR'],
            tarfn,
            tdir+prog,
            mock = mock)
        elif ('main' in indir.lower()) or ('holding' in indir.lower()):
            if verbose:
                log.info('main survey')
            if (not reproducing) or (targver == '1.1.1'):
                if verbose:
                    log.info('if reproducing is True, targver must be 1.1.1')
                    log.info(f'targver  = {targver}')
                    log.info(f'reproducing = {reproducing}')
                altcreate_mtl(tilef,
                mtldir+prog,        
                gaiadr,
                fht['PMCORR'],
                tarfn,
                tdirMain.format(targver)+prog,
                survey = 'main',
                mock = mock)
            #tdirMain+prog,
            #LGN 06/11/25: Adding compatibility with reproducing tests for DARK1B
            elif targver == '1.0.0' or targver == '3.0.0':
                if verbose:
                    log.info('targver must be 1.0.0 (or at least not 1.1.1) and reproducing must be True')
                    log.info(f'targver  = {targver}')
                    log.info(f'reproducing = {reproducing}')
                    
                shutil.copyfile(indir+ts+'-targ.fits', tarfn)

        else:
            log.critical('invalid input directory. must contain either sv3, main, or holding')
            raise ValueError('indir must contain either sv3, main, or holding')


    #LGN 06/11/25: Moving this outside of loop to work with DARK1B program
    if not os.path.exists(outdir):
        log.info('running makedirs. making {0}'.format(outdir))
        os.makedirs(outdir)
    
    if getosubp:
        if tileid == 315:
            log.info('special tile 315 case triggered')
            log.info('tileid = {0}'.format(tileid))
            
            otar = Table.read(indir+ts+'-targ.fits')
            otar.keep_columns(['TARGETID','PRIORITY','SUBPRIORITY'])
            
            ntar = Table.read(tarfn)
            ntar.remove_columns(['SUBPRIORITY', 'PRIORITY'])
            ntar = join(ntar,otar,keys=['TARGETID'])
            ntar.write(tarfn,format='fits', overwrite=True)
        else:
            
            otar = Table.read(indir+ts+'-targ.fits')
            otar.keep_columns(['TARGETID','SUBPRIORITY'])

            ntar = Table.read(tarfn)
            ntar.remove_columns(['SUBPRIORITY'])
            ntar = join(ntar,otar,keys=['TARGETID'])
            ntar.write(tarfn,format='fits', overwrite=True)
    
    fo = open(outdir+'fa-'+ts+'.sh','w')
    fo.write('#!/bin/bash\n\n')
    fo.write('source /global/common/software/desi/desi_environment.sh main\n')

    if faver == None:
        faver = float(fht['FA_VER'][:3])
        if faver == 2.4:
            fo.write('export SKYBRICKS_DIR=${DESI_ROOT}/target/skybricks/v2\n')

        if faver < 2.4:
            if int(indir[-7:-1]) > 210413:
                fo.write("module swap fiberassign/2.3.0\n") #inspection of results revealed tiles that used 2.2.dev* after 20210413 are reproduced using 2.3.0 and those before using 2.2.0
            else:
                fo.write("module swap fiberassign/"+fht['FA_VER'][:3]+'.0'+"\n")
        elif faver >= 2.4 and faver < 4.0:
            fo.write("module swap fiberassign/"+fht['FA_VER']+"\n")
        elif faver >= 4.0 and faver < 5.0:
            #fo.write("module swap fiberassign/5.0.0\n")
            fo.write("export PATH=/global/cfs/cdirs/desi/users/raichoor/fiberassign-rerun-main/fiberassign_main_godesi23.10/bin:$PATH\n")
            fo.write("export PYTHONPATH=/global/cfs/cdirs/desi/users/raichoor/fiberassign-rerun-main/fiberassign_main_godesi23.10/py:$PYTHONPATH\n")
            fo.write("export SKYHEALPIXS_DIR=$DESI_ROOT/target/skyhealpixs/v1\n")
        else:
            #with fiberassign refactor as of May 2025, it is no longer necessary to specify a module swap for faver > 5.00
            assert faver >= 5.0
    else:
        fo.write("module swap fiberassign/"+str(faver)+"\n")
        faver = float(faver[:3])
    fo.write("fba_run")
    fo.write(" --targets "+tarfn)
    if scnd:
        fo.write(" "+scndf)
    if too:
        fo.write(" "+toof)
    fo.write(" --sky "+skyf)
    fo.write(" --footprint "+tilef)
    rundate= fht['RUNDATE']
    if rundate == '2021-04-10T21:28:37':
        rundate = '2021-04-10T20:00:00'
    fo.write(" --rundate "+rundate)
    fo.write(" --fieldrot "+np.format_float_positional(fht['FIELDROT']))
    fo.write(" --dir "+outdir)
    fo.write(" --sky_per_petal 40 --standards_per_petal 10")
    if overwriteFA:
        fo.write(" --overwrite")
    #fo.write(" --by_tile true")
    if faver >= 2.4:
        fo.write(" --sky_per_slitblock 1")
    if faver >= 3:
        fo.write(" --ha "+str(fht['FA_HA']))
        fo.write(" --margin-gfa 0.4 --margin-petal 0.4 --margin-pos 0.05")
    if faver >=5:
        fo.write(" --fafns_for_stucksky "+fa_fn)
    fo.close()    

#     if float(fht['FA_VER'][:3]) < 2.4:
#         fo.write("module swap fiberassign/2.3.0\n")
#     else:
#         fo.write("module swap fiberassign/"+fht['FA_VER']+"\n")
#     fo.write("fba_run")
#     fo.write(" --targets "+tarfn+" "+scndf)
#     if too:
#         fo.write(" "+toof)
#     fo.write(" --sky "+skyf)
#     fo.write(" --footprint "+tilef)
#     fo.write(" --rundate "+fht['RUNDATE'])
#     fo.write(" --fieldrot "+str(fht['FIELDROT']))
#     fo.write(" --dir "+outdir)
#     #fo.write(" --by_tile true")
#     if float(fht['FA_VER'][:3]) >= 3:
#         fo.write(" --ha "+str(fht['FA_HA']))
#     fo.close()    


def altcreate_mtl(
    tilesfn,
    mtldir,            
    gaiadr,
    pmcorr,
    outfn,
    targdir,
    survey='sv3',
    mtltime=None,#I think we will just want this to be the latest for the re/alt runs    
    pmtime_utc_str=None,
    add_plate_cols=True,
    verbose=False,
    mock=False
    #tmpoutdir=tempfile.mkdtemp(),
):
    """
    Mostly copied from fiberassign.fba_launch_io.create_mtl
    Create a (primary or secondary) target fits file, based on MTL ledgers (and complementary columns from desitarget targets files).
    
    Args:
        tilesfn: path to a tiles fits file (string)
        mtldir: folder with ledger files        
        targdir: desitarget targets folder (or file name if secondary) (string)        
        gaiadr: Gaia dr ("dr2" or "edr3")
        pmcorr: apply proper-motion correction? ("y" or "n")
        outfn: fits file name to be written (string)
        survey: sv3 or main
        mtltime: MTL isodate (string formatted as yyyy-mm-ddThh:mm:ss+00:00); this needs be considered carefully for alt mtls
        tmpoutdir (optional, defaults to a temporary directory): temporary directory where
                write_targets will write (creating some sub-directories)
        pmtime_utc_str (optional, defaults to None): UTC time use to compute
                new coordinates after applying proper motion since REF_EPOCH
                (string formatted as "yyyy-mm-ddThh:mm:ss+00:00")
        add_plate_cols (optional, defaults to True): adds a PLATE_RA and PLATE_DEC columns (boolean)
        log (optional, defaults to Logger.get()): Logger object
        step (optional, defaults to ""): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start(optional, defaults to time()): start time for log (in seconds; output of time.time()
    Notes:
        if pmcorr="y", then pmtime_utc_str needs to be set; will trigger an error otherwise.
        for sv3-backup, we remove BACKUP_BRIGHT targets.
        TBD : if secondary targets, we currently disable the inflate_ledger(), as it
                seems to not currently work.
                hence if secondary and pmcorr="y", the code will crash, as the 
                GAIA_ASTROMETRIC_EXCESS_NOISE column will be missing; though we do not
                expect this configuration to happen, so it should be fine for now.
        TBD: the PLATE_{RA,DEC,REF_EPOCH} columns currently simply are copy of RA,DEC,REF_EPOCH
        TBD:    but it prepares e.g. to add chromatic offsets.
        20210526 : implementation of using subpriority=False in write_targets
                    to avoid an over-writting of the SUBPRIORITY; AJR changed to True reproduce SV3
    """
    tiles = fitsio.read(tilesfn)
    tileIDs = tiles['TILEID']
    # AR mtl: read mtl
    if (315 in tileIDs) and (len(tiles) == 1):
        log.info('special handling of tile 315 for SV3')
        assert(survey.lower() == 'sv3')
        d0 = io.read_targets_in_tiles(
            mtldir,
            tiles,
            quick=False,
            mtl=True,
            unique=False,
            isodate=None,
        )
        mtltime = np.unique(d0[d0['ZTILEID'] == 314]['TIMESTAMP'])

        assert(mtltime.shape[0] == 1)
        
        mtltime = str(mtltime[0])
        
        d = io.read_targets_in_tiles(
            mtldir,
            tiles,
            quick=False,
            mtl=True,
            unique=True,
            isodate=mtltime,
            verbose=verbose,
        )
    elif (315 in tileIDs) and (len(tiles) > 1):
        log.critical('315 in tiles but multiple tiles provided')
        log.critical(tileIDs)
        log.critical(tiles)
        raise ValueError('When processing tile 315, code should strip out processing of all other tiles. ')
    else:
        # LGN: Adding handling for bright1b/dark1b tiles
        if 'dark1b' or 'bright1b' in mtldir:
            log.info('Running with maketwostyle=True')
            is_ext = True
            # LGN Formatting the path to bright or dark ledgers
            # LGN Passing list of directories to read_targets_in_tiles
            mtldir_short = os.path.join(
                os.path.dirname(mtldir), 
                os.path.basename(mtldir).replace('1b', '')
            )
            mtldirs = [mtldir_short, mtldir]
            
            d = io.read_targets_in_tiles(
                mtldirs,
                tiles,
                quick=False,
                mtl=True,
                unique=True,
                isodate=mtltime,
                verbose=verbose,
                tabform='ascii.ecsv',
                maketwostyle = is_ext
            )
        
        else:
            is_ext = False
        
            d = io.read_targets_in_tiles(
                mtldir,
                tiles,
                quick=False,
                mtl=True,
                unique=True,
                isodate=mtltime,
                verbose=verbose,
                tabform='ascii.ecsv',
                maketwostyle = is_ext
            )
        
    
    try:
        log.info('shape of read_targets_in_tiles output')
        log.info(d.shape)
        ntargs = d.shape  
    except:
        log.info('len of read_targets_in_tiles output post failure of shape')
        log.info(len(d))
        ntargs = len(d)
    assert(ntargs)
    # AR mtl: removing by hand BACKUP_BRIGHT for sv3/BACKUP
    # AR mtl: using an indirect way to find if program=backup,
    # AR mtl:   to avoid the need of an extra program argument
    # AR mtl:   for sv3, there is no secondary-backup, so no ambiguity
    if (survey == "sv3") & ("backup" in mtldir):
        from desitarget.sv3.sv3_targetmask import mws_mask

        keep = (d["SV3_MWS_TARGET"] & mws_mask["BACKUP_BRIGHT"]) == 0
        d = d[keep]

    #tcol = ['SV3_DESI_TARGET','SV3_BGS_TARGET','SV3_MWS_TARGET','SV3_SCND_TARGET']
    #for col in tcol:
    #    columns.append(col) 
    if not mock:
        log.info('len(d)= {0}'.format(len(d)))
        
        #LGN, we now only run inflate_ledger for SV3
        if survey == "sv3":
            columns = [key for key in minimal_target_columns if key not in d.dtype.names]
            d = inflate_ledger(d, targdir, columns=columns, header=False, strictcols=False, quick=True)    # AR adding PLATE_RA, PLATE_DEC, PLATE_REF_EPOCH ?
        
        log.info('len(d)= {0}'.format(len(d)))
        if add_plate_cols:
            d = Table(d)
            d["PLATE_RA"] = d["RA"]
            d["PLATE_DEC"] = d["DEC"]
            d["PLATE_REF_EPOCH"] = d["REF_EPOCH"]
            d = d.as_array()

        # AR mtl: PMRA, PMDEC: convert NaN to zeros
        d = force_finite_pm(d)

        # AR mtl: update RA, DEC, REF_EPOCH using proper motion?
        if pmcorr == "y":
            if pmtime_utc_str is None:
                sys.exit(1)
            d = update_nowradec(d, gaiadr, pmtime_utc_str)
        else:
            d = force_nonzero_refepoch(
                d, gaia_ref_epochs[gaiadr]
            )
    d = Table(d)

    outfndir = '/'.join(outfn.split('/')[:-1])
    if not os.path.exists(outfndir):
        os.makedirs(outfndir, exist_ok=True)
    try:
        log.info('shape of read_targets_in_tiles output pre writing')
        log.info(d.shape)  
    except:
        log.info('len of read_targets_in_tiles output pre writing post failure of shape')
        log.info(len(d))
    d.write(outfn,format='fits', overwrite=True)
    del d
    return True
    # AR mtl: write fits
    #n, tmpfn = io.write_targets(tmpoutdir, d, indir=mtldir, indir2=targdir, survey=survey, subpriority=True)
    #_ = mv_write_targets_out(tmpfn, tmpoutdir, outfn)

    # AR mtl: update header if pmcorr = "y"
    if pmcorr == "y":
        fd = fitsio.FITS(outfn, "rw")
        fd["TARGETS"].write_key("COMMENT", "RA,DEC updated with PM for AEN objects")
        fd["TARGETS"].write_key("COMMENT", "REF_EPOCH updated for all objects")
        fd.close()
    
"""
copying functions from fba_launch_io.py just so these are stable and in one place; don't
actually want to have to load proper version of fiberassign just for this
"""

def mv_write_targets_out(infn, targdir, outfn):
    """
    Moves the file created by desitarget.io.write_targets
    and removes folder created by desitarget.io.write_targets    
    
    Args:
        infn: filename output by desitarget.io.write_targets
        targdir: folder provided as desitarget.io.write_targets input
        outfn: desired renaming of infn
        log (optional, defaults to Logger.get()): Logger object
        step (optional, defaults to ""): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start(optional, defaults to time()): start time for log (in seconds; output of time.time()
    """
    # AR renaming
    _ = shutil.move(infn, outfn)
    # AR removing folders
    if targdir[-1] != "/":
        targdir = "{}/".format(targdir)
    tmpdirs = infn.replace(targdir, "").split("/")[:-1]
    for i in range(len(tmpdirs))[::-1]:
        os.rmdir(os.path.join(*[targdir] + tmpdirs[: i + 1]))

def get_nowradec(ra, dec, pmra, pmdec, parallax, ref_year, pmtime_utc_str, scnd=False):
    """
    Apply proper motion correction
    
    Args:
        ra: numpy array of RAs (deg)
        dec: numpy array of DECs (deg)
        pmra: numpy array of projected proper-motion in RA (mas/year)
        pmdec: numpy array of projected proper-motion in DEC (mas/year)
        parallax: numpy array of parallax (mas)
        ref_year: reference epoch (e.g. 2015.5 for Gaia/DR2)
        pmtime_utc_str: date to update position to (format: YYYY-MM-DDThh:mm:ss+00:00)
        scnd (optional, defaults to False): secondary target? (boolean; if True, sets parallax=0)
    Returns:
        ra: numpy array of RAs updated to pmtime_utc_str (deg)
        dec: numpy array of DECs updated to pmtime_utc_str (deg)
    Notes:
        Courtesy of DL; adapted from legacypipe.survey
        Originally named radec_at_mjd()
    """
    # AR pmtime_utc : UTC time of the new ref_epoch; "%Y-%m-%dT%H:%M:%S%z", e.g. "2021-04-21T00:00:00+00:00"
    # AR scnd=True -> parallax is set to 0, i.e. not used
    """
    Units:
    - matches Gaia DR1/DR2
    - pmra,pmdec are in mas/yr.
      pmra is in angular speed (ie, has a cos(dec) factor)
    - parallax is in mas.
    Returns: RA,Dec
    """
    equinox = 53084.28  # mjd of the spring equinox in 2004
    equinox_jyear = Time(equinox, format="mjd").jyear
    axistilt = 23.44  # degrees
    arcsecperrad = 3600.0 * 180.0 / np.pi
    # AR pmtime
    pmtime_utc = datetime.strptime(pmtime_utc_str, "%Y-%m-%dT%H:%M:%S%z")
    pmtime_utc_jyear = Time(pmtime_utc).jyear
    pmtime_utc_mjd = Time(pmtime_utc).mjd

    def xyztoradec(xyz):
        assert len(xyz.shape) == 2
        ra = np.arctan2(xyz[:, 1], xyz[:, 0])  # AR added "np." in front of arctan2...
        ra += 2 * np.pi * (ra < 0)
        norm = np.sqrt(np.sum(xyz ** 2, axis=1))
        dec = np.arcsin(xyz[:, 2] / norm)
        return np.rad2deg(ra), np.rad2deg(dec)

    def radectoxyz(ra_deg, dec_deg):  # AR changed inputs from ra,dec to ra_deg,dec_deg
        ra = np.deg2rad(ra_deg)
        dec = np.deg2rad(dec_deg)
        cosd = np.cos(dec)
        return np.vstack((cosd * np.cos(ra), cosd * np.sin(ra), np.sin(dec))).T

    dt = pmtime_utc_jyear - ref_year
    cosdec = np.cos(np.deg2rad(dec))
    dec = dec + dt * pmdec / (3600.0 * 1000.0)
    ra = ra + (dt * pmra / (3600.0 * 1000.0)) / cosdec
    parallax = np.atleast_1d(parallax)
    # AR discards parallax for scnd=True
    if scnd == True:
        parallax *= 0.0
    I = np.flatnonzero(parallax)
    if len(I):
        suntheta = 2.0 * np.pi * np.fmod(pmtime_utc_jyear - equinox_jyear, 1.0)
        # Finite differences on the unit sphere -- xyztoradec handles
        # points that are not exactly on the surface of the sphere.
        axis = np.deg2rad(axistilt)
        scale = parallax[I] / 1000.0 / arcsecperrad
        xyz = radectoxyz(ra[I], dec[I])
        xyz[:, 0] += scale * np.cos(suntheta)
        xyz[:, 1] += scale * np.sin(suntheta) * np.cos(axis)
        xyz[:, 2] += scale * np.sin(suntheta) * np.sin(axis)
        r, d = xyztoradec(xyz)
        ra[I] = r
        dec[I] = d
    return ra, dec
    
def force_finite_pm(
    d, pmra_key="PMRA", pmdec_key="PMDEC"
):
    """
    Replaces NaN PMRA, PMDEC by 0    
    
    Args:
        d: array with at least proper-motion columns
        pmra_key (optional, defaults to PMRA): column name for PMRA
        pmdec_key (optional, defaults to PMDEC): column name for PMDEC
        log (optional, defaults to Logger.get()): Logger object
        step (optional, defaults to ""): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start(optional, defaults to time()): start time for log (in seconds; output of time.time()

    Returns:
        d: same as input d, but NaN proper motions replaced by 0
        
    """
    for key in [pmra_key, pmdec_key]:
        keep = ~np.isfinite(d[key])
        if keep.sum() > 0:
            d[key][keep] = 0.0
            
    return d


def force_nonzero_refepoch(
    d,
    force_ref_epoch,
    ref_epoch_key="REF_EPOCH",
    pmra_key="PMRA",
    pmdec_key="PMDEC",
):
    """
    Replaces 0 by force_ref_epoch in ref_epoch
    
    Args:
        d: array with at least proper-motion columns
        force_ref_epoch: float, ref_epoch to replace 0 by
        ref_epoch_key (optional, defaults to REF_EPOCH): column name for the ref_epoch
        pmra_key (optional, defaults to PMRA): column name for PMRA
        pmdec_key (optional, defaults to PMDEC): column name for PMDEC
        log (optional, defaults to Logger.get()): Logger object
        step (optional, defaults to ""): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start(optional, defaults to time()): start time for log (in seconds; output of time.time()        
    Returns:
        d: same as input d, but 0 ref_epochs replaced by force_ref_epoch

    Notes:
        Will exit with error if ref_epoch=0, but pmra or pmdec != 0
        
    """
    keep = d[ref_epoch_key] == 0
    n = ((d[pmra_key][keep] != 0) | (d[pmra_key][keep] != 0)).sum()
    if n > 0:
        sys.exit(1)
    d[ref_epoch_key][keep] = force_ref_epoch
    return d


def update_nowradec(
    d,
    gaiadr,
    pmtime_utc_str,
    ra_key="RA",
    dec_key="DEC",
    pmra_key="PMRA",
    pmdec_key="PMDEC",
    parallax_key="PARALLAX",
    ref_epoch_key="REF_EPOCH",
    gaiag_key="GAIA_PHOT_G_MEAN_MAG",
    gaiaaen_key="GAIA_ASTROMETRIC_EXCESS_NOISE",
    scnd=False,
):
    """
    Update (RA, DEC, REF_EPOCH) using proper motion
    
    Args:
        d: array with at least proper-motion columns
        pmtime_utc_str: date to update position to (format: YYYY-MM-DDThh:mm:ss+00:00)
        gaiadr: Gaia dr ("dr2" or "edr3")
        ra_key (optional, defaults to RA): column name for RA
        dec_key (optional, defaults to DEC): column name for DEC
        pmra_key (optional, defaults to PMRA): column name for PMRA
        pmdec_key (optional, defaults to PMDEC): column name for PMDEC
        parallax_key (optional, defaults to PARALLAX): column name for PARALLAX
        ref_epoch_key (optional, defaults to REF_EPOCH): column name for the REF_EPOCH
        gaia_key (optional, defaults to GAIA_PHOT_G_MEAN_MAG): column name for Gaia g-mag
        gaiaaen_key (optional, defaults to GAIA_ASTROMETRIC_EXCESS_NOISE): column name for Gaia GAIA_ASTROMETRIC_EXCESS_NOISE
        scnd (optional, defaults to False): secondary target? (boolean);
              if False, update for REF_EPOCH>0 + AEN only
              if True, update for REF_EPOCH>0 + finite(PMRA,PMDEC) ; forces PARALLAX=0
        log (optional, defaults to Logger.get()): Logger object
        step (optional, defaults to ""): corresponding step, for fba_launch log recording
            (e.g. dotiles, dosky, dogfa, domtl, doscnd, dotoo)
        start(optional, defaults to time()): start time for log (in seconds; output of time.time()        
        
    Returns:
        d: same as input, but with RA, DEC updated to pmtime_utc_str

    Notes:
        REF_EPOCH is updated for *all* objects
    """
    # AR
    pmtime_utc = datetime.strptime(pmtime_utc_str, "%Y-%m-%dT%H:%M:%S%z")
    pmtime_utc_jyear = Time(pmtime_utc).jyear
    # AR computing positions at pmtime_utc_str using Gaia PMRA, PMDEC
    nowra, nowdec = get_nowradec(
        d[ra_key],
        d[dec_key],
        d[pmra_key],
        d[pmdec_key],
        d[parallax_key],
        d[ref_epoch_key],
        pmtime_utc_str,
        scnd=scnd,
    )
    if scnd == True:
        # AR secondary: REF_EPOCH>0
        keep = d["REF_EPOCH"] > 0
    else:
        # AR targets with REF_EPOCH>0 and passing the AEN criterion
        keep = d["REF_EPOCH"] > 0
        # AR gaia_psflike arguments changed at desitarget-0.58.0
        if desitarget.__version__ < "0.58.0":
            keep &= gaia_psflike(d[gaiag_key], d[gaiaaen_key])
        else:
            keep &= gaia_psflike(d[gaiag_key], d[gaiaaen_key], dr=gaiadr)
    # AR storing changes to report extrema in the log
    dra = nowra - d[ra_key]
    ddec = nowdec - d[dec_key]
    # AR updating positions to pmtime_utc_str for targets passing the AEN criterion
    d[ra_key][keep] = nowra[keep]
    d[dec_key][keep] = nowdec[keep]
    # AR updating REF_EPOCH for *all* objects (for PlateMaker)
    d[ref_epoch_key] = pmtime_utc_jyear
    return d
    

