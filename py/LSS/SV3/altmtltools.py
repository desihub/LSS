from desitarget import io 
from desitarget import mtl
import astropy.io.fits as pf
from astropy.table import Table,join,unique,vstack
from desitarget import mtl
from desitarget.targets import initial_priority_numobs
from desitarget.targetmask import obsconditions, obsmask
from desitarget.targetmask import desi_mask
from desitarget.mtl import get_mtl_dir, get_mtl_tile_file_name,get_mtl_ledger_format
import numpy as np
import matplotlib.pyplot as plt
import sys
from numpy import random as rand
import desitarget
import numpy.lib.recfunctions as rfn
import fitsio
import subprocess
import os
from desiutil.log import get_logger
from desitarget.mtl import get_zcat_dir, get_ztile_file_name, tiles_to_be_processed
from LSS.SV3.fatools import get_fba_fromnewmtl
import LSS.SV3.fatools as fatools
from desitarget.mtl import make_zcat,survey_data_model,update_ledger
from desitarget.targets import decode_targetid
log = get_logger()

print(desitarget.__file__)
print(mtl.__file__)
print(io.__file__)
print(pf.__file__)
print(fatools.__file__)

os.environ['DESIMODEL'] = '/global/common/software/desi/cori/desiconda/current/code/desimodel/master'

zcatdatamodel = np.array([], dtype=[
    ('RA', '>f8'), ('DEC', '>f8'), ('TARGETID', '>i8'),
    ('NUMOBS', '>i4'), ('Z', '>f8'), ('ZWARN', '>i8'), ('ZTILEID', '>i4')
    ])

mtltilefiledm = np.array([], dtype=[
    ('TILEID', '>i4'), ('TIMESTAMP', 'U25'),
    ('VERSION', 'U14'), ('PROGRAM', 'U6'), ('ZDATE', 'U8')
    ])
def createFAmap(FAReal, FAAlt, debug = False):
    TIDReal = FAReal['TARGETID']
    TIDAlt = FAAlt['TARGETID']
    FibReal = FAReal['FIBER']
    FibAlt = FAAlt['FIBER']
    
    Real2Alt = {}
    Alt2Real = {}
    if debug:
        inc1 = 0
        inc2 = 0
    negMisMatch = []
    for tr, fr in zip(TIDReal, FibReal):
        taMatch = TIDAlt[FibAlt == fr]
        assert(len(taMatch) == 1)
        if debug:
            try:
                assert(tr == taMatch[0])
            except:
                inc1+=1
        Real2Alt[tr] = taMatch[0]
    
    for ta, fa in zip(TIDAlt, FibAlt):
        trMatch = TIDReal[FibReal == fa]
        try:
            assert(len(trMatch) == 1)
        except:
            if ta < 0:
                #print('no match for negative ta {0}'.format(ta))
                negMisMatch.append(ta)
                continue
            else:
                print(ta)
                assert(0)
        if debug:
            try:
                assert(ta == trMatch[0])
            except:
                inc2+=1
        Alt2Real[ta] = trMatch[0]
    print('no matches for negative tas {0}'.format(negMisMatch))
    if debug:
        print(inc1)
        print(inc2)
    return Alt2Real, Real2Alt



def makeAlternateZCat(zcat, real2AltMap, alt2RealMap, debug = False):
    from collections import Counter
    zcatids = zcat['TARGETID']
    altZCat = Table(zcat)
    if debug:
        failures = 0
        negativeIDs = 0
    for n, i in zip(zcatids, range(len(zcatids))):
        cond = (n == zcatids)
        if debug and (n < 0):
            negativeIDs +=1   
        #try:
        altid = real2AltMap[n]
        altZCat['TARGETID'][i] = altid
    if debug:
        print('negIDs')
        print(negativeIDs)
        print('failures')
        print(failures)
        print('testctr')
    d =  Counter(altZCat['TARGETID'])  
    res = [ k for k, v in d.items() if v > 1]
    print(res)
    if len(res):
        print('how many pre dup cuts')
        print(zcatids.shape)
        cond2 = np.ones(zcatids.shape, dtype=bool)
        for i in res:
            print('test')
            print(np.sum(zcatids == i))
            cond2 = cond2 & (altcatids != i)
        print("how many post dup cuts")
        print(np.sum(cond2))
    else:
        print("supposedly, no duplicates")
    return altZCat

def checkMTLChanged(MTLFile1, MTLFile2):
    MTL1 = io.read_mtl_ledger(MTLFile1, unique = True)
    MTL2 = io.read_mtl_ledger(MTLFile2, unique = True)
    NDiff = 0
    NDiff2 = 0
    NDiff3 = 0
    for tar1 in MTL1:
        tar2 = MTL2[MTL2['TARGETID'] == tar1['TARGETID']]

        if tar1['NUMOBS'] != tar2['NUMOBS']:
            NDiff +=1

        if tar1['TIMESTAMP'] != tar2['TIMESTAMP']:
            NDiff2 +=1
            
        if tar1['SUBPRIORITY'] != tar2['SUBPRIORITY']:
            NDiff3 +=1

    print('Number targets with different NUMOBS')
    print(NDiff)
    print('Number targets with different TIMESTAMP')
    print(NDiff2)
    print('Number targets with different SUBPRIORITY')
    print(NDiff3)

def trimToMTL(notMTL, MTL, debug = False):
    # JL trims a target file, which possesses all of the information in an MTL, down
    # JL to the columns allowed in the MTL data model. 
    allNames = notMTL.dtype.names
    MTLNames = MTL.dtype.names
    for n in allNames:
        if n in MTLNames:
            if debug:
                print('allowed')
                print(n)
            continue
        else:
            if debug:
                print('killed')
                print(n)
            notMTL = rfn.drop_fields(notMTL, n)
    return notMTL



def initializeAlternateMTLs(initMTL, outputMTL, nAlt = 2, seed = 314159, obscon = 'DARK'):
    rand.seed(seed)

    allentries = Table.read(initMTL) 
    
    meta = allentries.meta
    firstTS = allentries[0]["TIMESTAMP"] 
    initialentries = allentries[allentries["TIMESTAMP"] == firstTS]
    subpriorsInit = initialentries["SUBPRIORITY"]
    for n in range(nAlt):
        outputMTLDir = outputMTL.format(n)
        subpriors = initialentries['SUBPRIORITY']
        shuffler = rand.permutation(len(subpriors))
        assert(np.std(subpriorsInit - subpriors[shuffler]) > 0.001)
        initialentries['SUBPRIORITY'] = subpriors[shuffler]
            
        io.write_mtl(outputMTLDir, initialentries, survey='sv3', obscon=obscon, extra=meta, nsidefile=meta['FILENSID'], hpxlist = [meta['FILEHPX']])
        
def loop_alt_ledger(obscon, survey='main', zcatdir=None, mtldir=None,
                altmtlbasedir=None, ndirs = 3, numobs_from_ledger=True, 
                secondary=False, singletile = None, debugOrig = False, 
                    getosubp = False, quickRestart = False, redoFA = False,
                    multiproc = False, nproc = None):
    """Execute full MTL loop, including reading files, updating ledgers.

    Parameters
    ----------
    obscon : :class:`str`
        A string matching ONE obscondition in the desitarget bitmask yaml
        file (i.e. in `desitarget.targetmask.obsconditions`), e.g. "DARK"
        Governs how priorities are set when merging targets.
    survey : :class:`str`, optional, defaults to "main"
        Used to look up the correct ledger, in combination with `obscon`.
        Options are ``'main'`` and ``'svX``' (where X is 1, 2, 3 etc.)
        for the main survey and different iterations of SV, respectively.
    zcatdir : :class:`str`, optional, defaults to ``None``
        Full path to the "daily" directory that hosts redshift catalogs.
        If this is ``None``, look up the redshift catalog directory from
        the $ZCAT_DIR environment variable.
    mtldir : :class:`str`, optional, defaults to ``None``
        Full path to the directory that hosts the MTL ledgers and the MTL
        tile file. If ``None``, then look up the MTL directory from the
        $MTL_DIR environment variable.
    altmtlbasedir : :class:`str`, optional, defaults to ``None``
        Formattable path to a directory that hosts alternate MTL ledgers  
        If ``None``, then look up the MTL directory from the
        $ALT_MTL_DIR environment variable. This will fail since that variable
        is not currently set in the desi code setup.
    ndirs : :class:`int`, optional, defaults to ``3``
        Number of alternate MTLs to process within altmtlbasedir
    numobs_from_ledger : :class:`bool`, optional, defaults to ``True``
        If ``True`` then inherit the number of observations so far from
        the ledger rather than expecting it to have a reasonable value
        in the `zcat.`
    secondary : :class:`bool`, optional, defaults to ``False``
        If ``True`` then process secondary targets instead of primaries
        for passed `survey` and `obscon`.
    quickRestart : :class:`bool`, optional, defaults to ``False``
        If ``True`` then copy original alternate MTLs from 
        altmtlbasedir/Univ*/survey/obscon/orig and  
    redoFA : :class:`bool`, optional, defaults to ``False``
        If ``True`` then automatically redo fiberassignment regardless of
        existence of fiberassign file in alternate fiberassign directory
    multiproc : :class:`bool`, optional, defaults to ``False``
        If ``True`` then run a single MTL update in a directory specified by
        nproc.
    nproc : :class:`int`, optional, defaults to None
        If multiproc is ``True`` this must be specified. Integer determines 
        directory of alternate MTLs to update.

    Returns
    -------
    :class:`str`
        The directory containing the ledger that was updated.
    :class:`str`
        The name of the MTL tile file that was updated.
    :class:`str`
        The name of the ZTILE file that was used to link TILEIDs to
        observing conditions and to determine if tiles were "done".
    :class:`~numpy.array`
        Information for the tiles that were processed.

    Notes
    -----
    - Assumes all of the relevant ledgers have already been made by,
      e.g., :func:`~LSS.SV3.altmtltools.initializeAlternateMTLs()`.
    """
    if debugOrig:
        iterloop = range(1)
        ndirs = 1
    elif multiproc:
        iterloop = range(nproc, nproc+1)
        ndirs = nproc+1
    else:
        iterloop = range(ndirs)

    if quickRestart:
        from shutil import copyfile, move
        from glob import glob as ls
        print('Restarting MTLs from copies')
        for nRestart in iterloop:
            altmtldirRestart = altmtlbasedir + '/Univ{0:03d}/'.format(nRestart)
            if os.path.exists(altmtldirRestart + 'mtl-done-tiles.ecsv'):
                move(altmtldirRestart + 'mtl-done-tiles.ecsv',altmtldirRestart + 'mtl-done-tiles.ecsv.old')
            restartMTLs = ls(altmtldirRestart +'/' + survey + '/' + obscon + '/' + '/orig/*')
            #print(altmtldirRestart +'/' + survey + '/' + obscon + '/' + '/orig/*')
            #print(restartMTLs)
            for fn in restartMTLs:
                copyfile(fn, altmtldirRestart +'/' + survey + '/' + obscon + '/' + fn.split('/')[-1])
    # ADM first grab all of the relevant files.
    # ADM grab the MTL directory (in case we're relying on $MTL_DIR).
    mtldir = get_mtl_dir(mtldir)
    # ADM construct the full path to the mtl tile file.
    mtltilefn = os.path.join(mtldir, get_mtl_tile_file_name(secondary=secondary))
    
    
    # ADM construct the relevant sub-directory for this survey and
    # ADM set of observing conditions..
    form = get_mtl_ledger_format()
    resolve = True
    msg = "running on {} ledger with obscon={} and survey={}"
    if secondary:
        log.info(msg.format("SECONDARY", obscon, survey))
        resolve = None
    else:
        log.info(msg.format("PRIMARY", obscon, survey))
    
    # ADM grab the zcat directory (in case we're relying on $ZCAT_DIR).
    zcatdir = get_zcat_dir(zcatdir)
    # ADM And contruct the associated ZTILE filename.
    ztilefn = os.path.join(zcatdir, get_ztile_file_name())
    
    if altmtlbasedir is None:
        print('This may automatically find the alt mtl dir in the future but fails now. Bye.')
        assert(0)

    for n in iterloop:
        print('')
        print('')
        print('')
        print('**NUMBER OF DIRECTORY*******')
        print(n)
        print('')
        print('')
        print('')
        print('')
        if debugOrig:
            altmtldir = altmtlbasedir
        else:
            altmtldir = altmtlbasedir + '/Univ{0:03d}/'.format(n)
        altmtltilefn = os.path.join(altmtldir, get_mtl_tile_file_name(secondary=secondary))

        althpdirname = io.find_target_files(altmtldir, flavor="mtl", resolve=resolve,
                                     survey=survey, obscon=obscon, ender=form)
        
        # ADM grab an array of tiles that are yet to be processed.
        tiles = tiles_to_be_processed(zcatdir, altmtltilefn, obscon, survey)

        # ADM stop if there are no tiles to process.
        if len(tiles) == 0:
            if n != ndirs - 1:
                continue
            else:
                return althpdirname, mtltilefn, ztilefn, tiles
        if not (singletile is None):
            tiles = tiles[tiles['TILEID'] == singletile]
        sorttiles = np.sort(tiles, order = 'ZDATE')

        for t in sorttiles:
            date = t['ZDATE']
            
            ts = str(t['TILEID']).zfill(6)
            print('processing tile {0}'.format(ts))
            FAOrigName = '/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'+ts[:3]+'/fiberassign-'+ts+'.fits.gz'
            fhtOrig = fitsio.read_header(FAOrigName)
            fadate = fhtOrig['RUNDATE']
            fadate = ''.join(fadate.split('T')[0].split('-'))

            fbadirbase = altmtldir + '/fa/' + survey.upper() +  '/' + fadate + '/'
            if getosubp:
                FAAltName = altmtldir + '/fa/' + survey.upper() +  '/' + fadate + '/orig/fba-' + ts+ '.fits'
                fbadir = altmtldir + '/fa/' + survey.upper() +  '/' + fadate + '/orig/'
            else:
                
                FAAltName = altmtldir + '/fa/' + survey.upper() +  '/' + fadate + '/fba-' + ts+ '.fits'
                fbadir = fbadirbase
            print('PYTHONPATH preFA')
            print(os.getenv('PYTHONPATH'))
            print('PATH preFA')
            print(os.getenv('PATH'))
            if True or redoFA or (not os.path.exists(FAAltName)):
                try:
                    get_fba_fromnewmtl(ts,mtldir=altmtldir + survey.lower() + '/',outdir=fbadirbase, getosubp = getosubp)
                    command_run = (['bash', fbadir + 'fa-' + ts + '.sh'])
                    result = subprocess.run(command_run, capture_output = True, shell = True)
                except:
                    print('get_fba_fromnewmtl failed in univ {0}. Trying again'.format(n))
                    print('mtldir = {0}'.format(altmtldir + survey.lower()))
                    print('fbadirbase = {0}'.format(fbadirbase))
                    get_fba_fromnewmtl(ts,mtldir=altmtldir + survey.lower() + '/',outdir=fbadirbase, getosubp = getosubp)
                    command_run = (['bash', fbadir + 'fa-' + ts + '.sh'])
                    result = subprocess.run(command_run, capture_output = True, shell = True)
            print('PYTHONPATH postFA')
            print(os.getenv('PYTHONPATH'))
            print('PATH postFA')
            print(os.getenv('PATH'))
            OrigFA = pf.open(FAOrigName)[1].data
            AltFA = pf.open(FAAltName)[1].data
            
            # ADM create the catalog of updated redshifts.
            zcat = make_zcat(zcatdir, [t], obscon, survey)

            # ADM insist that for an MTL loop with real observations, the zcat
            # ADM must conform to the data model. In particular, it must include
            # ADM ZTILEID, and other columns addes for the Main Survey. These
            # ADM columns may not be needed for non-ledger simulations.
            # ADM Note that the data model differs with survey type.
            zcatdm = survey_data_model(zcatdatamodel, survey=survey)
            if zcat.dtype.descr != zcatdm.dtype.descr:
                msg = "zcat data model must be {} not {}!".format(
                    zcatdm.dtype.descr, zcat.dtype.descr)
                log.critical(msg)
                raise ValueError(msg)
            # ADM useful to know how many targets were updated.
            _, _, _, _, sky, _ = decode_targetid(zcat["TARGETID"])
            ntargs, nsky = np.sum(sky == 0), np.sum(sky)
            msg = "Update state for {} targets".format(ntargs)
            msg += " (the zcats also contain {} skies with +ve TARGETIDs)".format(nsky)
            log.info(msg)
            
            

            
            A2RMap, R2AMap = createFAmap(OrigFA, AltFA)

            altZCat = makeAlternateZCat(zcat, R2AMap, A2RMap)

            # ADM update the appropriate ledger.
            print('althpdirname = {0}'.format(althpdirname))

            update_ledger(althpdirname, altZCat, obscon=obscon.upper(),
                          numobs_from_ledger=numobs_from_ledger)
            print(os.getenv('PYTHONPATH'))
            # JL write single processed tile to the MTL tile file. Since adding
            # JL alternate fiber assignment takes longer than simply updating
            # JL mtls, it is possible/likely to be interrupted in the middle 
            # JL and will need to restart update using this information
            io.write_mtl_tile_file(altmtltilefn,np.array([t]))
        # ADM for the main survey "holding pen" method, ensure the TIMESTAMP
        # ADM in the mtl-done-tiles file is always later than in the ledgers.
        if survey == "main":
            sleep(1)
            tiles["TIMESTAMP"] = get_utc_date(survey=survey)

    return althpdirname, altmtltilefn, ztilefn, tiles