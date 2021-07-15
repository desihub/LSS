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
log = get_logger()


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
                secondary=False, singletile = None, debugOrig = False, getosubp = False):
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
      e.g., :func:`~desitarget.mtl.make_ledger()`.
    """
    # ADM first grab all of the relevant files.
    # ADM grab the MTL directory (in case we're relying on $MTL_DIR).
    mtldir = get_mtl_dir(mtldir)
    # ADM construct the full path to the mtl tile file.
    mtltilefn = os.path.join(mtldir, get_mtl_tile_file_name(secondary=secondary))
    
    print(mtldir)
    print(mtltilefn)
    
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
    
    print(zcatdir)
    print(ztilefn)
    
    if altmtlbasedir is None:
        print('This may automatically find the alt mtl dir in the future but fails now. Bye.')
        assert(0)
    if debugOrig:
        ndirs = 1
    for n in range(ndirs):
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
            print(date)
            print(ts)
            FAOrigName = '/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'+ts[:3]+'/fiberassign-'+ts+'.fits.gz'
            fhtOrig = fht = fitsio.read_header(FAOrigName)
            fadate = fhtOrig['RUNDATE']
            fbadirbase = altmtldir + '/fa/' + survey.upper() +  '/' + date + '/'
            if getosubp:
                FAAltName = altmtldir + '/fa/' + survey.upper() +  '/' + date + '/orig/fba-' + ts+ '.fits'
                fbadir = altmtldir + '/fa/' + survey.upper() +  '/' + date + '/orig/'
            else:
                FAAltName = altmtldir + '/fa/' + survey.upper() +  '/' + date + '/fba-' + ts+ '.fits'
                fbadir = fbadirbase
            if not os.path.exists(FAAltName):
                
                get_fba_fromnewmtl(ts,mtldir=altmtldir + survey.lower() + '/',outdir=fbadirbase, getosubp = getosubp)

                command_run = (['bash', fbadir + 'fa-' + ts + '.sh'])
                result = subprocess.run(command_run, capture_output = True)
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
            update_ledger(althpdirname, altZCat, obscon=obscon.upper(),
                          numobs_from_ledger=numobs_from_ledger)

        # ADM for the main survey "holding pen" method, ensure the TIMESTAMP
        # ADM in the mtl-done-tiles file is always later than in the ledgers.
        if survey == "main":
            sleep(1)
            tiles["TIMESTAMP"] = get_utc_date(survey=survey)

        # ADM write the processed tiles to the MTL tile file.
        io.write_mtl_tile_file(altmtltilefn, tiles)

    return althpdirname, altmtltilefn, ztilefn, tiles