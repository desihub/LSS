import os
from astropy.table import Table, join, vstack
import fitsio
import numpy as np
from desitarget.io import read_targets_in_tiles
import desimodel.focalplane

#Function that test if a directory exist or not
def test_dir(value):
    if not os.path.exists(value):
        try:
            os.makedirs(value, 0o755)
            print('made %s'%value)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

#Function that reads tile information
def read_info_tiles(tilef, mtld, pr):
    if os.path.isfile(tilef):
        ta = Table.read(tilef)
    else:
        #construct a table with the needed tile information
        tilel = []
        ral = []
        decl = []
        mtlt = []
        fal = []
        obsl = []
        pl = []
        hal = []
        #for tile,pro in zip(mtld['TILEID'],mtld['PROGRAM']):
        for tile in mtld['TILEID']:
                ts = str(tile).zfill(6)
                fht = fitsio.read_header('/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'+ts[:3]+'/fiberassign-'+ts+'.fits.gz')
                tilel.append(tile)
                ral.append(fht['TILERA'])
                decl.append(fht['TILEDEC'])
                mtlt.append(fht['MTLTIME'])
                fal.append(fht['RUNDATE'])
                obsl.append(fht['FIELDROT'])
                hal.append(fht['FA_HA'])
                #pl.append(pro)
                pl.append(pr)
        ta = Table()
        ta['TILEID'] = tilel
        ta['RA'] = ral
        ta['DEC'] = decl
        ta['MTLTIME'] = mtlt
        ta['RUNDATE'] = fal
        ta['FIELDROT'] = obsl
        ta['PROGRAM'] = pl
        ta['FA_HA'] = hal
        #if pd == 'dark':
        ta['OBSCONDITIONS'] = 15
        ta['IN_DESI'] = 1
        ta.write(tilef, format='fits')
    return ta

#CREATE A DICTIONARY WITH TILEID AND THE DIRECTORY OF THE ALTMTL FBA RUN
##############################################################################################
def create_tile_altmtldir(mockrea, id_, ta):
    list_runFA = {}
    #'/global/cscratch1/sd/acarnero/alt_mtls_masterScriptTest_256dirs_rea{MOCKREA}/Univ{UNIV}/mtl-done-tiles.ecsv'.format(MOCKREA=mockrea, UNIV=id_))
    for tile in ta['TILEID']:
        ts = str(tile).zfill(6)
        faf_d = '/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'+ts[:3]+'/fiberassign-'+ts+'.fits.gz'
        fht = fitsio.read_header(faf_d)
        stamp = fht['RUNDATE'].split('T')[0].replace('-','')
        list_runFA[tile] = stamp
    return list_runFA

#Concatenate fiberassigment FASSIGN hdu and match to masterTarget
def combtile_specmock(tiles, fbaRun, list_runFA, masterTarget, fout=''):

    print('***************************************************************************************************')
    print('*Entering my combtile_specmock to concatenate fiberassigment FASSIGN hdu and match to masterTarget*')
    print('***************************************************************************************************')
    s = 0
    n = 0
    tmask = np.ones(len(tiles)).astype('bool')

    tars = Table.read(masterTarget)

    for tile in tiles[tmask]['TILEID']:
        print('reading combtile_speckmock', tile)
        faf = fbaRun.format(stamp=list_runFA[tile], ts = str(tile).zfill(6))

        tt = Table.read(faf,hdu='FASSIGN')
        tt['TILEID'] = tile
        tt = join(tt, tars, keys=['TARGETID'])

        if s == 0:
            ttn = tt
            s = 1
        else:
            ttn = vstack([ttn,tt], metadata_conflicts = 'silent')

        print('----')
        ttn.sort('TARGETID')
        n += 1

    ttn = add_tilelocid(ttn)
    ttn.write(fout,format='fits', overwrite=True)
    print('Ending combtile_speckmock, read ', n, ' tiles')
    print('---------------------------------------------------------------------------------------------------')
    return ttn

#Add TILELOCID TO Table
def add_tilelocid(input_arr):
    input_arr['TILELOCID'] = 10000*input_arr['TILEID'] + input_arr['LOCATION']
    return input_arr



def combtiles_wdup(tiles, mdir='',fout='', mtl_done=None, tarcol=['RA','DEC','TARGETID','SV3_DESI_TARGET','SV3_BGS_TARGET','SV3_MWS_TARGET','SUBPRIORITY','PRIORITY_INIT','TARGET_STATE','TIMESTAMP','ZWARN','PRIORITY'], isodate=None, univ='001', mockrea='000'):

    print('***************************************************************************************************')
    print('*Entering my combtiles_wdup to concatenate fiberassigment FAVAIL hdu and match to masterTarget*')
    print('***************************************************************************************************')

    s = 0
    n = 0
    tmask = np.ones(len(tiles)).astype('bool')

    if mtl_done is not None:
        infp = Table.read(mtl_done)
        for tile in tiles[tmask]['TILEID']:
            ts = str(tile).zfill(6)
            faf_d = '/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'+ts[:3]+'/fiberassign-'+ts+'.fits.gz'
            fht = fitsio.read_header(faf_d)
            wt = tiles['TILEID'] == tile
            tars = read_targets_in_tiles(mdir, tiles[wt], mtl=True, isodate=isodate, columns=tarcol)
           
            zdate = infp[(infp['TILEID']==tile)]['ZDATE'][0]
            stamp = fht['RUNDATE'].split('T')[0].replace('-','')
            faf = os.path.join('/global/cscratch1/sd/acarnero/alt_mtls_masterScriptTest_256dirs_rea{MOCKREA}/Univ{UNIV}/fa/SV3'.format(UNIV=univ, MOCKREA=mockrea),stamp,'fba-'+ts+'.fits')
            
            tt = Table.read(faf,hdu='FAVAIL')
            tars = join(tars,tt,keys=['TARGETID'])
            tars['TILEID'] = tile
            tars['ZWARN'].name = 'ZWARN_MTL'
            if s == 0:
                tarsn = tars
                s = 1
            else:
                tarsn = vstack([tarsn,tars],metadata_conflicts='silent')
            tarsn.sort('TARGETID')
            n += 1
            print(tile,n,len(tiles[tmask]),len(tarsn))

    else:
        for tile in tiles[tmask]['TILEID']:
            ts = str(tile).zfill(6)
            faf = '/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'+ts[:3]+'/fiberassign-'+ts+'.fits.gz'
            fht = fitsio.read_header(faf)
            wt = tiles['TILEID'] == tile
            #tars = read_targets_in_tiles(mdir,tiles[wt],mtl=True,isodate=fht['MTLTIME'])
            tars = read_targets_in_tiles(mdir,tiles[wt],mtl=True,isodate=fht['MTLTIME'],columns=tarcol)
            #tars.keep_columns(tarcols)
            #tars = tars[[b for b in tarcol]]

            tt = Table.read(faf,hdu='POTENTIAL_ASSIGNMENTS')
            tars = join(tars,tt,keys=['TARGETID'])
            tars['TILEID'] = tile
            tars['ZWARN'].name = 'ZWARN_MTL'
            if s == 0:
                tarsn = tars
                s = 1
            else:
                tarsn = vstack([tarsn,tars],metadata_conflicts='silent')
            tarsn.sort('TARGETID')
            n += 1
            print(tile,n,len(tiles[tmask]),len(tarsn))
    #tarsn.keep_columns(['RA','DEC','TARGETID','SV3_DESI_TARGET','SV3_BGS_TARGET','SV3_MWS_TARGET','SUBPRIORITY',
    #    'PRIORITY_INIT','TARGET_STATE','TIMESTAMP','ZWARN_MTL','PRIORITY','FIBER','LOCATION','TILEID'])
    
    tarsn = add_tilelocid(tarsn)
    tarsn.write(fout,format='fits', overwrite=True)

    return tarsn

#count_tiles_better(real_specf, filename_tarspecwdup, pdir, specrel=specrel)
def count_tiles_better(real_specf, dr, pdir, specrel='fuji', rann=0, fibcol='COADD_FIBERSTATUS'):

#def count_tiles_better(fs,dr,pd,rann=0,specrel='daily',fibcol='COADD_FIBERSTATUS'):
    '''
    from files with duplicates that have already been sorted by targetid, quickly go
    through and get the multi-tile information
    dr is either 'dat' or 'ran'
    returns file with TARGETID,NTILE,TILES,TILELOCIDS
    '''

    #fs = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/'+specrel+'/datcomb_'+pd+'_specwdup_Alltiles.fits')
    #wf = fs['FIBERSTATUS'] == 0
    real_specf = add_tilelocid(real_specf)
    wf = real_specf[fibcol] == 0
    gtl = np.unique(real_specf['TILELOCID'][wf])

    if dr == 'dat':
        fj = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/'+specrel+'/datcomb_'+pd+'_tarspecwdup_Alltiles.fits')
        #outf = '/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/datcomb_'+pd+'ntileinfo.fits'
    elif dr == 'ran':
        fj = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/'+specrel+'/rancomb_'+str(rann)+pd+'wdupspec_Alltiles.fits')
        #outf = '/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/random'+str(rann)+'/rancomb_'+pd+'ntileinfo.fits'
    else:
        fj = fitsio.read(dr)

    wg = np.isin(fj['TILELOCID'],gtl)
    fjg = fj[wg]

    tids = np.unique(fjg['TARGETID'])
    nloc = []#np.zeros(len(np.unique(f['TARGETID'])))
    nt = []
    tl = []
    tli = []
    ti = 0
    i = 0
    while i < len(fjg):
        tls  = []
        tlis = []
        nli = 0

        while fjg[i]['TARGETID'] == tids[ti]:
            nli += 1
            tls.append(fjg[i]['TILEID'])
            tlis.append(fjg[i]['TILELOCID'])
            i += 1
            if i == len(fjg):
                break
        nloc.append(nli)
        tlsu = np.unique(tls)
        tlisu = np.unique(tlis)
        nt.append(len(tlsu))
        tl.append("-".join(tlsu.astype(str)))
        tli.append("-".join(tlisu.astype(str)))

        if ti%100000 == 0:
            print(ti)
        ti += 1
    tc = Table()
    tc['TARGETID'] = tids
    tc['NTILE'] = nt
    tc['TILES'] = tl
    tc['TILELOCIDS'] = tli

    return tc

#Reading the target file, create a file per tile, with duplicates of the targets
def randomtiles_allSV3(tiles, mytargets, directory_output='.'):
    '''
    tiles should be a table containing the relevant info
    '''
    trad = desimodel.focalplane.get_tile_radius_deg()*1.1 #make 10% greater just in case
    
    rt = fitsio.read(mytargets)
    print('loaded random file', mytargets)
    for i in range(0,len(tiles)):

        #print('length of tile file is (expected to be 1):'+str(len(tiles)))
        tile = tiles['TILEID'][i]
        fname = os.path.join(directory_output, 'tilenofa-'+str(tile)+'.fits')
        if os.path.isfile(fname):
            print(fname +' already exists')
        else:
            tdec = tiles['DEC'][i]
            decmin = tdec - trad
            decmax = tdec + trad
            wdec = (rt['DEC'] > decmin) & (rt['DEC'] < decmax)
            print(len(rt[wdec]))
            inds = desimodel.footprint.find_points_radec(tiles['RA'][i], tdec,rt[wdec]['RA'], rt[wdec]['DEC'])
            print('got indexes')
            rtw = rt[wdec][inds]
            rmtl = Table(rtw)
    
            '''
            tot += len(rmtl['TARGETID'])
            #rmtl['TARGETID'] = np.arange(len(rmtl))
            print(len(rmtl['TARGETID'])) #checking this column is there

            rmtl['SV3_SCND_TARGET'] = np.zeros(len(rmtl),dtype=int)
            rmtl['SV3_BGS_TARGET'] = np.zeros(len(rmtl),dtype=int)
            rmtl['SV3_MWS_TARGET'] = np.zeros(len(rmtl),dtype=int)

    
            rmtl['DESI_TARGET'] = np.ones(len(rmtl),dtype=int)*2
            rmtl['SV3_DESI_TARGET'] = np.ones(len(rmtl),dtype=int)*2
            rmtl['NUMOBS_INIT'] = np.zeros(len(rmtl),dtype=int)
            rmtl['NUMOBS_MORE'] = np.ones(len(rmtl),dtype=int)
            rmtl['PRIORITY'] = np.ones(len(rmtl),dtype=int)*3400
            rmtl['OBSCONDITIONS'] = np.ones(len(rmtl),dtype=int)*516#tiles['OBSCONDITIONS'][i]
            rmtl['SUBPRIORITY'] = np.random.random(len(rmtl))
            

            colums_to_remove = ['ZWARN', 'NZ', 'TRUEZ', 'RSDZ', 'BGS_TARGET', 'MWS_TARGET', 'DESI_TARGET']
            for col in colums_to_remove:
                if col in rmtl.keys():
                    del rmtl[col]
            '''
            rmtl.write(fname,format='fits', overwrite=True)
            print('added columns, wrote to '+fname)


def combran_wdup(tiles,rann,randir,tp,sv3dir,specf,keepcols=[]):

    s = 0
    td = 0
    #tiles.sort('ZDATE')
    print(len(tiles))
    #delcols = ['DESI_TARGET','BGS_TARGET','MWS_TARGET','SUBPRIORITY','OBSCONDITIONS','PRIORITY_INIT',\
    #'NUMOBS_INIT','SCND_TARGET','NUMOBS_MORE','NUMOBS','Z','ZWARN','TARGET_STATE','TIMESTAMP','VERSION','PRIORITY']

    outf = os.path.join(randir+str(rann),'rancomb_'+tp+'wdup_Alltiles.fits')

#    if os.path.isfile(outf):
#        fgu = Table.read(outf)
        #tarsn.keep_columns(['RA','DEC','TARGETID''LOCATION','FIBER','TILEID'])
#        s = 1
#        tdone = np.unique(fgu['TILEID'])
#        tmask = ~np.isin(tiles['TILEID'],tdone)
#    else:
    
    tmask = np.ones(len(tiles)).astype('bool')
    for tile in tiles[tmask]['TILEID']:
        ffa = os.path.join(randir+str(rann),'fba-'+str(tile).zfill(6)+'.fits')
        ffna = os.path.join(randir+str(rann),'tilenofa-'+str(tile)+'.fits')
        if os.path.isfile(ffa):
            fa = Table.read(ffa,hdu='FAVAIL')

            ffna = Table.read(ffna)
            fgun = join(fa,ffna,keys=['TARGETID'])
            #fgun.remove_columns(delcols)

            td += 1
            fgun['TILEID'] = int(tile)
            fgun.keep_columns(['RA','DEC','TARGETID','LOCATION','FIBER','TILEID', 'SV3_DESI_TARGET'])
            if s == 0:
                fgu = fgun
                s = 1
            else:
                fgu = vstack([fgu,fgun],metadata_conflicts='silent')
            fgu.sort('TARGETID')
            print(tile,td, len(tiles), len(fgun),len(fgu))
        else:
            print('did not find '+ffa)

    if len(tiles[tmask]['TILEID']) > 0:
        fgu.write(outf,format='fits', overwrite=True)
    #specf = Table.read(sv3dir+'datcomb_'+tp+'_specwdup_Alltiles.fits')

    #    specf = add_tilelocid(specf) #['TILELOCID'] = 10000*specf['TILEID'] +specf['LOCATION']
    specf.keep_columns(keepcols)
    #specf.keep_columns(['ZWARN','LOCATION','TILEID','TILELOCID','FIBERSTATUS','FIBERASSIGN_X','FIBERASSIGN_Y','PRIORITY','DELTA_X','DELTA_Y','EXPTIME','PSF_TO_FIBER_SPECFLUX','TSNR2_ELG_B','TSNR2_LYA_B','TSNR2_BGS_B','TSNR2_QSO_B','TSNR2_LRG_B','TSNR2_ELG_R','TSNR2_LYA_R','TSNR2_BGS_R','TSNR2_QSO_R','TSNR2_LRG_R','TSNR2_ELG_Z','TSNR2_LYA_Z','TSNR2_BGS_Z','TSNR2_QSO_Z','TSNR2_LRG_Z','TSNR2_ELG','TSNR2_LYA','TSNR2_BGS','TSNR2_QSO','TSNR2_LRG'])
    fgu = join(fgu,specf,keys=['LOCATION','TILEID','FIBER'])
    fgu.sort('TARGETID')
    outf = os.path.join(sv3dir,'rancomb_'+str(rann)+tp+'wdupspec_Alltiles.fits')
    fgu.write(outf,format='fits', overwrite=True)
    '''
    print(outf)
    if os.path.isfile('tmp.fits'):
        os.system('rm tmp.fits')
    fd = fitsio.FITS('tmp.fits', "rw")
    fd.write(np.array(fgu),extname='POTENTIAL_ASSINGMENT')
    fd['POTENTIAL_ASSINGMENT'].write_comment("concatenation of SV3 POTENTIAL_ASSIGNMENT information for randoms, joined to columns in targe files")
    fd['POTENTIAL_ASSINGMENT'].write_history("updated on "+datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
    fd.close()
    os.system('mv tmp.fits '+fout)
    '''
    #fgu.write(outf,format='fits', overwrite=True)
