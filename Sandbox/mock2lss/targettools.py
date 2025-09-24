from astropy.table import Table, join, unique
import numpy as np
from LSS import common_tools as common
import os

def mkfulldat(zf, imbits, tdir, tp, bit, outf, ftiles, azf='', desitarg='SV3_DESI_TARGET', specver='guadalupe', notqso='', qsobit=4, bitweightfile=None):
    '''
    zf is the name of the file containing all of the combined spec and target info compiled already
    imbits is the list of imaging mask bits to mask out
    tdir is the directory for the targets
    tp is the target type
    bit is the SV3_{type}_MASK bit to use for select the correct target type
    outf is the full path + name for the output file
    ftiles is the name of the file containing information on, e.g., how many tiles each target was available on
    azf is the file name for OII flux info (relevant for ELGs only)
    desitarg is the column to use for the target type cut (all use SV3_DESI_TARGET except BGS_BRIGHT)
    specver is the version of the pipeline used for the redshift info; only 'daily' exists for now
    '''
    
    
    #from desitarget.mtl import inflate_ledger
    if tp[:3] == 'BGS' or tp[:3] == 'MWS':
        pd = 'bright'        
        tscol = 'TSNR2_BGS'
    else:    
        pd = 'dark'
        tscol = 'TSNR2_ELG'
    #load in the appropriate dark/bright combined spec file and use to denote the tileid + location that had good observations:
    #fs = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/'+specver+'/datcomb_'+pd+'_specwdup_Alltiles.fits')

    if specver == 'daily':
        #fbcol = 'FIBERSTATUS'
        print('no longer supported')
        return False

    #read in the big combined data file
    dz = Table.read(zf) 
    dz['RSDZ'].name = 'Z'
    dz['Z'][(dz['Z']==1.e20)] = 999999
    dz[tscol][(dz[tscol]==1.e20)] = 999999

    #find the rows that satisfy the target type
    wtype = ((dz[desitarg] & bit) > 0)
    if notqso == 'notqso':
        print('removing QSO targets')
        wtype &= ((dz[desitarg] & qsobit) == 0)
    #find the rows that are 'good' tilelocid
    
    print(len(dz[wtype]))
    print('length after selecting type '+str(len(dz)))
    #print(len(dz[wg]))
    #down-select to target type of interest and good tilelocid
    dz = dz[wtype]#&wg]
    
    wz = dz['ZWARN'] != 999999 #this is what the null column becomes
    wz &= dz['ZWARN']*0 == 0 #just in case of nans
    wz &= dz['COADD_FIBERSTATUS'] == 0
    fs = dz[wz]
    print('number of good obs '+str(len(fs)))
    #fs = common.cut_specdat(dz)
    gtl = np.unique(fs['TILELOCID'])
    wg = np.isin(dz['TILELOCID'],gtl)
    dz['GOODHARDLOC'] = np.zeros(len(dz)).astype('bool')
    dz['GOODHARDLOC'][wg] = 1

    #THIS IS NOT NEEDED FOR MOCKS
    '''
    print('joining to full imaging')
    ftar = Table.read('/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/'+pd+'_targets.fits')
    ftar.keep_columns(['TARGETID','EBV','FLUX_G','FLUX_R','FLUX_Z','FLUX_IVAR_G','FLUX_IVAR_R','FLUX_IVAR_Z','MW_TRANSMISSION_G','MW_TRANSMISSION_R',\
            'MW_TRANSMISSION_Z','FRACFLUX_G','FRACFLUX_R','FRACFLUX_Z','FRACMASKED_G','FRACMASKED_R','FRACMASKED_Z','FRACIN_G','FRACIN_R',\
            'FRACIN_Z','NOBS_G','NOBS_R','NOBS_Z','PSFDEPTH_G','PSFDEPTH_R','PSFDEPTH_Z','GALDEPTH_G','GALDEPTH_R','GALDEPTH_Z','FLUX_W1',\
            'FLUX_W2','FLUX_IVAR_W1','FLUX_IVAR_W2','MW_TRANSMISSION_W1','MW_TRANSMISSION_W2','ALLMASK_G','ALLMASK_R','ALLMASK_Z','FIBERFLUX_G',\
            'FIBERFLUX_R','FIBERFLUX_Z','FIBERTOTFLUX_G','FIBERTOTFLUX_R','FIBERTOTFLUX_Z','WISEMASK_W1','WISEMASK_W2','MASKBITS',\
            'RELEASE','BRICKID','BRICKNAME','BRICK_OBJID','MORPHTYPE','PHOTSYS','SHAPE_R'])
    dz = join(dz,ftar,keys=['TARGETID'])
    print('length after join to full targets (should be same) '+str(len(dz)))
    '''

    #apply imaging veto mask
    #NOT APPLY IN MOCKS    dz = common.cutphotmask(dz,imbits)
    
    #load in file with information about where repeats occurred and join it
    dtl = Table.read(ftiles)
    dtl.keep_columns(['TARGETID','NTILE','TILES','TILELOCIDS'])
    dz = join(dz,dtl,keys='TARGETID')
    
    #find the rows where we have spectroscopic observations
    wz = dz['ZWARN'] != 999999 #this is what the null column becomes
    wz &= dz['ZWARN']*0 == 0 #just in case of nans
    
    
    #mark them as having LOCATION_ASSIGNED
    dz['LOCATION_ASSIGNED'] = np.zeros(len(dz)).astype('bool')
    dz['LOCATION_ASSIGNED'][wz] = 1
    #find the TILELOCID that were assigned and mark them as so
    tlids = np.unique(dz['TILELOCID'][wz])
    #test that all with goodhardloc have z
    gin = np.isin(gtl,tlids)
    print('gtl in tlids, should be all',np.sum(gin),len(gtl))
    wtl = np.isin(dz['TILELOCID'],tlids)
    dz['TILELOCID_ASSIGNED'] = 0
    dz['TILELOCID_ASSIGNED'][wtl] = 1
    print('number of unique targets at assigned tilelocid:')
    print(len(np.unique(dz[wtl]['TARGETID'])))

    #get OII flux info for ELGs
    if tp == 'ELG' or tp == 'ELG_HIP':
        if azf != '':
            '''
            arz = fitsio.read(azf,columns=['TARGETID','LOCATION','TILEID','OII_FLUX','OII_FLUX_IVAR','SUBSET','DELTACHI2'])
            st = []
            for i in range(0,len(arz)):
                st.append(arz['SUBSET'][i][:4])
            st = np.array(st)
            #wg = arz[fbcol] == 0
            wg = st == "thru"
            arz = arz[wg]
            o2c = np.log10(arz['OII_FLUX'] * np.sqrt(arz['OII_FLUX_IVAR']))+0.2*np.log10(arz['DELTACHI2'])
            w = (o2c*0) != 0
            w |= arz['OII_FLUX'] < 0
            o2c[w] = -20
            #arz.keep_columns(['TARGETID','LOCATION','TILEID','o2c','OII_FLUX','OII_SIGMA'])#,'Z','ZWARN','TSNR2_ELG'])    
            arz = Table(arz)
            arz['o2c'] = o2c
            dz = join(dz,arz,keys=['TARGETID','LOCATION','TILEID'],join_type='left',uniq_col_name='{col_name}{table_name}',table_names=['', '_OII'])
   
            dz.remove_columns(['SUBSET','DELTACHI2_OII'])
            print('check length after merge with OII strength file:' +str(len(dz)))
            #join changes order, so get wz again
            '''
            wz = dz['ZWARN'] != 999999 #this is what the null column becomes
            wz &= dz['ZWARN']*0 == 0 #just in case of nans

    if tp[:3] == 'QSO':
        if azf != '':
            '''
            arz = Table(fitsio.read(azf))
            arz.keep_columns(['TARGETID','LOCATION','TILEID','Z','ZERR','Z_QN'])
            arz['TILEID'] = arz['TILEID'].astype(int)
            print(arz.dtype.names)
            #arz['TILE'].name = 'TILEID'
            dz = join(dz,arz,keys=['TARGETID','TILEID','LOCATION'],join_type='left',uniq_col_name='{col_name}{table_name}',table_names=['','_QF'])
            dz['Z'].name = 'Z_RR' #rename the original redrock redshifts
            dz['Z_QF'].name = 'Z' #the redshifts from the quasar file should be used instead
            #join changes order, so get wz again
            '''
            wz = dz['ZWARN'] != 999999 #this is what the null column becomes
            wz &= dz['ZWARN']*0 == 0 #just in case of nans

    
    #sort and then cut to unique targetid; sort prioritizes observed targets and then TSNR2
    wnts = dz[tscol]*0 != 0
    wnts |= dz[tscol] == 999999
    dz[tscol][wnts] = 0
    print(np.max(dz[tscol]))
    dz['sort'] = dz['LOCATION_ASSIGNED']*np.clip(dz[tscol],0,200)*dz['GOODHARDLOC']+dz['TILELOCID_ASSIGNED']*dz['GOODHARDLOC']+dz['GOODHARDLOC']
    print('sort min/max',np.min(dz['sort']),np.max(dz['sort']))
    dz.sort('sort')
    dz = unique(dz,keys=['TARGETID'],keep='last')

    '''
    if tp == 'ELG' or tp == 'ELG_HIP':
        print('number of masked oII row (hopefully matches number not assigned) '+ str(np.sum(dz['o2c'].mask)))
    if tp == 'QSO':
        print('number of good z according to qso file '+str(len(dz)-np.sum(dz['Z'].mask)))
    '''
#THERE IS NO NANS IN MOCK    dz['Z'] = dz['Z'].filled(999999)
    selm = dz['Z'] == 999999
    print('999999s for Z',len(dz[selm]))
    print('length after cutting to unique targetid '+str(len(dz)))
    print('LOCATION_ASSIGNED numbers')
    print(np.unique(dz['LOCATION_ASSIGNED'],return_counts=True))
   
    print('TILELOCID_ASSIGNED numbers')
    print(np.unique(dz['TILELOCID_ASSIGNED'],return_counts=True))

    
    #get completeness based on unique sets of tiles "comp_tile"
    tll,compa = common.comp_tile(dz)
    comp_dicta = dict(zip(tll, compa))
    fcompa = []
    for tl in dz['TILES']:
        fcompa.append(comp_dicta[tl]) 
    dz['COMP_TILE'] = np.array(fcompa)
    wc0 = dz['COMP_TILE'] == 0
    print('number of targets in 0 completeness regions '+str(len(dz[wc0])))   
    
    #write out comp_tile info
    tll = np.array(tll).astype(dz['TILES'].dtype)
    co = Table()
    co['TILES'] = tll
    co['COMP_TILE'] = compa
    cof = outf.strip('_full_noveto.dat.fits')+'_comp_tile.fits'
    print('writing comp_tile completeness to '+cof)
    co.write(cof,overwrite=True,format='fits')    

    #get counts at unique TILELOCID
    locl,nlocl = np.unique(dz['TILELOCID'],return_counts=True)
    #do same after cutting to only the data with location_assigned
    wz = dz['LOCATION_ASSIGNED'] == 1
    dzz = dz[wz]
    loclz,nloclz = np.unique(dzz['TILELOCID'],return_counts=True)
    natloc = ~np.isin(dz['TILELOCID'],loclz)
    print('number of unique targets left around unassigned locations is '+str(np.sum(natloc)))
#    locs = np.copy(dz['TILELOCID'])
# 
# 
    print('reassigning TILELOCID for duplicates and finding rosette')
    #re-assigning "naked" targets; if we gave a targetid a tilelocid that was not assigned 
    #by the same target was available at a location that was assigned, we re-assign its tilelocid
    nch = 0
    nbl = 0
    tlids = dz['TILELOCIDS']
    ros = np.zeros(len(dz))
    rosr = np.zeros(len(dz))
    for ii in range(0,len(dz['TILEID'])): #not sure why, but this only works when using loop for Table.read but array option works for fitsio.read
        ti = dz[ii]['TILEID']
        rosn = tile2rosette(ti) #get rosette id
        rosr[ii] = calc_rosr(rosn,dz[ii]['RA'],dz[ii]['DEC']) #calculates distance in degrees from rosette center
        ros[ii] = rosn
     

    dz['rosette_number'] = ros
    dz['rosette_r'] = rosr
    print('rosette number and the number on each rosette')
    print(np.unique(dz['rosette_number'],return_counts=True))

    loco,fzo = common.comp_tileloc(dz)
    pd = dict(zip(loco,fzo))
    probl = np.zeros(len(dz))
    for i in range(0,len(dz)):
        probl[i] = pd[dz['TILELOCID'][i]]
    dz['FRACZ_TILELOCID'] = probl

    #write out FRACZ_TILELOCID info
    #loco = np.array(loco).astype(dz['TILELOCID'].dtype)
    co = Table()
    co['TILELOCID'] = loco
    co['FRACZ_TILELOCID'] = fzo
    cof = outf.strip('_full_noveto.dat.fits')+'_comp_tileloc.fits'
    print('writing comp_tileloc completeness to '+cof)
    co.write(cof,overwrite=True,format='fits')    


    print('sum of 1/FRACZ_TILELOCID, 1/COMP_TILE, length of input, number of inputs with obs; no longer rejecting unobserved loc, so wont match')
    print(np.sum(1./dz[wz]['FRACZ_TILELOCID']),np.sum(1./dz[wz]['COMP_TILE']),len(dz),len(dz[wz]))
    
    oct = np.copy(dz['COMP_TILE'])
    if bitweightfile is not None:
        fb = Table.read(bitweightfile)
        dz = join(dz,fb,keys=['TARGETID'])
    wz = dz['LOCATION_ASSIGNED'] == 1 #join re-ordered array, reget mask for assigned locations and check comp_tile
    print('length after join with bitweight file and sum of 1/comp_tile',len(dz),np.sum(1./dz[wz]['COMP_TILE']),len(dz[wz]))
    #print('check comp_tile array',np.array_equal(oct,dz['COMP_TILE']))

    

    comments = ["SV3 'full' LSS catalog for data without any vetos applied","entries are for all targetid that showed up in POTENTIAL_ASSIGNMENTS"]
    common.write_LSS(dz,outf,comments)


def tile2rosette(tile):
    if tile < 433:
        return (tile-1)//27
    else:
        if tile >= 433 and tile < 436:
            return 13
        if tile >= 436 and tile < 439:
            return 14
        if tile >= 439 and tile < 442:
            return 15
        if tile >= 442 and tile <=480:
            return (tile-442)//3
            
        if tile > 480:
            return tile//30    
    return 999999 #shouldn't be any more?

def calc_rosr(rosn,ra,dec):
    #given rosetter number and ra,dec, calculate distance from center 
    roscen = {0:(150.100,2.182),1:(179.6,0),2:(183.1,0),3:(189.9,61.8),4:(194.75,28.2)\
    ,5:(210.0,5.0),6:(215.5,52.5),7:(217.8,34.4),8:(216.3,-0.6),9:(219.8,-0.6)\
    ,10:(218.05,2.43),11:(242.75,54.98),12:(241.05,43.45),13:(245.88,43.45),14:(252.5,34.5)\
    ,15:(269.73,66.02),16:(194.75,24.7),17:(212.8,-0.6),18:(269.73,62.52),19:(236.1,43.45)}
    ra = ra*np.pi/180.
    dec = dec*np.pi/180.
    rac,decc = roscen[rosn]
    rac = rac*np.pi/180.
    decc = decc*np.pi/180.
    cd = np.sin(dec)*np.sin(decc)+np.cos(dec)*np.cos(decc)*np.cos(rac-ra)
    ad = np.arccos(cd)*180./np.pi
    if ad > 2.5:
        print(rosn,ra,dec,rac,decc)
    return ad

def apply_veto(fin, fout, ebits=None, zmask=False, maxp=3400):
    '''
    fl is a string with the path to the file name to load
    fout is a string with the path to the outpur file
    ebits are the new imaging mask bits to apply
    zmask is whether or not to apply any zmask
    maxp is the maximum priority to keep in the data files
    '''
    ff = Table.read(fin)
    print('length of input '+str(len(ff)))
    seld = ff['GOODHARDLOC'] == 1
    print('length after cutting to good locations '+str(len(ff[seld])))
    if '.dat' in fin:
        seld &= ff['PRIORITY_INIT'] <= maxp
        print('length after cutting locations with priority_init > '+str(maxp)+': '+str(len(ff[seld])))
    if '.ran' in fin:
        seld &= ff['ZPOSSLOC'] == 1
        print('length after cutting locations where target type could not be observed: '+str(len(ff[seld])))
        seld &= ff['PRIORITY'] <= maxp
        print('length after cutting locations with priority > '+str(maxp)+': '+str(len(ff[seld])))


    ff = ff[seld]

    ''' DO NOT APPLY TO MOCKS
    if ebits is not None:
        print('number before imaging mask '+str(len(ff)))
        if ebits == 'lrg_mask':
            sel = ff['lrg_mask'] == 0
            ff = ff[sel]
        else:
            ff = cutphotmask(ff,ebits)
        print('number after imaging mask '+str(len(ff)))
    '''

    if zmask:
        whz = ff['Z'] < 1.6
        ff = ff[whz]

        fzm = Table.read('/global/homes/m/mjwilson/desi/DX2DROPOUT/radial_mask.fits')
        zma = []
        for z in ff['Z']:
            zind = int(z/1e-6)
            zma.append(fzm[zind]['RADIAL_MASK'])
        zma = np.array(zma)
        wm = zma == 0
        ff = ff[wm]

    if '.dat' in fin:
#NOT NEEDED IN MOCKS        ff['Z'].name = 'Z_not4clus'
        print('updating completenes')
        compa = []
        tll = []
        ti = 0
        ff.sort('TILES')
        nts = len(np.unique(ff['TILES']))
        tlsl = ff['TILES']
        tlslu = np.unique(tlsl)
        laa = ff['LOCATION_ASSIGNED']
        print('TILELOCID_ASSIGNED',np.unique(ff['TILELOCID_ASSIGNED'],return_counts=True))

        #for tls in np.unique(dz['TILES']): #this is really slow now, need to figure out a better way
        i = 0
        while i < len(ff):
            tls  = []
            tlis = []
            nli = 0
            nai = 0

            while tlsl[i] == tlslu[ti]:
                nli += 1
                nai += laa[i]
                i += 1
                if i == len(ff):
                    break

            if ti%1000 == 0:
                print('at tiles '+str(ti)+' of '+str(nts))

            cp = nai/nli#no/nt
            #print(tls,cp,no,nt)
            compa.append(cp)
            tll.append(tlslu[ti])
            ti += 1
        comp_dicta = dict(zip(tll, compa))
        fcompa = []
        for tl in ff['TILES']:
            fcompa.append(comp_dicta[tl])
        ff['COMP_TILE'] = np.array(fcompa)
        wz = ff['ZWARN'] != 999999
        wz &= ff['ZWARN']*0 == 0
        wz &= ff['ZWARN'] != 1.e20
        print('sum of 1/FRACZ_TILELOCID, 1/COMP_TILE, and length of input; should approximately match')
        print(np.sum(1./ff[wz]['FRACZ_TILELOCID']),np.sum(1./ff[wz]['COMP_TILE']),len(ff))

    if '.ran' in fin:
        print('area is '+str(len(ff)/2500))
    comments = ["'full' LSS catalog without after vetos for priority, good hardware and imaging quality","entries are for targetid that showed up in POTENTIAL_ASSIGNMENTS"]
    common.write_LSS(ff, fout, comments)


def mkclusdat(fl, weightmd='tileloc', zmask=False, tp='', dchi2=9, tsnrcut=80, rcut=None, ntilecut=0, ccut=None, ebits=None, nreal=128):
    '''
    fl is the root of the input/output file
    weighttileloc determines whether to include 1/FRACZ_TILELOCID as a completeness weight
    zmask determines whether to apply a mask at some given redshift
    tp is the target type
    dchi2 is the threshold for keeping as a good redshift
    tnsrcut determines where to mask based on the tsnr2 value (defined below per tracer)
    '''
    ff = Table.read(fl+'full.dat.fits')

    wzm = ''
    if zmask:
        wzm = 'zmask_'
    if rcut is not None:
        wzm += 'rmin'+str(rcut[0])+'rmax'+str(rcut[1])+'_'
    if ntilecut > 0:
        wzm += 'ntileg'+str(ntilecut)+'_'
    if ccut is not None:
        wzm += ccut+'_' #you could change this to however you want the file names to turn out

    if ccut == 'main':
        if tp != 'LRG':
            print('this is only defined for LRGs!' )
        else:
            lrgmaintar = Table.read('/global/cfs/cdirs/desi/survey/catalogs/main/LSS/LRGtargetsDR9v1.1.1.fits',columns=['TARGETID'])
            sel = np.isin(ff['TARGETID'],lrgmaintar['TARGETID'])
            print('numbers before/after cut:')
            print(len(ff),len(ff[sel]))
            ff = ff[sel]
            ff.write(fl+wzm+'full.dat.fits',format='fits',overwrite='True')

##NO NEED FOR MOCKS    ff['Z_not4clus'].name = 'Z'
    '''
    This is where redshift failure weights go
    '''

    ff['WEIGHT_ZFAIL'] = np.ones(len(ff))

    '''
    #The LRGs just have this fairly ad hoc model that AJR fit in the notebook, definitely needs refinement/automation
    if tp == 'LRG':
        fibfluxz = ff['FIBERFLUX_Z']/ff['MW_TRANSMISSION_Z']
        coeff = [117.46,-60.91,11.49,-0.513] #from polyfit, 3rd to zeroth order in 1/fiberflu
        efs = coeff[-1]+coeff[-2]*(1/fibfluxz)+coeff[-3]*(1/fibfluxz)**2.+coeff[-4]*(1/fibfluxz)**3.
        ems = erf((ff['TSNR2_LRG']-13.2)/39.7)*.9855
        ff['WEIGHT_ZFAIL'] = 1./(1. -(1.-ems)*efs)
    '''
    '''
    One could plug in imaging systematic weights here
    Probably better to put it here so that full file only gets written out once and includes
    all of the weights
    '''


    outf = fl+wzm+'clustering.dat.fits'
    wz = ff['ZWARN'] == 0
    print('length before cutting to objects with redshifts '+str(len(ff)))
    print('length after cutting to zwarn == 0 '+str(len(ff[wz])))

    if tp == 'QSO':
        #good redshifts are currently just the ones that should have been defined in the QSO file when merged in full
        wz = ff['Z']*0 == 0
        wz &= ff['Z'] != 999999
        wz &= ff['Z'] != 1.e20
        wz &= ff['ZWARN'] != 999999
        wz &= ff['TSNR2_QSO'] > tsnrcut

    if tp == 'ELG' or tp == 'ELG_HIP':
##NO IN MOCKS        wz = ff['o2c'] > dchi2
        wz = ff['ZWARN']*0 == 0
        wz &= ff['ZWARN'] != 999999
        print('length after oII cut '+str(len(ff[wz])))
        wz &= ff['LOCATION_ASSIGNED'] == 1
        print('length after also making sure location assigned '+str(len(ff[wz])))
        wz &= ff['TSNR2_ELG'] > tsnrcut
        print('length after tsnrcut '+str(len(ff[wz])))

    if tp == 'LRG':
        print('applying extra cut for LRGs')
        # Custom DELTACHI2 vs z cut from Rongpu
        #drz = (10**(3 - 3.5*ff['Z']))
        #mask_bad = (drz>30) & (ff['DELTACHI2']<30)
        #mask_bad |= (drz<30) & (ff['DELTACHI2']<drz)
        #mask_bad |= (ff['DELTACHI2']<10)
        #wz &= ff['Z']<1.4
        #wz &= (~mask_bad)
        wz &= ff['ZWARN']*0 == 0
        wz &= ff['ZWARN'] != 999999
        wz &= ff['ZWARN'] != 1.e20

##THIS IS SUBSTITUTE BY NEXT        selg = ssr_tools.LRG_goodz(ff)
        wz &= ff['ZWARN']==0
        wz &= ff['Z']<1.5
        #do not apply  wz &= ff['DELTACHI2']>15

        #wz &= ff['DELTACHI2'] > dchi2
        print('length after Rongpu cut '+str(len(ff[wz])))
        wz &= ff['TSNR2_ELG'] > tsnrcut
        print('length after tsnrcut '+str(len(ff[wz])))

    if tp[:3] == 'BGS':
        print('applying extra cut for BGS')
        wz &= ff['DELTACHI2'] > dchi2
        print('length after dchi2 cut '+str(len(ff[wz])))
        wz &= ff['TSNR2_BGS'] > tsnrcut
        print('length after tsnrcut '+str(len(ff[wz])))



    ff = ff[wz]
    print('length after cutting to good z '+str(len(ff)))
    print('minimum,maximum Z',min(ff['Z']),max(ff['Z']))
    ff['WEIGHT'] = ff['WEIGHT_ZFAIL']
    ff['WEIGHT_COMP'] = np.ones(len(ff))
    if weightmd == 'tileloc':
        ff['WEIGHT_COMP'] = 1./ff['FRACZ_TILELOCID']

    if weightmd == 'probobs' :
        print('sources with zero value in PROB_OBS', len(ff[(ff['PROB_OBS']==0)]))
        ff = ff[(ff['PROB_OBS']!=0)]
        ff['WEIGHT_COMP'] *= 1./ff['PROB_OBS']


        #nassign = nreal*ff['PROB_OBS']+1 #assignment in actual observation counts
        #ff['WEIGHT_COMP'] *= (nreal+1)/nassign#1./ff['PROB_OBS']
        print(np.min(ff['WEIGHT']),np.max(ff['WEIGHT']))

    ff['WEIGHT'] *= ff['WEIGHT_COMP']
    if zmask:
        whz = ff['Z'] < 1.6
        ff = ff[whz]

        fzm = Table.read('/global/homes/m/mjwilson/desi/DX2DROPOUT/radial_mask.fits')
        zma = []
        for z in ff['Z']:
            zind = int(z/1e-6)
            zma.append(fzm[zind]['RADIAL_MASK'])
        zma = np.array(zma)
        wm = zma == 0
        ff = ff[wm]
    #apply any cut on rosette radius
    if rcut is not None:
        wr = ff['rosette_r'] > rcut[0]
        wr &= ff['rosette_r'] <  rcut[1]
        print('length before rosette radius cut '+str(len(ff)))
        ff = ff[wr]
        print('length after rosette radius cut '+str(len(ff)))
    #apply cut on ntile
    if ntilecut > 0:
        print('length before ntile cut '+str(len(ff)))
        wt = ff['NTILE'] > ntilecut
        ff = ff[wt]
        print('length after ntile cut '+str(len(ff)))
    if ccut == 'notQSO':
        wc = (ff['SV3_DESI_TARGET'] & sv3_targetmask.desi_mask['QSO']) ==  0
        print('length before cutting to not QSO '+str(len(ff)))
        ff = ff[wc]
        print('length after cutting to not QSO '+str(len(ff)))
    if ccut == 'zQSO':
        wc = ff['SPECTYPE'] ==  'QSO'
        print('length before cutting to spectype QSO '+str(len(ff)))
        ff = ff[wc]
        print('length after cutting to spectype QSO '+str(len(ff)))

    #select down to specific columns below and then also split N/S
    #if tp[:3] == 'BGS':
    #    kl = ['RA','DEC','Z','WEIGHT','TARGETID','NTILE','rosette_number','rosette_r','TILES','WEIGHT_ZFAIL','FRACZ_TILELOCID']
    #else:
    kl = ['RA','DEC','Z','WEIGHT','TARGETID','NTILE','COMP_TILE','rosette_number','rosette_r','TILES','WEIGHT_ZFAIL','FRACZ_TILELOCID','PROB_OBS','BITWEIGHTS']

    if tp[:3] == 'BGS':
        ff['flux_r_dered'] = ff['FLUX_R']/ff['MW_TRANSMISSION_R']
        kl.append('flux_r_dered')
        print(kl)

    ff.keep_columns(kl)#,'PROB_OBS'
    print('minimum,maximum weight')
    print(np.min(ff['WEIGHT']),np.max(ff['WEIGHT']))

    comments = ["SV3 'clustering' LSS catalog for data, all regions","entries are only for data with good redshifts"]
    common.write_LSS(ff, outf, comments)



def mkfullran(fs,indir,randir,rann,imbits,outf,tp,pd,bit,desitarg='SV3_DESI_TARGET',tsnr= 'TSNR2_ELG',notqso='',qsobit=4,fbcol='COADD_FIBERSTATUS',maxp=103400):
    '''
    indir is directory with inputs
    rann is the random file number (0-17)
    imbits are the maskbits for the imaging veto mask
    outf is the name (including full path) of the output file
    tp is the target type
    pd is the program, dark or bright
    bit is the bit to use to select to the target type
    randir doesn't get used anymore
    desitarg is the column to use to select the target type
    tsnr is the tsnr2 used for this sample
    '''
    
    #first, need to find locations to veto based on data
    #the same is done in mkfulldat
    #fs = fitsio.read(indir+'datcomb_'+pd+'_specwdup_Alltiles.fits')
    wf = fs[fbcol] == 0
    stlid = 10000*fs['TILEID'] +fs['LOCATION']
    gtl = np.unique(stlid[wf])
    #gtl now contains the list of good locations
    #we now want to load in the bigger data file with all the target info
    #we use it to find the locations where observations of the given type were not possible and then mask them
    zf = os.path.join(indir,'datcomb_'+pd+'_tarspecwdup_Alltiles.fits')
    dz = Table.read(zf) 
    wtype = ((dz[desitarg] & bit) > 0)
    if notqso == 'notqso':
        wtype &= ((dz[desitarg] & qsobit) == 0)

    
    dz = dz[wtype]#&wg]

    dz['RSDZ'].name = 'Z'
    dz['Z'][(dz['Z']==1.e20)] = 999999
    dz[tsnr][(dz[tsnr]==1.e20)] = 999999

    #print('length after selecting type and fiberstatus == 0 '+str(len(dz)))
    lznp = common.find_znotposs(dz)

    #lznp will later be used to veto
    #load in random file
    zf = os.path.join(indir,'rancomb_'+str(rann)+pd+'wdupspec_Alltiles.fits')
    dz = Table.read(zf)

    print(rann, 'before cuttin rancomb to desi target', len(dz))
    wk = ((dz[desitarg] & bit) > 0)
    if notqso == 'notqso':
        wk &= ((dz[desitarg] & qsobit) == 0)
    #####
    dz = dz[wk]

    print(rann, 'after cuttin rancomb to desi target', len(dz))
    wg = np.isin(dz['TILELOCID'],gtl)
    dz['GOODHARDLOC'] = np.zeros(len(dz)).astype('bool')
    dz['GOODHARDLOC'][wg] = 1

    wk = ~np.isin(dz['TILELOCID'],lznp)
    dz['ZPOSSLOC'] = np.zeros(len(dz)).astype('bool')
    dz['ZPOSSLOC'][wk] = 1

    #load in tileloc info for this random file and join it
    zfpd = os.path.join(indir,'rancomb_'+str(rann)+pd+'_Alltilelocinfo.fits')
    dzpd = Table.read(zfpd)
    dz = join(dz,dzpd,keys=['TARGETID'])
    print('length before cutting to good positions '+str(len(dz)))
    #cut to good and possible locations
    #wk = ~np.isin(dz['TILELOCID'],lznp)
    #wk &= np.isin(dz['TILELOCID'],gtl)
    #dz = dz[wk]    
    print('length after cutting to good positions '+str(len(dz)))
    '''THIS IS NOT NEEDED TO MOCK
    #get all the additional columns desired from original random files through join
    tarf = Table.read(os.path.join(randir, 'alltilesnofa.fits'))
    delcols = ['RA','DEC','SV3_DESI_TARGET','BGS_TARGET','MWS_TARGET','SUBPRIORITY','OBSCONDITIONS','NUMOBS_INIT',
        'NUMOBS_MORE','ZWARN','PRIORITY']
    tarf.remove_columns(delcols)
    print('columns left in alltilesnofa',tarf.columns)
    dz = join(dz,tarf,keys=['TARGETID'])
    '''

    dz.remove_columns(['PRIORITY'])
    data_specwdup = Table.read(os.path.join(indir,'datcomb_'+pd+'_specwdup_Alltiles.fits'))
    data_specwdup.keep_columns(['TILEID', 'PRIORITY', 'LOCATION','FIBER'])
    dz = join(dz,data_specwdup,keys=['TILEID', 'LOCATION', 'FIBER'],join_type='left')

    dz[tsnr][(dz[tsnr]==1.e20)] = 999999

    #apply imaging vetos
#NO NEED FOR MOCKS    dz = common.cutphotmask(dz,imbits)
    print('length after cutting to based on imaging veto mask '+str(len(dz)))
    #pl = np.copy(dz['PRIORITY']).astype(float)#dz['PRIORITY']
    #sp = pl <= 0
    #pl[sp] = .1
    dz['GOODPRI'] = np.zeros(len(dz)).astype('bool')
    sel = dz['PRIORITY'] <= maxp
    dz['GOODPRI'][sel] = 1
    

    dz['sort'] =  dz['GOODPRI']*dz['GOODHARDLOC']*dz['ZPOSSLOC']*(1+dz[tsnr])
    #dz[tsnr]*dz['GOODHARDLOC']*dz['ZPOSSLOC']+dz['GOODHARDLOC']*dz['ZPOSSLOC']+dz['GOODHARDLOC']*dz['ZPOSSLOC']/pl
    #sort by tsnr, like done for data, so that the highest tsnr are kept
    dz.sort('sort') 
    dz = unique(dz,keys=['TARGETID'],keep='last')
    print('length after cutting to unique TARGETID '+str(len(dz)))
    dz['rosette_number'] = 0
    dz['rosette_r'] = np.zeros(len(dz))
    for ii in range(0,len(dz)):
        rosn = tile2rosette(dz[ii]['TILEID'])
        rosd = calc_rosr(rosn,dz[ii]['RA'],dz[ii]['DEC']) #calculates distance in degrees from the rosette center
        dz[ii]['rosette_number'] = rosn
        dz[ii]['rosette_r'] = rosd
    print(np.unique(dz['NTILE']))
    if int(rann) < 10:
        cof = Table.read(outf[:-23]+'_comp_tile.fits')
    else:
        cof = Table.read(outf[:-24]+'_comp_tile.fits')  
    comp_dicta = dict(zip(cof['TILES'], cof['COMP_TILE']))
    fcompa = []
    tls = dz['TILES']
    ctls = cof['TILES']
    ctiles = np.zeros(len(dz))
    tlsd = np.isin(tls,cof['TILES'])
    print('number of tiles groups in randoms not in data '+str(len(np.unique(tls[~tlsd]))))
    for i in range(0,len(tls)):
        if tlsd[i]:#np.isin(tl,ctls):
            ctiles[i] = comp_dicta[tls[i]]
        #    fcompa.append(comp_dicta[tl]) 
        #else:
        #    fcompa.append(0)
    dz['COMP_TILE'] = ctiles#np.array(fcompa)
    wc0 = dz['COMP_TILE'] == 0
    print('number of randoms in 0 completeness regions '+str(len(dz[wc0])))   
    
    
    comments = ["SV3 'full' LSS catalog for random # "+str(rann)+" without any vetos applied","entries are for all targetid that showed up in POTENTIAL_ASSIGNMENTS"]
    common.write_LSS(dz,outf,comments)

def mkclusran(fl, rann, rcols=['Z','WEIGHT'], zmask=False, tsnrcut=80, tsnrcol='TSNR2_ELG', rcut=None, ntilecut=0, ccut=None, ebits=None):
    '''
    fl is the root of our catalog file names
    rann is the random number
    rcols are the columns that we randomly select from the data file
    zmask is whether or not we mask out certain redshift
    tsnrcut is the tsnr2 value below which we discard data
    tsnrcol is the specific column used for the tsnrcut
    '''
    
    wzm = ''
    if zmask:
        wzm = 'zmask_'
    if rcut is not None:
        wzm += 'rmin'+str(rcut[0])+'rmax'+str(rcut[1])+'_'
    if ntilecut > 0:
        wzm += 'ntileg'+str(ntilecut)+'_'    
    if ccut is not None:
        wzm += ccut+'_'  #you could change this to however you want the file names to turn out

    #load in data clustering catalog
    fcd = Table.read(fl+wzm+'clustering.dat.fits')
    #load in full random file
    ffr = Table.read(fl+str(rann)+'_full.ran.fits')

    #mask mask on tsnr
    wz = ffr[tsnrcol] > tsnrcut
    ffc = ffr[wz]
    print('length after,before tsnr cut:')
    print(len(ffc),len(ffr))
    #apply any cut on rosette radius
    if rcut is not None:
        wr = ffc['rosette_r'] > rcut[0]
        wr &= ffc['rosette_r'] <  rcut[1]
        print('length before rosette radius cut '+str(len(ffc)))
        ffc = ffc[wr]
        print('length after rosette radius cut '+str(len(ffc)))
    #apply cut on ntile
    if ntilecut > 0:
        print('length before ntile cut '+str(len(ffc)))
        wt = ffc['NTILE'] > ntilecut
        ffc = ffc[wt]
        print('length after ntile cut '+str(len(ffc)))    


    #randomly sample data rows to apply redshifts, weights, etc. to randoms
    inds = np.random.choice(len(fcd),len(ffc))
    dshuf = fcd[inds]
    kl =  ['RA','DEC','TARGETID','NTILE','COMP_TILE','rosette_number','rosette_r','TILES'] + rcols 

    for col in rcols: 
        ffc[col] = dshuf[col] 
    #cut to desired small set of columns and write out files, splitting N/S as well
#    wn = ffc['PHOTSYS'] == 'N'
    
    ffc.keep_columns(kl)  
    outf =  fl+wzm+str(rann)+'_clustering.ran.fits' 

    comments = ["SV3 'clustering' LSS catalog for random #"+str(rann)+", all regions","columns that are not ra,dec are sampled from data with good redshifts"]
    common.write_LSS(ffc,outf,comments)


    #ffc.write(outf,format='fits', overwrite=True)
    '''NOT NEEDED FOR MOCKS
    outfn =  fl+wzm+'N_'+str(rann)+'_clustering.ran.fits' 
    fcdn = Table.read(fl+wzm+'N_clustering.dat.fits')
    ffcn = ffc[wn]
    inds = np.random.choice(len(fcdn),len(ffcn))
    dshuf = fcdn[inds]
    for col in rcols: 
        ffcn[col] = dshuf[col]     

    comments = ["SV3 'clustering' LSS catalog for random #"+str(rann)+", BASS/MzLS region","columns that are not ra,dec are sampled from data with good redshifts"]
    common.write_LSS(ffcn,outfn,comments)


    outfs =  fl+wzm+'S_'+str(rann)+'_clustering.ran.fits' 
    fcds = Table.read(fl+wzm+'S_clustering.dat.fits')
    ffcs = ffc[~wn]
    inds = np.random.choice(len(fcds),len(ffcs))
    dshuf = fcds[inds]
    for col in rcols: 
        ffcs[col] = dshuf[col]     

    comments = ["SV3 'clustering' LSS catalog for random #"+str(rann)+", DECaLS region","columns that are not ra,dec are sampled from data with good redshifts"]
    common.write_LSS(ffcs,outfs,comments)
    '''
