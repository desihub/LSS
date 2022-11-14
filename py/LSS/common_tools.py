import numpy as np
import fitsio
from astropy.table import Table,join,vstack
import datetime
import os
import sys

from desitarget.targetmask import obsmask, obsconditions, zwarn_mask

from LSS.tabulated_cosmo import TabulatedDESI
cosmo = TabulatedDESI()
dis_dc = cosmo.comoving_radial_distance

def dl(z):   # Luminosity distance from now to z
    return dis_dc(z)*(1.+z)

def dm(z):
    return 5.*np.log10(dl(z)) + 25.


#functions that shouldn't have any dependence on survey go here

def cut_specdat(dz,badfib=None):
    selz = dz['ZWARN'] != 999999
    selz &= dz['ZWARN']*0 == 0 #just in case of nans
    fs = dz[selz]

    #first, need to find locations to veto based data
    nodata = fs["ZWARN_MTL"] & zwarn_mask["NODATA"] != 0
    num_nod = np.sum(nodata)
    print('number with no data '+str(num_nod))
    badqa = fs["ZWARN_MTL"] & zwarn_mask.mask("BAD_SPECQA|BAD_PETALQA") != 0
    num_badqa = np.sum(badqa)
    print('number with bad qa '+str(num_badqa))
    nomtl = nodata | badqa
    wfqa = ~nomtl
    #veto fibers later determined to have poor success rates
    if badfib is not None:
        bad = np.isin(fs['FIBER'],badfib)
        print('number at bad fibers '+str(sum(bad)))
        wfqa &= ~bad
    return fs[wfqa]


def cutphotmask(aa,bits):
    print(str(len(aa)) +' before imaging veto' )
    keep = (aa['NOBS_G']>0) & (aa['NOBS_R']>0) & (aa['NOBS_Z']>0)
    for biti in bits:
        keep &= ((aa['MASKBITS'] & 2**biti)==0)
    aa = aa[keep]
    print(str(len(aa)) +' after imaging veto' )
    return aa

def find_znotposs_tloc(dz,priority_thresh=10000):
    #dz should contain the potential targets of a given type, after cutting bad fibers
    #priority_thresh cuts repeat observations (so, e.g., 3000 work for dark time main survey)
    tileids = np.unique(dz['TILEID'])
    ual = []
    ufl = []
    for tile in tileids:
        sel = dz['TILEID'] == tile
        dzs = dz[sel]
        sela = dzs['ZWARN'] != 999999
        sela &= dzs['ZWARN']*0 == 0
        tlida = np.unique(dzs[sela]['TILELOCID']) #tilelocids with good assignments
        #print(tile,len(sela),len(dzs),np.sum(sela))
        tida = np.unique(dzs[sela]['TARGETID']) #targetids with good assignments
        ua = ~np.isin(dzs['TARGETID'],tida) #columns that were not assigned
        #ua &= dzs['NUMOBS'] == 0
        ua &= dzs['PRIORITY'] > priority_thresh #columns corresponding to targets with priority indicating observations unfinished
        ua &= ~np.isin(dzs['TILELOCID'],tlida) #columns corresponding to tilelocid that were not assigned
        #combination then gives the tilelocid of unassigned targets that have not finished observation; these must have been blocked by something
        uatlids = np.unique(dzs[ua]['TILELOCID'])
        ual.append(uatlids)
        selp = dzs['PRIORITY'] > priority_thresh
        tlids_gp = np.unique(dzs[selp]['TILELOCID'])
        tlids_all = np.unique(dzs['TILELOCID'])
        tlids_full = tlids_all[~np.isin(tlids_all,tlids_gp)]
        print('done with tile '+str(tile),str(len(uatlids)),str(len(tlids_full)))
        ufl.append(tlids_full)
    print('concatenating')
    ualt = np.concatenate(ual)
    uflt = np.concatenate(ufl)
    return ualt,uflt


def find_znotposs(dz):

    dz.sort('TARGETID')
    tidnoz = []
    tids = np.unique(dz['TARGETID'])
    ti = 0
    i = 0

    print('finding targetids that were not observed')
    while i < len(dz):
        za = 0

        while dz[i]['TARGETID'] == tids[ti]:
            if dz[i]['ZWARN'] != 999999:
                za = 1
                #break
            i += 1
            if i == len(dz):
                break
        if za == 0:
            tidnoz.append(tids[ti])

        if ti%30000 == 0:
            print(ti)
        ti += 1


    selnoz = np.isin(dz['TARGETID'],tidnoz)
    tidsb = np.unique(dz[selnoz]['TILELOCID'])
    #dz = dz[selnoz]
    dz.sort('TILELOCID')
    tids = np.unique(dz['TILELOCID'])
    print('number of targetids with no obs '+str(len(tidnoz)))
    tlidnoz = []
    lznposs = []

    ti = 0
    i = 0

    while i < len(dz):
        za = 0

        while dz[i]['TILELOCID'] == tids[ti]:
            if dz[i]['ZWARN'] != 999999:
                za = 1
                #break
            i += 1
            if i == len(dz):
                break
        if za == 0:
            tlidnoz.append(tids[ti])
            #if np.isin(tids[ti],tidsb):
            #    lznposs.append(tids[ti])

        if ti%30000 == 0:
            print(ti,len(tids))
        ti += 1
    #the ones to veto are now the join of the two
    wtbtlid = np.isin(tlidnoz,tidsb)
    tlidnoz = np.array(tlidnoz)
    lznposs = tlidnoz[wtbtlid]
    print('number of locations where assignment was not possible because of priorities '+str(len(lznposs)))
    return lznposs

def comp_tile(dz):
    compa = []
    tll = []
    ti = 0
    print('getting completenes')
    #sorting by tiles makes things quicker with while statements below
    dz.sort('TILES')
    nts = len(np.unique(dz['TILES']))
    tlsl = dz['TILES']
    tlslu = np.unique(tlsl)
    laa = dz['LOCATION_ASSIGNED']

    i = 0
    while i < len(dz):
        tls  = []
        tlis = []
        nli = 0
        nai = 0

        while tlsl[i] == tlslu[ti]:
            nli += 1 #counting unique targetids within the given TILES value
            nai += laa[i] #counting the number assigned
            i += 1
            if i == len(dz):
                break

        if ti%1000 == 0:
            print('at tiles '+str(ti)+' of '+str(nts))
        cp = nai/nli #completeness is number assigned over number total
        compa.append(cp)
        tll.append(tlslu[ti])
        ti += 1
    return tll,compa

def comp_tileloc(dz):

    locl,nlocl = np.unique(dz['TILELOCID'],return_counts=True)
    wz = dz['LOCATION_ASSIGNED'] == 1
    dzz = dz[wz]

    loclz,nloclz = np.unique(dzz['TILELOCID'],return_counts=True)

    print('getting fraction assigned for each tilelocid')
    #should be one (sometimes zero, though) assigned target at each tilelocid and we are now counting how many targets there are per tilelocid
    #probability of assignment is then estimated as 1/n_tilelocid
    nm = 0
    nmt =0

    loco = []
    fzo = []
    nlist = 0
    nlistg1 = 0
    for i in range(0,len(locl)):
        if i%10000 == 0:
            print('at row '+str(i))
        nt = nlocl[i]
        loc = locl[i]
        w = loclz == loc
        nz = 0
        if len(loclz[w]) == 1:
            nz = nloclz[w][0] #these are supposed all be 1...
        else:
            nm += 1.
            nmt += nt
        if len(loclz[w]) > 1:
            print('why is len(loclz[w]) > 1?') #this should never happen
        loco.append(loc)
        frac = nz/nt
        #if type(frac) != float:
        #    if len(frac) > 1:
        #        nlistg1 += 1
        #    frac = frac[0]
        #    nlist += 1

        fzo.append(frac)
    print(str(nlist)+ ' were type list for some reason; '+str(nlistg1)+ ' had length greater than 1')
    print('number of fibers with no observation, number targets on those fibers')
    print(nm,nmt)


    return loco,fzo


def mknz(fcd,fcr,fout,bs=0.01,zmin=0.01,zmax=1.6,randens=2500.):
    '''
    fcd is the full path to the catalog file in fits format with the data; requires columns Z and WEIGHT
    fcr is the full path to the random catalog meant to occupy the same area as the data; assumed to come from the imaging randoms that have a density of 2500/deg2
    fout is the full path to the file name
    bs is the bin width for the n(z) calculation
    zmin is the lower edge of the first bin
    zmax is the upper edge of the last bin
    '''
    #cd = distance(om,1-om)
    ranf = fitsio.read_header(fcr,ext=1) #should have originally had 2500/deg2 density, so can convert to area
    area = ranf['NAXIS2']/randens
    print('area is '+str(area))

    df = fitsio.read(fcd)

    nbin = int((zmax-zmin)/bs)
    zhist = np.histogram(df['Z'],bins=nbin,range=(zmin,zmax),weights=df['WEIGHT'])
    outf = open(fout,'w')
    outf.write('#area is '+str(area)+'square degrees\n')
    outf.write('#zmid zlow zhigh n(z) Nbin Vol_bin\n')
    for i in range(0,nbin):
        zl = zhist[1][i]
        zh = zhist[1][i+1]
        zm = (zh+zl)/2.
        voli = area/(360.*360./np.pi)*4.*np.pi/3.*(dis_dc(zh)**3.-dis_dc(zl)**3.)
        nbarz =  zhist[0][i]/voli
        outf.write(str(zm)+' '+str(zl)+' '+str(zh)+' '+str(nbarz)+' '+str(zhist[0][i])+' '+str(voli)+'\n')
    outf.close()

def addnbar(fb,nran=18,bs=0.01,zmin=0.01,zmax=1.6,P0=10000,add_data=True,ran_sw='',ranmin=0):
    '''
    fb is the root of the file name, including the path
    nran is the number of random files to add the nz to
    bs is the bin size of the nz file (read this from file in future)
    zmin is the lower edge of the minimum bin (read this from file in future)
    zmax is the upper edge of the maximum bin (read this from file in the future)
    '''

    nzd = np.loadtxt(fb.replace(ran_sw,'')+'_nz.txt').transpose()[3] #column with nbar values
    fn = fb.replace(ran_sw,'')+'_clustering.dat.fits'
    #ff = fitsio.FITS(fn,'rw')
    #fd = Table(ff['LSS'].read())
    #fd = fitsio.read(fn) #reading in data with fitsio because it is much faster to loop through than table
    fd = Table(fitsio.read(fn))
    zl = fd['Z']
    nl = np.zeros(len(zl))
    for ii in range(0,len(zl)):
        z = zl[ii]
        zind = int((z-zmin)/bs)
        if z > zmin and z < zmax:
            nl[ii] = nzd[zind]
    mean_comp = len(fd)/np.sum(fd['WEIGHT'])
    print('mean completeness '+str(mean_comp))
    #del fd
    #ft = Table.read(fn)
    #ft['NZ'] = nl
    fd['NZ'] = nl
    #ff['LSS'].insert_column('NZ',nl)
    print(np.min(nl),np.max(nl))

    fkpl = 1./(1+nl*P0*mean_comp)
    #ft['WEIGHT_FKP'] = 1./(1+ft['NZ']*P0)
    if add_data:
        fd['WEIGHT_FKP'] = fkpl
        write_LSS(fd,fn)
    #fd = np.array(fd)
    #ff['LSS'].insert_column('WEIGHT_FKP',fkpl)
    #ff['LSS'].write(fd)
    #ff['LSS'].write_history("added NZ and WEIGHT_FKP columns on "+datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
    #ff.close()
    #ft.write(fn,format='fits',overwrite=True)
    print('done with data')
    for rann in range(ranmin,nran):
        fn = fb+'_'+str(rann)+'_clustering.ran.fits'
        #ff = fitsio.FITS(fn,'rw')
        #fd = ff['LSS'].read()
        fd = Table(fitsio.read(fn))
        #fd = fitsio.read(fn) #reading in data with fitsio because it is much faster to loop through than table
        zl = fd['Z']
        nl = np.zeros(len(zl))
        for ii in range(0,len(zl)):
            z = zl[ii]
            zind = int((z-zmin)/bs)
            if z > zmin and z < zmax:
                nl[ii] = nzd[zind]
        #del fd
        #ft = Table.read(fn)
        #ft['NZ'] = nl
        #ff['LSS'].insert_column('NZ',nl)
        fd['NZ'] = nl
        fkpl = 1./(1+nl*P0*mean_comp)
        fd['WEIGHT_FKP'] = fkpl
        write_LSS(fd,fn)
        #ff['LSS'].insert_column('WEIGHT_FKP',fkpl)
        #fd = np.array(fd)
        #ff['LSS'].write(fd)
        #ff['LSS'].write_history("added NZ and WEIGHT_FKP columns on "+datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        #ff.close()
        #ft['WEIGHT_FKP'] = 1./(1+ft['NZ']*P0)
        #ft.write(fn,format='fits',overwrite=True)
        print('done with random number '+str(rann))
    return True

def add_dered_flux(data,fcols=['G','R','Z','W1','W2']):
    #data should be table with fcols flux columns existing
    for col in fcols:
        data['flux_'+col.lower()+'_dered'] = data['FLUX_'+col]/data['MW_TRANSMISSION_'+col]
    return data

def add_ke(dat,zcol='Z'):#,n_processes=100):
    from multiprocessing import Pool
    #dat should be table with flux_g_dered and flux_r_dered
    #from kcorr package https://github.com/SgmAstro/DESI, needs to be added to path
    #
    #ke_code_root = '/global/homes/a/ajross/desicode/DESI_ke'
    ke_code_root = os.environ['LSSCODE']+'/LSS/py/LSS/DESI_ke'
    sys.path.append(ke_code_root)
    os.environ['CODE_ROOT'] = ke_code_root
    print(os.environ['CODE_ROOT'])
    from   smith_kcorr     import GAMA_KCorrection
    from   rest_gmr        import smith_rest_gmr,rest_gmr
    from   tmr_ecorr       import tmr_ecorr, tmr_q

    kcorr_r   = GAMA_KCorrection(band='R')
    kcorr_g   = GAMA_KCorrection(band='G')

    r_dered = 22.5 - 2.5*np.log10(dat['flux_r_dered'])
    g_dered = 22.5 - 2.5*np.log10(dat['flux_g_dered'])
    gmr = g_dered-r_dered

#     chunk_size = len(dat)//n_processes
#     list = []
#     for i in range(0,n_processes):
#         list.append(0)
#     def _wrapper(N):
#         mini = N*chunk_size
#         maxi = mini+chunk_size
#         if maxi > len(dat):
#             maxi = len(dat)
#         idx = np.arange(mini,maxi)
#         data = Table()
#         data['idx'] = idx
#         data['REST_GMR_0P1'], rest_gmr_0p1_warn = smith_rest_gmr(dat[zcol][mini:maxi], gmr[mini:maxi])
#         list[N] = data
#         #return data
#
#     with Pool(processes=n_processes+1) as pool:
#         #res = pool.map(_wrapper, np.arange(n_processes))
#         pool.map(_wrapper, np.arange(n_processes))
#
#     res = vstack(list)#vstack(res)
#     res.sort('idx')
#     res.remove_column('idx')
#     print(len(res),len(dat))
    selz = dat[zcol] > 0
    selz &= dat[zcol] < 4.5
    cols = ['REST_GMR_0P1','KCORR_R0P1','KCORR_G0P1','KCORR_R0P0','KCORR_G0P0','REST_GMR_0P0','EQ_ALL_0P0','EQ_ALL_0P1','ABSMAG_RP1','ABSMAG_RP0']
    for col in cols:
        dat[col] = np.zeros(len(dat))
    dat['REST_GMR_0P1'][selz], rest_gmr_0p1_warn = smith_rest_gmr(dat[zcol][selz], gmr[selz])
    dat['KCORR_R0P1'][selz] = kcorr_r.k(dat[zcol][selz], dat['REST_GMR_0P1'][selz])
    dat['KCORR_G0P1'][selz] = kcorr_g.k(dat[zcol][selz], dat['REST_GMR_0P1'][selz])
    dat['KCORR_R0P0'][selz] = kcorr_r.k_nonnative_zref(0.0, dat[zcol][selz], dat['REST_GMR_0P1'][selz])
    dat['KCORR_G0P0'][selz] = kcorr_g.k_nonnative_zref(0.0, dat[zcol][selz], dat['REST_GMR_0P1'][selz])
    dat['REST_GMR_0P0'][selz] = gmr[selz] - (dat['KCORR_G0P0'][selz] - dat['KCORR_R0P0'][selz])
    dat['EQ_ALL_0P0'][selz]   = tmr_ecorr(dat[zcol][selz], dat['REST_GMR_0P0'][selz], aall=True)
    dat['EQ_ALL_0P1'][selz]   = tmr_ecorr(dat[zcol][selz], dat['REST_GMR_0P1'][selz], aall=True)
    dat['ABSMAG_RP1'][selz] = r_dered[selz] -dm(dat[zcol][selz])-dat['KCORR_R0P1'][selz]-dat['EQ_ALL_0P1'][selz]
    dat['ABSMAG_RP0'][selz] = r_dered[selz] -dm(dat[zcol][selz])-dat['KCORR_R0P0'][selz]-dat['EQ_ALL_0P0'][selz]
    return dat
    #abg = g_dered -dm(data['Z'])


def join_etar(fn,tracer,tarver='1.1.1'):
    tr = tracer
    if tr == 'BGS_BRIGHT':
        tr = 'BGS_ANY'
    etar_fn = '/global/cfs/cdirs/desi/survey/catalogs/main/LSS/'+tr+'targets_pixelDR9v'+tarver+'.fits'
    ef = Table(fitsio.read(etar_fn))
    ef.remove_columns(['BRICKNAME','RA','DEC','PHOTSYS','DESI_TARGET'])
    df = fitsio.read(fn)
    print(len(df))
    df = join(df,ef,keys=['TARGETID'])
    print(len(df),'should match above')
    comments = ['Adding imaging mask column']
    write_LSS(df,fn,comments)


def add_veto_col(fn,ran=False,tracer_mask='lrg',rann=0,tarver='targetsDR9v1.1.1',redo=False):
    mask_fn = '/global/cfs/cdirs/desi/survey/catalogs/main/LSS/'+tracer_mask.upper()+tarver+'_'+tracer_mask+'imask.fits'
    if ran:
        mask_fn = '/global/cfs/cdirs/desi/survey/catalogs/main/LSS/randoms-1-'+str(rann)+tracer_mask+'imask.fits'
    maskf = fitsio.read(mask_fn)
    df = fitsio.read(fn)
    if np.isin(tracer_mask+'_mask',list(df.dtype.names)):
        print('mask column already in '+fn)
        if redo:
            df = Table(df)
            df.remove_columns([tracer_mask+'_mask'])
            print('will replace '+tracer_mask)
        else:
            return True
    else:
        print('adding '+tracer_mask)
    print(len(df))
    df = join(df,maskf,keys=['TARGETID'])
    print(len(df),'should match above')
    comments = ['Adding imaging mask column']
    write_LSS(df,fn,comments)

def apply_veto(fin,fout,ebits=None,zmask=False,maxp=3400):
    '''
    fl is a string with the path to the file name to load
    fout is a string with the path to the outpur file
    ebits are the new imaging mask bits to apply
    zmask is whether or not to apply any zmask
    maxp is the maximum priority to keep in the data files
    '''
    ff = Table(fitsio.read(fin))#+'full_noveto.'+dr+'.fits')
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

    if ebits is not None:
        print('number before imaging mask '+str(len(ff)))
        if ebits == 'lrg_mask':
            sel = ff['lrg_mask'] == 0
            ff = ff[sel]
        else:
            ff = cutphotmask(ff,ebits)
        print('number after imaging mask '+str(len(ff)))

    if zmask:
        whz = ff['Z'] < 1.6
        ff = ff[whz]

        fzm = fitsio.read('/global/homes/m/mjwilson/desi/DX2DROPOUT/radial_mask.fits')
        zma = []
        for z in ff['Z']:
            zind = int(z/1e-6)
            zma.append(fzm[zind]['RADIAL_MASK'])
        zma = np.array(zma)
        wm = zma == 0
        ff = ff[wm]

    if '.dat' in fin:
        ff['Z'].name = 'Z_not4clus'
        print('updating completeness')
        compa = []
        tll = []
        ti = 0
        ff.sort('TILES')
        nts = len(np.unique(ff['TILES']))
        tlsl = ff['TILES']
        tlslu = np.unique(tlsl)
        laa = ff['LOCATION_ASSIGNED']
        print('TILELOCID_ASSIGNED',np.unique(ff['TILELOCID_ASSIGNED'],return_counts=True))

        # for tls in np.unique(dz['TILES']): #this is really slow now, need to figure out a better way
        i = 0
        while i < len(ff):
            tls = []
            tlis = []
            nli = 0
            nai = 0

            while tlsl[i] == tlslu[ti]:
                nli += 1
                nai += laa[i]
                i += 1
                if i == len(ff):
                    break

            if ti % 1000 == 0:
                print('at tiles ' + str(ti) + ' of ' + str(nts))

            cp = nai / nli #no/nt
            # print(tls,cp,no,nt)
            compa.append(cp)
            tll.append(tlslu[ti])
            ti += 1
        comp_dicta = dict(zip(tll, compa))
        fcompa = []
        for tl in ff['TILES']:
            fcompa.append(comp_dicta[tl])
        ff['COMP_TILE'] = np.array(fcompa)
        wz = ff['ZWARN'] != 999999
        wz &= ff['ZWARN'] * 0 == 0
        wz &= ff['ZWARN'] != 1.e20
        print('sum of 1/FRACZ_TILELOCID, 1/COMP_TILE, and length of input; should approximately match')
        print(np.sum(1. / ff[wz]['FRACZ_TILELOCID']), np.sum(1. / ff[wz]['COMP_TILE']), len(ff))

    if '.ran' in fin:
        print('area is ' + str(len(ff) / 2500))
    comments = ["'full' LSS catalog without after vetos for priority, good hardware and imaging quality","entries are for targetid that showed up in POTENTIAL_ASSIGNMENTS"]
    write_LSS(ff, fout, comments)

#     tmpfn = fout+'.tmp'
#     if os.path.isfile(tmpfn):
#         os.system('rm '+tmpfn)
#     fd = fitsio.FITS(tmpfn, "rw")
#     fd.write(np.array(ff),extname='LSS')
#     fd['LSS'].write_history("created (or over-written) on "+datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
#     fd.close()
#     os.system('mv '+tmpfn+' '+fout)
    #ff.write(fout,overwrite=True,format='fits')


def write_LSS(ff, outf, comments=None):
    '''
    ff is the structured array/Table to be written out as an LSS catalog
    outf is the full path to write out
    comments is a list of comments to include in the header
    '''
    tmpfn = outf + '.tmp'
    if os.path.isfile(tmpfn):
        os.system('rm ' + tmpfn)
    fd = fitsio.FITS(tmpfn, "rw")
    fd.write(np.array(ff), extname='LSS')
    if comments is not None:
        for comment in comments:
            fd['LSS'].write_comment(comment)
    fd['LSS'].write_history("updated on " + datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
    fd.close()
    print('closed fits file')
    os.system('mv ' + tmpfn + ' ' + outf)
    print('moved output to ' + outf)


def create_sky(dirname, columns=None, format_output='fits', release='1.1.1', version='main', program='dark', dr='dr9', nfiles=10, mpicomm=None):
    """ It is a real nightmare to work with all the fits file in hpdirname_list more than 800 for each directory. Rewrite it in a smaller number of files. """
    from glob import glob
    import mpytools as mpy
    from mpi4py import MPI
    from fiberassign.fba_launch_io import get_desitarget_paths

    if mpicomm is None:
        print('This code works under MPI paradigm throught cosmodesi/mpytools Software. Please pass a mpi communicator as input: mpicomm=MPI.COMM_WORLD')
        sys.exit(2)

    if columns is None:
        columns = ["RA", "DEC", "TARGETID", "DESI_TARGET", "BGS_TARGET", "MWS_TARGET", "SUBPRIORITY", "OBSCONDITIONS", "PRIORITY_INIT", "NUMOBS_INIT"]

    # Now load the sky target files.  These are main-survey files that we will
    # force to be treated as the survey type of the other target files.
    mydirs = get_desitarget_paths(release, version, program, dr=dr)
    skydirs = [mydirs["sky"]]
    if os.path.isdir(mydirs["skysupp"]):
        skydirs.append(mydirs["skysupp"])

    for hpdirname in skydirs:
        fns = []
        if mpicomm.rank == 0:
            fns = glob(os.path.join(hpdirname, "*fits"))
        fns = list(mpy.bcast(fns, mpiroot=0))

        # Create mlpytools.Catalog
        targets = mpy.Catalog.read(fns, mpicomm=mpicomm)

        if format_output == 'fits':
            # first save it in nfiles fits files.
            # Create filenames to write targets in less files.
            basename = '-'.join(os.path.basename(fns[0]).split('-')[:-2])
            fns_to_save = [os.path.join(dirname, f'{basename}-{i}.fits') for i in range(nfiles)]

            start = MPI.Wtime()
            targets[columns].write(fns_to_save, filetype='fits')
            mpicomm.Barrier()
            if mpicomm.rank == 0: print(f'Write {basename} in {nfiles} done in {MPI.Wtime() - start:2.2f} s.')
        elif format_output == 'bigfile':
            # Save it on one bigfile (more efficient ?)
            start = MPI.Wtime()
            targets[columns].write(os.path.join(dirname, f'{basename}'), filetype='bigfile')
            mpicomm.Barrier()
            if mpicomm.rank == 0: print(f'Write {basename} in {nfiles} done in {MPI.Wtime() - start:2.2f} s.')
        else:
            print('Only fits or bigfile format are extected...')
            sys.exit(2)


def combtiles_pa_wdup(tiles, fbadir, outdir, tarf, addcols=['TARGETID', 'RA', 'DEC'], fba=True, tp='dark', ran='ran'):
    if ran == 'dat':
        # addcols.append('PRIORITY')
        addcols.append('PRIORITY_INIT')
        addcols.append('DESI_TARGET')
    s = 0
    td = 0
    # tiles.sort('ZDATE')
    print(len(tiles))
    outf = outdir + '/' + ran + 'comb_' + tp + 'wdup.fits'
    if fba:
        pa_hdu = 'FAVAIL'
    tl = []
    for tile in tiles['TILEID']:
        if fba:
            ffa = fbadir + '/fba-' + str(tile).zfill(6) + '.fits'
        if os.path.isfile(ffa):
            fa = Table(fitsio.read(ffa, ext=pa_hdu))
            sel = fa['TARGETID'] >= 0
            fa = fa[sel]
            td += 1
            fa['TILEID'] = int(tile)
            tl.append(fa)
            print(td, len(tiles))
        else:
            print('did not find ' + ffa)
    dat_comb = vstack(tl)
    print(len(dat_comb))
    tar_in = fitsio.read(tarf, columns=addcols)
    dat_comb = join(dat_comb, tar_in, keys=['TARGETID'])
    print(len(dat_comb))

    dat_comb.write(outf, format='fits', overwrite=True)
    print('wrote ' + outf)
    return dat_comb

def combtiles_assign_wdup(tiles,fbadir,outdir,tarf,addcols=['TARGETID','RSDZ','TRUEZ','ZWARN','PRIORITY'],fba=True,tp='dark'):

    s = 0
    td = 0
    #tiles.sort('ZDATE')
    print(len(tiles))
    outf = outdir+'/datcomb_'+tp+'assignwdup.fits'
    if fba:
        pa_hdu = 'FASSIGN'
    tl = []
    for tile in tiles['TILEID']:
        if fba:
            ffa = fbadir+'/fba-'+str(tile).zfill(6)+'.fits'
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

    dat_comb.write(outf,format='fits', overwrite=True)
    print('wrote '+outf)
    return dat_comb

def addNS(tab):
    '''
    given a table that already includes RA,DEC, add PHOTSYS column denoting whether
    the data is in the DECaLS ('S') or BASS/MzLS ('N') photometric region
    '''
    wra = (tab['RA'] > 100-tab['DEC'])
    wra &= (tab['RA'] < 280 +tab['DEC'])
    tab['PHOTSYS'] = 'S'
    seln = dec > 32.375
    seln &= wra
    tab['PHOTSYS'][seln] = 'N'
    return tab
