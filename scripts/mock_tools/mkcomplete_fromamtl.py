def get_parent_clus_fromfull(tracer,realization,reg,base_dir='/global/cfs/cdirs/desi/mocks/cai/LSS/DA2/mocks/',
                             mockver='GLAM-Uchuu_v1'):
    tardir = base_dir+mockver
    tar_data_fn = tardir+'/forFA'+str(realization)+'.fits'
    tar_data = read_file(tar_data_fn,columns=['TARGETID','RSDZ'])
    tar_data.rename_column('RSDZ','Z')
    LSSdir = tardir+'/altmtl'+str(realization)+'/loa-v1/mock'+str(realization)+'/LSScats/'
    pars = params(tracer)
    zmin = pars.zmin
    zmax = pars.zmax
    if tracer == 'QSO':
        P0 = 6000
        zmin = 0.8
        zmax = 3.5
    if tracer[:3] == 'LRG':
        P0 = 10000
        zmin = 0.4
        zmax = 1.1
    if tracer[:3] == 'ELG':
        P0 = 4000
        zmin = 0.8
        zmax = 1.6
    in_data_fn = LSSdir+tracer+'_full_HPmapcut.dat.h5'    
    print(in_data_fn)
    in_data = read_file(in_data_fn,columns=['TARGETID','RA','DEC','NTILE','ZWARN'])
    
    sel_con = in_data['TARGETID'] < 419430400000000 #remove contaminants
    in_data = in_data[sel_con]
    ld = len(in_data)
    tids, in_ind, orig_ind = np.intersect1d(in_data['TARGETID'], tar_data['TARGETID'], return_indices=True)
    in_data = in_data[in_ind]
    print('input length was '+str(ld)+' and new length is '+str(len(in_data))+', should not have changed')
    tar_data = tar_data[orig_ind]
    in_data['Z'] = tar_data['Z']
    if 'ELG' in tracer:
        rans = np.random.random(len(in_data))
        subzhz = rans < 0.76
        subzlz = rans < 0.96
        zsplit = in_data['Z'] < 1.49
        subz = subzhz 
        subz |= (subzlz & zsplit)
        
        print('for zfail, keeping '+str(np.sum(subz)/len(in_data)))
        in_data =in_data[subz]
              
    cntl = np.zeros(8)
    sela = in_data['ZWARN'] != 999999
    in_data['WEIGHT'] = np.ones(len(in_data))
    for nt in range(0,8):
        seld = in_data['NTILE'] == nt
        comp = np.sum(sela&seld)/np.sum(seld)		
        cntl[nt] = comp
    in_data['WEIGHT'] = cntl[in_data['NTILE']]
        
    print('completeness as a function of NTILE is '+str(cntl))
    sel_ngc = common.splitGC(in_data)
    if reg == 'NGC':
        in_data = in_data[sel_ngc]
    if reg == 'SGC':
        in_data = in_data[~sel_ngc]
        
    nzf = np.loadtxt(LSSdir+tracer+'_'+reg+'_nz.txt').transpose()
    
    bs = nzf[2][0]-nzf[1][0]
    nzd = nzf[3] #column with nbar values
    zl = in_data['Z']
    nl = np.zeros(len(zl))
    zind = ((zl - zmin) / bs).astype(int)
    valid = (zl > zmin) & (zl < zmax) & (zind >= 0) & (zind < len(nzd))
    nl[valid] = nzd[zind[valid]]
    in_data['NX'] = nl*cntl[in_data['NTILE']]
    in_data['WEIGHT_FKP'] = 1/(1+P0*in_data['NX'])
    return in_data
