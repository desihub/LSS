def dorecon(cat,ran,output='_rec',niter=3,padding=200.,smooth=15.,nbins=512,nthreads=1,bias=1.4,f=.82):
    '''
    cat is the data table assumed to have names 'RA','DEC','Z','WEIGHT'
    ran is the random table with the same names
    the rest of the parameters control reconstruction algorithm parameters
    '''
    from recon import Recon

	cat_we = cat['WEIGHT']
	ran_we = ran['WEIGHT']
    opt_box = 1 #optimize box dimensions

    rec = Recon(cat['RA'], cat['DEC'], cat['Z'], cat_we, \
                ran['RA'], ran['DEC'], ran['Z'], ran_we, \
                nbins=nbins, smooth=smooth, f=f, bias=bias, \
                padding=padding, opt_box=opt_box, nthreads=nthreads)
    for i in range(niter):
        rec.iterate(i)
    rec.apply_shifts()
    rec.summary()

    cat['RA'], cat['DEC'], cat['Z'] = rec.get_new_radecz(rec.dat) 
    ran['RA'], ran['DEC'], ran['Z'] = rec.get_new_radecz(rec.ran) 
    cat.write(output+'.dat.fits', format='fits', overwrite=True)
    ran.write(output+'.ran.fits', format='fits', overwrite=True)
