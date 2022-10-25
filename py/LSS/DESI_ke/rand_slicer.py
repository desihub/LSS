from astropy.table import Table 

root = '/global/cscratch1/sd/mjwilson/desi/BGS/Sam/'

dat = Table.read(root + 'randoms_bd.fits')  

ss = dat[dat['Z'] >= 0.4]  
ss = ss[ss['Z'] <= 0.5]  

ss.write(root + 'randoms_slice.fits', format='fits', overwrite=True)
