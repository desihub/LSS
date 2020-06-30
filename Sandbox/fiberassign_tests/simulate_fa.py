'''
This runs the functions in tartools.py, as necessary, in order to simulate the fiber assignment
'''

#variables we might want to change

#observing conditions for initial mtl
obscon="DARK|GRAY"
obsconi = [1,2]

#footprint to consider
ramin = 0
ramax = 10
decmin = 0
decmax = 10

fraclya = 0.2 #fraction of quasar targets that we will want to observe 4 times

target_science_sample='/project/projectdirs/desi/users/ajross/dr8tar/target_science_sample.fits' # AJR wrote out the whole target sample here
target_sky_sample='/project/projectdirs/desi/users/ajross/dr8tar/target_sky_sample.fits'
bdir = '/global/cscratch1/sd/ajross/fiberassigntest/fiducialtargets/temp/' #where outputs get written

import tartools as tt

#toggle to run steps below or not

mkmtli = True #make initial MTL file
if mkmtli:
    tt.mkmtl(obscon=obscon,target_ra_min=ramin,target_ra_max=ramax,target_dec_min=decmin,target_dec_max=decmax,outdir=bdir,target_sample=target_science_sample)
    tt.mkmtl_sky(target_ra_min=ramin,target_ra_max=ramax,target_dec_min=decmin,target_dec_max=decmax,outdir=bdir,target_sample=target_sky_sample)
    tt.add_lya(frac=fraclya,indir=bdir)

mktiles = True #make the tile files
if mktiles:
   tt.mktilefile(obscon=obsconi,target_ra_min=ramin,target_ra_max=ramax,target_dec_min=decmin,target_dec_max=decmax,outdir=bdir)
        