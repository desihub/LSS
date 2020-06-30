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

passes = [0,1,2,3,4]

usedate = "2020-01-01T00:00:00"

target_science_sample='/project/projectdirs/desi/users/ajross/dr8tar/target_science_sample.fits' # AJR wrote out the whole target sample here
target_sky_sample='/project/projectdirs/desi/users/ajross/dr8tar/target_sky_sample.fits'
bdir = '/global/cscratch1/sd/ajross/fiberassigntest/fiducialtargets/temp/' #where outputs get written

sci_input='mtl_science.fits'

import tartools as tt

#toggle to run steps below or not

mkmtli = False #make initial MTL file
if mkmtli:
    tt.mkmtl(obscon=obscon,target_ra_min=ramin,target_ra_max=ramax,target_dec_min=decmin,target_dec_max=decmax,outdir=bdir,target_sample=target_science_sample)
    tt.mkmtl_sky(target_ra_min=ramin,target_ra_max=ramax,target_dec_min=decmin,target_dec_max=decmax,outdir=bdir,target_sample=target_sky_sample)
    tt.add_lya(frac=fraclya,indir=bdir)

mktiles = False #make the tile files
if mktiles:
   tt.mktilefile(obscon=obsconi,target_ra_min=ramin,target_ra_max=ramax,target_dec_min=decmin,target_dec_max=decmax,outdir=bdir)

runsurvey = True
if runsurvey:
    for ps in passes:
        footprint = 'tile_'+str(ps)+'.fits'
        run_assignment(footprint, assign_date = usedate, indir=bdir)
        obs, hist_tgassign, hist_tgavail, hist_tgconsid, hist_tgfrac = assignment_counts(footprint, science_input=sci_input, fba_dir='fiberassign/',indir=bdir)
        oldf = 'mtl_science_pass'+str(ps)+'.fits'
        update_mtl(obs,oldf=oldf,science_input=sci_input,indir=bdir )

       