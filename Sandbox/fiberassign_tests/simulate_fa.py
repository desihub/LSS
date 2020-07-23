'''
This runs the functions in tartools.py, as necessary, in order to simulate the fiber assignment
'''

#variables we might want to change

#observing conditions for initial mtl
obscon="DARK|GRAY"
obsconi = [1,2]

#footprint to consider
ramin = 0
ramax = 40
decmin = 0
decmax = 40

fraclya = 1 #fraction of quasar targets that we will want to observe 4 times



usedate = "2020-01-01T00:00:00"

target_science_sample='/project/projectdirs/desi/users/ajross/dr8tar/target_science_sample.fits' # AJR wrote out the whole target sample here
target_sky_sample='/project/projectdirs/desi/users/ajross/dr8tar/target_sky_sample.fits'

fracgray = 0

if fracgray == 0:
    mode = 'fiducialtargets'

if fracgray == 0.5:
    mode = 'ELG5050'


if fracgray == 0.4:
    mode = 'ELG4060'

if fracgray == 0.3:
    mode = 'ELG3070'
    
if fracgray == 0.25:
    mode = 'ELG2575'      
    
if fracgray == 0.2:
    mode = 'ELG2080'    

if fracgray == 0.1:
    mode = 'ELG1090' 
    
passes = [0,1,2,3,4]

graylast = False
if graylast:
    passes = [1,2,3,4,0]
    if fracgray == 0:
        mode = 'graylastfid'   

grayind = False
if grayind == True:
    if fracgray == 0:
        mode = 'grayindfid' 

if fraclya == 1:
    mode = 'lya1'                

print('mode is '+mode)

bdir = '/global/cscratch1/sd/ajross/fiberassigntest/'+mode+'/temp/' #base directory for outputs and then downstream inputs

sci_input='mtl_science.fits'

fullfoot = False

import tartools as tt

import fitsio
from matplotlib import pyplot as plt

#toggle to run steps below or not

mkmtli = True #make initial MTL file
if mkmtli:
    tt.mkmtl(obscon=obscon,target_ra_min=ramin,target_ra_max=ramax,target_dec_min=decmin,target_dec_max=decmax,outdir=bdir,target_sample=target_science_sample)
    if fracgray != 0:
        tt.splitdarkgray(fracgray,indir=bdir)
    print('science mtl done')
    if fullfoot != True:
    	tt.mkmtl_sky(target_ra_min=ramin,target_ra_max=ramax,target_dec_min=decmin,target_dec_max=decmax,outdir=bdir,target_sample=target_sky_sample)
    	print('sky mtl done')
    tt.add_lya(frac=fraclya,indir=bdir)
    print('lya added to science mtl')

mktiles = True #make the tile files
if mktiles:
   tt.mktilefile(obscon=obsconi,target_ra_min=ramin,target_ra_max=ramax,target_dec_min=decmin,target_dec_max=decmax,outdir=bdir)

runsurvey = True
if runsurvey:
    for ps in passes:
        footprint = 'tile_'+str(ps)+'.fits'
        tt.run_assignment(footprint, assign_date = usedate, indir=bdir,fullfoot=fullfoot,fullsky=target_sky_sample)
        print('\n')
        print('FINISHED PASS '+str(ps))
        obs, hist_tgassign, hist_tgavail, hist_tgconsid, hist_tgfrac = tt.assignment_counts(footprint, science_input=sci_input, fba_dir='fiberassign/',indir=bdir)
        oldf = 'mtl_science_pass'+str(ps)+'.fits'
        up = True
        if ps == 0 and grayind == True:
            up = False
        if up == True:
            print('updating MTL after pass '+str(ps))
            tt.update_mtl(obs,oldf=oldf,science_input=sci_input,indir=bdir )

plotnumobs = True
if plotnumobs:
    f = fitsio.read(bdir+'mtl_science.fits')
    plt.scatter(f['RA'],f['DEC'],c=f['NUMOBS_MORE'],s=.1)       
    plt.show()

getstats = True
if getstats:
    tt.get_mtlstats(bdir,passes)    