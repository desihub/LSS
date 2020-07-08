'''
This essentially runs the same things as run_multipass one directory up, but with everything split
out in order to allow more flexibility and tesing
For now, alter the options at the top in order to run with different options
'''

import os

#variables we might want to change

#observing conditions for initial mtl
obscon="DARK|GRAY"
obsconi = [1,2]

#footprint to consider
ramin = 170
ramax = 180
decmin = 0
decmax = 10

fraclya = 0.2 #fraction of quasar targets that we will want to observe 4 times



usedate = "2020-01-01T00:00:00"

target_science_sample='/project/projectdirs/desi/users/ajross/dr8tar/target_science_sample.fits' # AJR wrote out the whole target sample here
target_sky_sample='/project/projectdirs/desi/users/ajross/dr8tar/target_sky_sample.fits'

#some hardcoded things that we need to un-hardcode at some point
usedate = "2020-01-01T00:00:00"

full_target_data='/project/projectdirs/desi/users/ajross/dr8tar/target_science_sample.fits' # AJR wrote out the whole target sample here using tartools.py mktar
sky_data_file='/project/projectdirs/desi/users/ajross/dr8tar/target_sky_sample.fits'

path_to_targets = '/project/projectdirs/desi/target/catalogs/dr8/0.39.0/targets/main/resolve/'

pixweight_file = "/project/projectdirs/desi/target/catalogs/dr8/0.31.1/pixweight/pixweight-dr8-0.31.1.fits"

outdir = '/global/cscratch1/sd/ajross/fiberassigntest/'

os.makedirs(outdir+'targets', exist_ok=True)
os.makedirs(outdir+'footprint', exist_ok=True)

obscon = 'DARK|GRAY'
if obscon == 'DARK|GRAY':
	str_obscon = 'dark_gray'

cap = 'NGC'
dr = 'dr8'

initial_mtl_file = write_initial_mtl_files(cap=cap,dr=dr, ra_min=ramin, ra_max=ramax, dec_min=decmin, dec_max=decmax,outdir=outdir,full_target_data=full_target_data,obscon=obscon,sky_data_file=sky_data_file)

initial_truth_file = "targets/subset_truth_"+dr+"_mtl_"+str_obscon+"_"+cap+".fits"

initial_sky_file = "targets/subset_"+dr+"_sky.fits"

write_initial_sky_file(initial_sky_file,sky_data_file=sky_data_file, ra_min=ramin, ra_max=ramax, dec_min=decmin, dec_max=decmax,outdir=outdir)

write_initial_truth_file(initial_truth_file,initial_mtl_file,outdir,pixweight_file)

    
