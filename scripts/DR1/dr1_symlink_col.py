import os
import glob
import sys
import numpy as np
import shutil
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--survey",default='Y1')
parser.add_argument("--pubrelease",default='dr1')
parser.add_argument("--specrel",default='iron')
parser.add_argument("--outroot",default=None)
args = parser.parse_args()


def my_ln(infile, publink, isdir=False):
    if isdir:
        if os.path.isdir(infile):
            os.system('ln -s ' + infile + ' ' + publink)
            #os.system('cp -R ' + infile + ' ' + publink)
            return 'ok'
        else:
            raise Exception(infile + ' does not exist, please revise your script')
    else:
        if os.path.isfile(infile):
            #os.system('cp ' + infile + ' ' + publink)
            os.system('ln -s ' + infile + ' ' + publink)
            return 'ok'
        else:
            raise Exception(infile + ' does not exist, please revise your script')


if args.outroot is None:
    args.outroot = os.getenv('SCRATCH')

indir = '/global/cfs/cdirs/desi/survey/catalogs/'+args.survey+'/LSS/'
pubdir = args.outroot+'/desi/survey/catalogs/'+args.pubrelease+'/'

#MAKE DIRECTORY TREE
temp_path = os.path.join(pubdir, 'LSS')
if not os.path.exists(temp_path):
    os.makedirs(temp_path)

pubdir = temp_path

lssdirs = ['altmtl','iron'] #,'random{RAN}']
num_ran = np.arange(18)

for lssdir in lssdirs:
    temp_path = os.path.join(pubdir,lssdir)
    if not os.path.exists(temp_path):
        os.makedirs(temp_path)
    for nran in num_ran:
        temp_path = os.path.join(pubdir,'random{RAN}').format(RAN=nran)
        if not os.path.exists(temp_path):
            os.makedirs(temp_path)

#MAKE DIRECTORY TREE FOR altmtl
num_rea = np.arange(128)
altmtldirs = ['BRIGHT', 'DARK']
#univdirs = ['fa','main']

for altmtldir in altmtldirs:
    temp_path = os.path.join(pubdir,'altmtl',altmtldir)
    if not os.path.exists(temp_path):
        os.makedirs(temp_path)

    for nrea in num_rea:
        temp_path_univ = os.path.join(pubdir,'altmtl', altmtldir, 'Univ{UNIV}').format(UNIV=str(nrea).zfill(3))
        if not os.path.exists(temp_path_univ):
            os.makedirs(temp_path_univ)
#        for univdir in univdirs:
#            temp_path = os.path.join(temp_path_univ, univdir)
#            if not os.path.exists(temp_path):
#                os.makedirs(temp_path)


#MAKE DIRECTORY TREE FOR iron
temp_path = os.path.join(pubdir,'iron', 'LSScats')
if not os.path.exists(temp_path):
    os.makedirs(temp_path)

newtemp_pubdir = os.path.join(pubdir,'iron', 'LSScats')

print('root directory is '+ newtemp_pubdir)
versions = ['v1.2','v1.5','v1.5pip']
for ver in versions:
    temp_path = os.path.join(newtemp_pubdir, ver)
    if not os.path.exists(temp_path):
        os.makedirs(temp_path)
    if 'pip' not in ver:
        temp_path = os.path.join(newtemp_pubdir, ver, 'hpmaps')
        if not os.path.exists(temp_path):
            os.makedirs(temp_path)

#Symbolic links from LSS dir
lssfnames = ['collisions-BRIGHT.fits','collisions-DARK.fits','tiles-BRIGHT.fits','tiles-DARK.fits']
for fname in lssfnames:

    symindir = os.path.join(indir, fname)
    symoutdir = os.path.join(pubdir, fname)
    my_ln(symindir, symoutdir)

#Symbolic links for randomX dirs
ranlssfnames = ['pota-BRIGHT.fits', 'pota-DARK.fits']
for nran in num_ran:
    
    symindir = os.path.join(indir, 'random{RAN}').format(RAN=nran)
    symoutdir = os.path.join(pubdir, 'random{RAN}').format(RAN=nran)
    for fname in ranlssfnames:
        my_ln(os.path.join(symindir, fname), os.path.join(symoutdir, fname))

#Symbolic links for altmtl
altmtl_keys = {'DARK':'JL_Y1Run2DARK', 'BRIGHT':'JL_Y1Run3BRIGHT'}
files_univ = ['mainsurvey-{PROG}obscon-TileTracker.ecsv', 'tiles-specstatus.ecsv']
for prog in altmtl_keys:
    symindir = os.path.join(indir, 'altmtl', altmtl_keys[prog])
    symoutdir = os.path.join(pubdir, 'altmtl', prog)
    for nrea in num_rea:
        for fname in files_univ:
            my_ln(os.path.join(symindir, 'Univ{UNIV}', fname).format(UNIV=str(nrea).zfill(3), PROG=prog), os.path.join(symoutdir, 'Univ{UNIV}', fname).format(UNIV=str(nrea).zfill(3), PROG=prog)) 
        dirs_univ = ['fa', 'main']
        for dname in dirs_univ:
            my_ln(os.path.join(symindir, 'Univ{UNIV}', dname).format(UNIV=str(nrea).zfill(3)), os.path.join(symoutdir, 'Univ{UNIV}').format(UNIV=str(nrea).zfill(3)), isdir=True)
for root, dirs, files in os.walk(os.path.join(pubdir), topdown=False):
    for dir_name in dirs:
        if dir_name.startswith("."):
            dir_path = os.path.join(root, dir_name)
            if os.path.islink(dir_path):
                print(f"Unlinking directory: {dir_path}")
                os.unlink(dir_path) 
            else:
                print(f"Removing directory: {dir_path}")
            
                shutil.rmtree(dir_path)  # Remove the directory and its contents


#Symbolic links for iron
files_in_iron = np.loadtxt('files_under_iron.txt', unpack=True, dtype=str)
symindir = os.path.join(indir, 'iron')
symoutdir = os.path.join(pubdir, 'iron')
for fname in files_in_iron:
    my_ln(os.path.join(symindir, fname), os.path.join(symoutdir, fname))

#Symbolic links for v1.X
filename = 'files_under_LSScats_{VER}.txt'
for ver in versions:
    files_in_ver = np.loadtxt(filename.format(VER=ver), unpack=True, dtype=str)
    #print(files_in_ver)
    for fname in files_in_ver:
        symindir = os.path.join(indir, 'iron', 'LSScats', ver)
        symoutdir = os.path.join(pubdir, 'iron', 'LSScats', ver)
        my_ln(os.path.join(symindir, fname), os.path.join(symoutdir, fname))

hpmaps_v12 = ['BGS_BRIGHT_mapprops_healpix_nested_nside256_N.fits', 'ELG_LOPnotqso_mapprops_healpix_nested_nside256_N.fits', 'LRG_mapprops_healpix_nested_nside256_N.fits', 'QSO_mapprops_healpix_nested_nside256_N.fits', 'BGS_BRIGHT_mapprops_healpix_nested_nside256_S.fits', 'ELG_LOPnotqso_mapprops_healpix_nested_nside256_S.fits', 'LRG_mapprops_healpix_nested_nside256_S.fits', 'QSO_mapprops_healpix_nested_nside256_S.fits']

symindir = os.path.join(indir, 'iron', 'LSScats', 'v1.2', 'hpmaps')
symoutdir = os.path.join(pubdir, 'iron', 'LSScats', 'v1.2', 'hpmaps')
for fname in hpmaps_v12:
    my_ln(os.path.join(symindir, fname), os.path.join(symoutdir, fname))

hpmaps_v15 = ['BGS_BRIGHT_mapprops_healpix_nested_nside256_N.fits', 'ELG_LOPnotqso_mapprops_healpix_nested_nside256_N.fits', 'LRG_mapprops_healpix_nested_nside256_N.fits', 'QSO_mapprops_healpix_nested_nside256_N.fits', 'BGS_BRIGHT_mapprops_healpix_nested_nside256_S.fits', 'ELG_LOPnotqso_mapprops_healpix_nested_nside256_S.fits', 'LRG_mapprops_healpix_nested_nside256_S.fits', 'QSO_mapprops_healpix_nested_nside256_S.fits', 'BGS_ANY_mapprops_healpix_nested_nside256_N.fits', 'BGS_ANY_mapprops_healpix_nested_nside256_S.fits']

symindir = os.path.join(indir, 'iron', 'LSScats', 'v1.5', 'hpmaps')
symoutdir = os.path.join(pubdir, 'iron', 'LSScats', 'v1.5', 'hpmaps')
for fname in hpmaps_v15:
    my_ln(os.path.join(symindir, fname), os.path.join(symoutdir, fname))




'''
fnames = ['LSS/collisions-BRIGHT.fits','LSS/collisions-DARK.fits','LSS/tiles-BRIGHT.fits','LSS/tiles-DARK.fits']
for fname in fnames:
    my_ln(indir+fname, pubdir+fname)
#    os.system('ln -s '+indir+fname+' ' +pubdir+fname)
my_ln(indir+'/LSS/altmtl/JL_Y1Run2DARK', pubdir+'/LSS/altmtl/DARK', isdir=True)
my_ln(indir+'/LSS/altmtl/JL_Y1Run3BRIGHT', pubdir+'/LSS/altmtl/BRIGHT', isdir=True)

#os.system('ln -s '+indir+'/LSS/altmtl/JL_Y1Run2DARK '+pubdir+'/LSS/altmtl/DARK')
#os.system('ln -s '+indir+'/LSS/altmtl/JL_Y1Run3BRIGHT '+pubdir+'/LSS/altmtl/BRIGHT')

progl = ['DARK','BRIGHT']
ranl = np.arange(0,18)

for rann in ranl:
    randir = pubdir+'/LSS/random'+str(rann)
    if not os.path.exists(randir):
        os.makedirs(randir)

    for prog in progl:
        fname = 'LSS/random'+str(rann)+'/pota-'+prog+'.fits'
        my_ln(indir+fname, pubdir+fname)
        #os.system('ln -s '+indir+fname+' ' +pubdir+fname)
        fname = 'LSS/'+args.specrel+'/rancomb_'+str(rann)+prog.lower()+'wdupspec_zdone.fits'
        my_ln(indir+fname, pubdir+fname)
        #os.system('ln -s '+indir+fname+' ' +pubdir+fname)

fname = 'LSS/'+args.specrel+'/unique_badfibers.txt'
#os.system('ln -s '+indir+fname+' ' +pubdir+fname)
my_ln(indir+fname, pubdir+fname)

tracers = ['BGS_ANY','BGS_BRIGHT','LRG','ELG_LOPnotqso','QSO']
for tr in tracers:
    fname = 'LSS/'+args.specrel+'/datcomb_'+tr+'_tarspecwdup_zdone.fits'
    my_ln(indir+fname, pubdir+fname)
    #os.system('ln -s '+indir+fname+' ' +pubdir+fname)

for prog in progl:
    fname = 'LSS/'+args.specrel+'/datcomb_'+prog.lower()+'_spec_zdone.fits'
    #os.system('ln -s '+indir+fname+' ' +pubdir+fname)
    my_ln(indir+fname, pubdir+fname)
    fname = 'LSS/'+args.specrel+'/datcomb_'+prog.lower()+'_zmtl_zdone.fits'
    #os.system('ln -s '+indir+fname+' ' +pubdir+fname)
    my_ln(indir+fname, pubdir+fname)

fname = 'LSS/'+args.specrel+'/emlin_catalog.fits'
#os.system('ln -s '+indir+fname+' ' +pubdir+fname)
my_ln(indir+fname, pubdir+fname)

for ver in versions:
    fnames = glob.glob(indir+'/LSS/'+args.specrel+'/LSScats/'+ver+'/*.fits')
    for fname in fnames:
        #os.system('ln -s '+fname +' ' +pubdir+fname.replace(indir,''))
        tfname = fname.replace(indir,'')
        my_ln(fname, pubdir+tfname)

    fnames = glob.glob(indir+'/LSS/'+args.specrel+'/LSScats/'+ver+'/*nz.txt')
    for fname in fnames:
        tfname = fname.replace(indir,'')
        #os.system('ln -s '+fname +' ' +pubdir+fname.replace(indir,''))
        my_ln(fname, pubdir+tfname)
    if 'pip' not in ver:
        dirname = '/LSS/'+args.specrel+'/LSScats/'+ver+'/hpmaps'
        #os.system('ln -s '+indir+dirname+' ' +pubdir+dirname)
        my_ln(indir+dirname, pubdir+dirname, isdir=True)
'''
