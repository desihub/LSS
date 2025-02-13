import os
import glob
import sys
import numpy as np

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
            os.system("rsync -av --ignore-existing --exclude '*.pickle' --exclude '*.sh' " + infile + ' ' + publink)
            #os.system('cp -R ' + infile + ' ' + publink)
            return 'ok'
        else:
            raise Exception(infile + ' does not exist, please revise your script')
    else:
        if os.path.isfile(infile):
            os.system('rsync -av --ignore-existing ' + infile + ' ' + publink)
            ##os.system('ln -s ' + infile + ' ' + publink)
            return 'ok'
        else:
            raise Exception(infile + ' does not exist, please revise your script')


if args.outroot is None:
    args.outroot = os.getenv('SCRATCH')

'''
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


for root, dirs, files in os.walk(os.path.join(pubdir, 'altmtl')):
    for file in files:
        if file.endswith(".pickle"):
            file_path = os.path.join(root, file)
            print(f"Deleting: {file_path}")
            os.remove(file_path)

#Symbolic links for iron
files_in_iron = np.loadtxt('files_under_iron.txt', unpack=True, dtype=str)
symindir = os.path.join(indir, 'iron')
symoutdir = os.path.join(pubdir, 'iron')
for fname in files_in_iron:
    my_ln(os.path.join(symindir, fname), os.path.join(symoutdir, fname))

#Symbolic links for v1.X
filename = 'files_under_LSScats_{VER}.txt'
#v1.2 are all physical files:
files_in_ver = np.loadtxt(filename.format(VER='v1.2'), unpack=True, dtype=str)
for fname in files_in_ver:
        symindir = os.path.join(indir, 'iron', 'LSScats', ver)
        symoutdir = os.path.join(pubdir, 'iron', 'LSScats', ver)
        my_ln(os.path.join(symindir, fname), os.path.join(symoutdir, fname))

#v1.5 and v1.5pip are both physical and symbolic links from previous versions
for ver in ['v1.5','v1.5pip']:
    files_in_ver = np.loadtxt(filename.format(VER=ver), unpack=True, dtype=str)
    #print(files_in_ver)
    for fname in files_in_ver:
        symindir = os.path.join(indir, 'iron', 'LSScats', ver)
        symoutdir = os.path.join(pubdir, 'iron', 'LSScats', ver)
        if os.path.islink(os.path.join(symindir, fname)):
            real_path = os.readlink(os.path.join(symindir, fname))
            real_path_abs = os.path.abspath(os.path.join(os.path.dirname(os.path.join(symindir, fname)), real_path))
            #print(real_path_abs)
            if os.path.islink(real_path_abs):
                real_path = os.readlink(real_path_abs)
                real_path_abs = os.path.abspath(os.path.join(os.path.dirname(real_path_abs), real_path))
                if os.path.islink(real_path_abs):
                    real_path = os.readlink(real_path_abs)
                    real_path_abs = os.path.abspath(os.path.join(os.path.dirname(real_path_abs), real_path))
#            
            my_ln(real_path_abs, os.path.join(symoutdir, fname))
        else:

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
    if os.path.islink(os.path.join(symindir, fname)):
        real_path = os.readlink(os.path.join(symindir, fname))
        real_path_abs = os.path.abspath(os.path.join(os.path.dirname(os.path.join(symindir, fname)), real_path))
        my_ln(real_path_abs, os.path.join(symoutdir, fname))
    else:
        my_ln(os.path.join(symindir, fname), os.path.join(symoutdir, fname))

my_ln('/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/v1.3pip/DARK_bitweights.fits', os.path.join(pubdir, 'iron', 'LSScats', 'v1.5pip', 'DARK_bitweights.fits'))
my_ln('/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/v1.3pip/BRIGHT_bitweights.fits', os.path.join(pubdir, 'iron', 'LSScats', 'v1.5pip', 'BRIGHT_bitweights.fits'))

'''
#QUASAR CATALOG 

indir = '/global/cfs/cdirs/desi/survey/catalogs/'+args.survey+'/QSO/'
pubdir = args.outroot+'/desi/survey/catalogs/'+args.pubrelease+'/'

#MAKE DIRECTORY TREE
temp_path = os.path.join(pubdir, 'QSO')
if not os.path.exists(temp_path):
    os.makedirs(temp_path)

temp_path = os.path.join(temp_path, 'iron')
if not os.path.exists(temp_path):
    os.makedirs(temp_path)

pubdir = temp_path

my_ln(os.path.join(indir, 'iron', 'QSO_cat_iron_cumulative_v0.fits'), os.path.join(pubdir, 'QSO_cat_iron_cumulative_v0.fits'))

#In v1.2, the actual clustering catalogs are in unblinded directory. You need to copy them as well
#Also remove from v1.2, v1.5 and v1.5pip any file that is not NGC or SGC

