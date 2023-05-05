import os
import pandas as pd
from astropy.table import Table
from astropy.io import fits
import numpy as np
import glob
#This script modifies SV3 products to meet the datamodel standards. It substitutes the files indicated in the specified directory, 
#Therefore, be careful. Suggestion is to work on a copy of the original files.

#This is the directory where we make the modifications. 
dir_ = '/global/cfs/cdirs/desi/survey/catalogs/edr_prepfor_public/fuji/sv3/v1'

#Here it reads all the directory tree structure
directories = [x[0] for x in os.walk(dir_)]


#This is the columns description file where it read existing columns and types, to compare with the files
filename_to_columns = '/global/homes/a/acarnero/codes/desidatamodel/py/desidatamodel/data/column_descriptions.csv'

#Here it reads the columns and save them in a list for future comparison
columns_official = pd.read_csv(filename_to_columns)
columns_official.Units.fillna('', inplace=True)
names_col_official = list(columns_official.Name)

#Here we set the define the different options of the code:
#The code can: 
#1) Check all files for columns that are not described in the columns description file
#2) Change the name of columns in the files. This was used to capitalize columns that were defined as lowercase in the original files
#3) Change the type of columns in the files. For example, ROSETTE_NUMBER is defined in the original files as floats, when they should be integers
#4) Remove columns from the original files.
#5) Define name to extensions in fits files (this is compulsory for all fits files run through desidatamodel)
#6) Add units to columns in fits files
#7) Add header comment redirecting to readthedocs
#8) Remove columns from specific files or directory
#9) Remove specific header keys from files

check_columns = True
change_names = False
change_types = False
remove_columns = False
define_extname = False
add_units = False
add_comment_header = False
add_desidr_header = False
remove_specific_columns = False
remove_specific_headers = False
check_extname = False

''' to check if there are columns not defined in description files '''
if check_columns:
    print('Checking the content of all files for missing columns in the description file')
    
    count = 0
    for root, dirs, files in os.walk(dir_):
        for file in files:
            if file.endswith(".fits"):
                filename = os.path.join(root,file)
                hdul = fits.open(filename)
                for i in range(len(hdul)-1):
                    hindex = i+1
                    data_to_read = hdul[hindex]
                    for sv3_old in data_to_read.columns.names:
                        
                        if sv3_old not in names_col_official:
                            count+=1
                            print('column ', sv3_old, ' in ',filename,' is not defined in column_descriptions')
                hdul.close()
            else:
                pass
    if count==0:
        print('All columns in the files are described in column_descriptions')
    else:
        print('Some columns are not described in column_description. Check with the WG if those need to be removed')
        print('or if they need to be add to column_descriptions')

'''check files extnames'''
if check_extname:
    for root, dirs, files in os.walk(dir_):
        for file in files:
            if file.endswith(".fits") and not file.startswith("fba"):
                filename = os.path.join(root,file)
                hdul = fits.open(filename)
                print(filename, hdul.info())

''' Change column names '''
if change_names:
    print('Changing names in files')

    sv3_columns_to_change = {'Z_not4clus':'Z', 'lrg_mask':'LRG_MASK', 'o2c':'O2C', 
                             'flux_g_dered':'FLUX_G_DERED', 'flux_r_dered':'FLUX_R_DERED',
                             'flux_z_dered':'FLUX_Z_DERED', 'flux_w1_dered':'FLUX_W1_DERED',
                             'flux_w2_dered':'FLUX_W2_DERED', 'rosette_number':'ROSETTE_NUMBER',
                             'rosette_r':'ROSETTE_R'}
    count = 0
    for root, dirs, files in os.walk(dir_):
        for file in files:
            if file.endswith(".fits"):
                filename = os.path.join(root,file)
                hdul = fits.open(filename,mode='update')
                for i in range(len(hdul)-1):
                    hindex = i+1
                    data_to_read = hdul[hindex]
                    for sv3_old in sv3_columns_to_change.keys():
                        if sv3_old in data_to_read.columns.names:
                            count += 1
                            data_to_read.columns[sv3_old].name = sv3_columns_to_change[sv3_old]
                hdul.close()
            else:
                pass
    if count==0:
        print('No changes applied')
    else:
        print(count,' columns changed')

''' Change column types '''
if change_types:
    print('Changing types in files')

    sv3_column_type_to_change = {'LOCATION_ASSIGNED':bool, 'ROSETTE_NUMBER':'int32'}
    count=0
    for root, dirs, files in os.walk(dir_):
        for file in files:
            if file.endswith(".fits"):
                filename = os.path.join(root,file)
                hdul = fits.open(filename,mode='update')
                for i in range(len(hdul)-1):
                    hindex = i+1
                    data_to_read = hdul[hindex]
                    for sv3_old in sv3_column_type_to_change.keys():
                        if sv3_old in data_to_read.columns.names and data_to_read.data[sv3_old].dtype != sv3_column_type_to_change[sv3_old]:
                            count += 1
                            data_to_read.data[sv3_old] = data_to_read.data[sv3_old].astype(sv3_column_type_to_change[sv3_old])
                hdul.close()
            else:
                pass
    if count==0:
        print('No changes applied')
    else:
        print(count,' columns changed')

''' Remove_columns '''
if remove_columns:
    print('Removing columns')

    columns_to_remove = ['ZERR_QF', 'TSNR2_LYA_QF', 'TSNR2_QSO_QF', 'Z_QN_QF', 'QSO_MASKBITS', 'elg_mask', 'sort', 'LOCFULL', 'FRACFLUX_G', 'FRACFLUX_R', 'FRACFLUX_Z', 'FRACMASKED_G', 'FRACMASKED_R', 'FRACMASKED_Z', 'FRACIN_G', 'FRACIN_R', 'FRACIN_Z', 'ALLMASK_G', 'ALLMASK_R', 'ALLMASK_Z', 'DELTACHI2_HP', 'TSNR2_ELG_HP', 'TSNR2_BGS_HP', 'TSNR2_QSO_HP', 'TSNR2_LRG_HP']
    count=0
    for root, dirs, files in os.walk(dir_):
        for file in files:
            if file.endswith(".fits"):
                filename = os.path.join(root,file)
                hdul = fits.open(filename)#, mode='update')
                cols_to_remove = []
                for i in range(len(hdul)-1):
                    hindex = i+1
                    data_to_read = hdul[hindex]
                    for colrem in columns_to_remove:
                        if colrem in data_to_read.columns.names:
                            cols_to_remove.append(colrem)
                hdul.close()
                if len(cols_to_remove) != 0:
                    count+=1
                    hdul = fits.open(filename, mode='update')
                    for i in range(len(hdul)-1):
                        hindex = i+1
                        data_to_read = hdul[hindex]
                        colsgood = []
                        for colindata,forindata in zip(data_to_read.columns.names,data_to_read.columns.formats):
                            if colindata not in cols_to_remove:
                                colsgood.append(fits.Column(name=colindata, array=data_to_read.data[colindata], format=forindata))
                        hdul[hindex] = fits.BinTableHDU.from_columns(colsgood)
                    hdul.close()
    if count==0:
        print('No files are affected')
    else:
        print(count, ' files have been affected')

''' Add units to columns in fits files '''
if add_units:
    print('Adding units to columns')

    count=0
    for root, dirs, files in os.walk(dir_):
        for file in files:
            if file.endswith(".fits"):
                filename = os.path.join(root,file)
                hdul = fits.open(filename, mode='update')
                for i in range(len(hdul)-1):
                    hindex = i+1
                    data_to_read = hdul[hindex]
                    for colindata in data_to_read.columns.names:
                        index_column = np.argwhere(np.array(columns_official.Name)==colindata)[0][0]
                        if data_to_read.columns[colindata].unit != columns_official.Units[index_column] and len(columns_official.Units[index_column])!=0:
                            print(file,colindata)
                            data_to_read.columns[colindata].unit = columns_official.Units[index_column]
                            count+=1
                hdul.close() 
    print(count, ' coumns affected')

''' Add comments to header '''
if add_comment_header:
    print('Adding Comment to header')

    count=0
    comment = "Visit https://desidatamodel.readthedocs.io/en/latest/ for a more complete description"
    for root, dirs, files in os.walk(dir_):
        for file in files:
            if file.endswith(".fits"):
                filename = os.path.join(root,file)
                hdul = fits.open(filename)
                hh_new = []
                for i in range(len(hdul)-1):
                    hindex = i+1
                    head_temp = hdul[hindex].header
                    to_update = True
                    if 'COMMENT' in list(head_temp.keys()):
                        if isinstance(head_temp['COMMENT'], list):
                            print('is list')
                            list_comments = head_temp['COMMENT']
                            for temp in list_comments:
                                if 'readthedocs' in temp:
                                    to_update = False
                                    continue
                        else:
                            if 'readthedocs' in str(head_temp['COMMENT']):
                                to_update = False
                    if to_update:
                        hh_new.append(hindex)
                    else:
                        print('comment exist already in file', filename)
                hdul.close()
                for i in hh_new:
                    count += 1
                    fits.setval(filename, 'COMMENT', value=comment, ext=i)
    print(count,' headers affected')


''' Remove specific columns in files or directory'''
if remove_specific_columns:
    print('Removing specific columns from file or directory')

    isDir = True
    

    columns_to_remove = {'/global/cfs/cdirs/desi/survey/catalogs/edr_prepfor_public/LSScats/full': ['REF_EPOCH','PARALLAX','PMRA','PMDEC','OBSCONDITIONS','NUMOBS_INIT','NUMOBS_MORE','NUMOBS','ZTILEID','VERSION']}
    count=0
    if isDir:
        for dir_check in columns_to_remove:
            for root, dirs, files in os.walk(dir_check):
                for file in files:
                    if file.endswith(".fits"):
                        filename = os.path.join(root,file)
                        hdul = fits.open(filename)#, mode='update')
                        cols_to_remove = []
                        for i in range(len(hdul)-1):
                            hindex = i+1
                            data_to_read = hdul[hindex]
                            for colrem in columns_to_remove[dir_check]:
                                if colrem in data_to_read.columns.names:
                                    cols_to_remove.append(colrem)
                        hdul.close()
                        if len(cols_to_remove) != 0:
                            count+=1
                            hdul = fits.open(filename, mode='update')
                            for i in range(len(hdul)-1):
                                hindex = i+1
                                data_to_read = hdul[hindex]
                                colsgood = []
                                for colindata,forindata in zip(data_to_read.columns.names,data_to_read.columns.formats):
                                    if colindata not in cols_to_remove:
                                        colsgood.append(fits.Column(name=colindata, array=data_to_read.data[colindata], format=forindata))
                                hdul[hindex] = fits.BinTableHDU.from_columns(colsgood)
                            hdul.close()
    else:
        for filename in columns_to_remove:
            hdul = fits.open(filename)#, mode='update')
            cols_to_remove = []
            for i in range(len(hdul)-1):
                hindex = i+1
                data_to_read = hdul[hindex]
                for colrem in columns_to_remove[filename]:
                    if colrem in data_to_read.columns.names:
                        cols_to_remove.append(colrem)
            hdul.close()
            if len(cols_to_remove) != 0:
                count+=1
                hdul = fits.open(filename, mode='update')
                for i in range(len(hdul)-1):
                    hindex = i+1
                    data_to_read = hdul[hindex]
                    colsgood = []
                    for colindata,forindata in zip(data_to_read.columns.names,data_to_read.columns.formats):
                        if colindata not in cols_to_remove:
                            colsgood.append(fits.Column(name=colindata, array=data_to_read.data[colindata], format=forindata))
                    hdul[hindex] = fits.BinTableHDU.from_columns(colsgood)
                hdul.close()

    if count==0:
        print('No files are affected')
    else:
        print(count, ' files have been affected')


'''Remove specidif header keys from files'''
if remove_specific_headers:
    print('removing headers from files')
    
    headers_to_remove = {'datcomb_dark_tarwdup_Alltiles.fits': ['TILEID', 'TILERA', 'TILEDEC', 'FIELDROT','FA_PLAN','FA_HA','FA_RUN',
                                                                'REQRA','REQDEC','FIELDNUM','FA_VER','FA_SURV','DEPNAM00','DEPVER00',
                                                                'DEPNAM01','DEPVER01','DEPNAM02','DEPVER02','DEPNAM03','DEPVER03',
                                                                'DEPNAM04','DEPVER04','DEPNAM05','DEPVER05','DEPNAM06','DEPVER06',
                                                                'DEPNAM07','DEPVER07','DEPNAM08','DEPVER08'], 
                         'rancomb_darkwdup_Alltiles.fits': ['TILEID', 'TILERA', 'TILEDEC', 'FIELDROT','FA_PLAN','FA_HA','FA_RUN',
                                                                'REQRA','REQDEC','FIELDNUM','FA_VER','FA_SURV','DEPNAM00','DEPVER00',
                                                                'DEPNAM01','DEPVER01','DEPNAM02','DEPVER02','DEPNAM03','DEPVER03',
                                                                'DEPNAM04','DEPVER04','DEPNAM05','DEPVER05','DEPNAM06','DEPVER06',
                                                                'DEPNAM07','DEPVER07','DEPNAM08','DEPVER08'],
                         'rancomb_brightwdup_Alltiles.fits': ['TILEID', 'TILERA', 'TILEDEC', 'FIELDROT','FA_PLAN','FA_HA','FA_RUN',
                                                                'REQRA','REQDEC','FIELDNUM','FA_VER','FA_SURV','DEPNAM00','DEPVER00',
                                                                'DEPNAM01','DEPVER01','DEPNAM02','DEPVER02','DEPNAM03','DEPVER03',
                                                                'DEPNAM04','DEPVER04','DEPNAM05','DEPVER05','DEPNAM06','DEPVER06',
                                                                'DEPNAM07','DEPVER07','DEPNAM08','DEPVER08']
                         }
    input_wspec_list = ['TILEID', 'TILERA', 'TILEDEC', 'FIELDROT','FA_PLAN','FA_HA','FA_RUN',
                                                                'REQRA','REQDEC','FIELDNUM','FA_VER','FA_SURV','DEPNAM00','DEPVER00',
                                                                'DEPNAM01','DEPVER01','DEPNAM02','DEPVER02','DEPNAM03','DEPVER03',
                                                                'DEPNAM04','DEPVER04','DEPNAM05','DEPVER05','DEPNAM06','DEPVER06',
                                                                'DEPNAM07','DEPVER07','DEPNAM08','DEPVER08','LONGSTRN','RRVER','TEMNAM00',
                                                                'TEMVER00','TEMNAM01','TEMVER01','TEMNAM02','TEMVER02','TEMNAM03','TEMVER03',
                                                                'TEMNAM04','TEMVER04','TEMNAM05','TEMVER05','TEMNAM06','TEMVER06','TEMNAM07',
                                                                'TEMVER07','TEMNAM08','TEMVER08','TEMNAM09','TEMVER09','DEPNAM09','DEPVER09',
                                                                'DEPNAM10','DEPVER10','DEPNAM11','DEPVER11','DEPNAM12','DEPVER12','DEPNAM13',
                                                                'DEPVER13','DEPNAM14','DEPVER14','DEPNAM15','DEPVER15','DEPNAM16','DEPVER16',
                                                                'DEPNAM17','DEPVER17','DEPNAM18','DEPVER18','SPGRP','SURVEY','PROGRAM']
    name_prep = 'rancomb_{NUMBER}{PROGRAM}wdupspec_Alltiles.fits'

    for i in range(18):
        for j in ['dark','bright']:
            headers_to_remove[name_prep.format(NUMBER=i, PROGRAM=j)] = input_wspec_list

    headers_to_remove['datcomb_bright_tarspecwdup_Alltiles.fits'] = input_wspec_list
    headers_to_remove['datcomb_dark_tarspecwdup_Alltiles.fits'] = input_wspec_list

    for root, dirs, files in os.walk(dir_):
        for file in files:
            if file in list(headers_to_remove.keys()):
                filename = os.path.join(root,file)
                hdul = fits.open(filename, mode='update')

                for header in headers_to_remove[file]:
                    if header in list(hdul[1].header.keys()):
                        del hdul[1].header[header]

                hdul.close()

''' Add extension name to given files '''
if define_extname:
    print('Setting extension name to files')

    extname_files = {'rancomb_bright_Alltilelocinfo.fits':[1,'TILELOC'], 'rancomb_dark_Alltilelocinfo.fits':[1,'TILELOC']}
    d=glob.glob('/global/cfs/cdirs/desi/survey/catalogs/edr_prepfor_public/LSScats/full/*')
    for dd in d:
        extname_files[dd.split('/')[-1]] = [1,'LSS']

    for root, dirs, files in os.walk(dir_):
        for file in files:
            if file in extname_files.keys(): 
                filename = os.path.join(root,file)
                fits.setval(filename, 'EXTNAME', value=extname_files[file][1], ext=extname_files[file][0])
                print(filename, ' EXTNAME changed')


''' Add DESIDR to header for EDR release '''
if add_desidr_header:
    print('Adding DESIDR to header')
    count=0
    for root, dirs, files in os.walk(dir_):
        for file in files:
            if file.endswith(".fits"):
                filename = os.path.join(root,file)
                hdul = fits.open(filename)
                hh_new = []
                for i in range(len(hdul)-1):
                    hindex = i+1
                    head_temp = hdul[hindex].header
                    to_update = True
                    if 'DESIDR' in list(head_temp.keys()):
                        continue
                    else:
                        hh_new.append(hindex)
                hdul.close()

                if len(hh_new) == 0:
                    continue
                else:
                    count += 1
                    for i in hh_new:
                        fits.setval(filename, 'DESIDR', value='edr', ext=i)

    print(count, 'files affected')
