import os
import pandas as pd
from astropy.table import Table
from astropy.io import fits
import numpy as np

#This script modifies SV3 products to meet the datamodel standards. It substitutes the files indicated in the specified directory, 
#Therefore, be careful. Suggestion is to work on a copy of the original files.

#This is the directory where we make the modifications. 
dir_ = '/pscratch/sd/a/acarnero/tobereleased'

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

check_columns = False
change_names = True
change_types = True
remove_columns = True
define_extname = True
add_units = True
add_comment_header = True

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

    columns_to_remove = ['ZERR_QF', 'TSNR2_LYA_QF', 'TSNR2_QSO_QF', 'Z_QN_QF', 'QSO_MASKBITS', 'elg_mask', 'sort', 'LOCFULL', 'FRACFLUX_G', 'FRACFLUX_R', 'FRACFLUX_Z', 'FRACMASKED_G', 'FRACMASKED_R', 'FRACMASKED_Z', 'FRACIN_G', 'FRACIN_R', 'FRACIN_Z', 'ALLMASK_G', 'ALLMASK_R', 'ALLMASK_Z', 'Z_HP', 'DELTACHI2_HP', 'TSNR2_ELG_HP', 'TSNR2_BGS_HP', 'TSNR2_QSO_HP', 'TSNR2_LRG_HP']
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

''' Add extension name to given files '''
if define_extname:
    print('Setting extension name to files')

    extname_files = {'rancomb_bright_Alltilelocinfo.fits':[1,'TILELOC'], 'rancomb_dark_Alltilelocinfo.fits':[1,'TILELOC']}
    for root, dirs, files in os.walk(dir_):
        for file in files:
            if file in extname_files.keys(): 
                filename = os.path.join(root,file)
                fits.setval(filename, 'EXTNAME', value=extname_files[file][1], ext=extname_files[file][0])
                print(filename, ' EXTNAME changed')

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


