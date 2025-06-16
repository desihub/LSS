import os
import pandas as pd
from astropy.table import Table
from astropy.io import fits
import numpy as np
import glob
#This script modifies SV3 products to meet the datamodel standards. It substitutes the files indicated in the specified directory, 
#Therefore, be careful. Suggestion is to work on a copy of the original files.

#This is the directory where we make the modifications. 
#dir_ = '/pscratch/sd/a/acarnero/dr1_datamodel3/desi/survey/catalogs/dr1'
dir_ = '/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/v1.5pip'
#Here it reads all the directory tree structure
directories = [x[0] for x in os.walk(dir_)]


#This is the columns description file where it read existing columns and types, to compare with the files
filename_to_columns = '/global/homes/a/acarnero/codes/desidatamodel/py/desidatamodel/data/column_descriptions.csv'

#Here it reads the columns and save them in a list for future comparison
columns_official = pd.read_csv(filename_to_columns)
#columns_official.Units.fillna('', inplace=True)
columns_official.fillna({'Units': ''}, inplace=True)
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

check_columns = False
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
check_headers = False


count = 0
for root, dirs, files in os.walk(dir_):
        if 'altmtl' in root:
                continue

        for file in files:
#            print('NAME ----', file)
            if file.endswith("dat.fits"):
                filename = os.path.join(root,file)
                hdul = fits.open(filename)[1]
                #for i in range(len(hdul)-1):
                #    hindex = i+1
                data_to_read = hdul#[hindex]
                if 'BITWEIGHTS' in data_to_read.columns.names or 'PROB_OBS' in data_to_read.columns.names:
                    print(filename)

#                        if y1_old not in names_col_official:
#                            if y1_old == 'Position':
                            #print('column ', y1_old, ' in ',filename,' is not defined in column_descriptions')
                                #print(y1_old,filename,set(data_to_read.data[y1_old]))
#                            if y1_old.upper() in names_col_official:
#                                upper_columns.append(y1_old)
#                                print('column ', y1_old, ' in ',filename,' is defined in capital letters')
#                            else:
#                                count+=1
#                                print('column ', y1_old, ' in ',filename,' is not defined in column_descriptions', data_to_read.data[y1_old].dtype)

#                                new_columns.append(y1_old)
#                hdul.close()
#            else:
#                print('this file is not fits', file)
#                pass


exit()

if check_headers:
    output_file = 'fits_metadata.txt'
    with open(output_file, 'w') as f:
        for root, dirs, files in os.walk(dir_):
            if 'altmtl' in root:
                continue

            for file in files:
                if file.endswith(".fits"):
                    filename = os.path.join(root,file)
                    hdu = fits.open(filename)[1]
                    header = hdu.header
                    non_column_keys = [key for key in header.keys() if not key.startswith(('TTYPE', 'TFORM', 'TUNIT'))]
                    f.write(f"{filename} {' '.join(non_column_keys)}\n")
    #                print(non_column_keys)
    #                for hh in hdu.header.keys():
    #                    print(hh)
    #            exit()
    
upper_columns = []
new_columns = []
''' to check if there are columns not defined in description files '''
if check_columns:
    print('Checking the content of all files for missing columns in the description file')
    
    count = 0
    for root, dirs, files in os.walk(dir_):
        for file in files:
            print('NAME ----', file)
            if file.endswith(".fits"):
                filename = os.path.join(root,file)
                hdul = fits.open(filename)
                for i in range(len(hdul)-1):
                    hindex = i+1
                    data_to_read = hdul[hindex]
                    for y1_old in data_to_read.columns.names:
                        
                        if y1_old not in names_col_official:
#                            if y1_old == 'Position': 
                            #print('column ', y1_old, ' in ',filename,' is not defined in column_descriptions')
                                #print(y1_old,filename,set(data_to_read.data[y1_old]))
                            if y1_old.upper() in names_col_official:
                                upper_columns.append(y1_old)
                                print('column ', y1_old, ' in ',filename,' is defined in capital letters')
                            else:
                                count+=1
                                print('column ', y1_old, ' in ',filename,' is not defined in column_descriptions', data_to_read.data[y1_old].dtype)
                        
                                new_columns.append(y1_old)
                hdul.close()
            else:
                print('this file is not fits', file)
                pass
    if count==0:
        print('All columns in the files are described in column_descriptions')
    else:
        print('Some columns are not described in column_description. Check with the WG if those need to be removed')
        print('or if they need to be add to column_descriptions')
    print('new columns:')
    print(set(new_columns))
    print('columns to be renamed with capital letters')
    print(set(upper_columns))


from astropy.io import fits

if check_extname:
    for root, dirs, files in os.walk(dir_):
        if 'altmtl' in root:
            continue

        for file in files:
            if file.endswith(".fits"):
                filename = os.path.join(root,file)
                hdu = fits.open(filename)[1]
                if 'EXTNAME' in hdu.header:
                    extname = hdu.header['EXTNAME']
                    if extname.strip() == '':  # Check if EXTNAME is an empty string
                        print(f"HDU {filename}: EXTNAME is present but empty.")
                else:
                    if 'tiles' in file:
                        fits.setval(filename, 'EXTNAME', value='TILES', ext=1)
                    elif 'emlin' in file:
                        fits.setval(filename, 'EXTNAME', value='EMLIN', ext=1)
                    else:
                        fits.setval(filename, 'EXTNAME', value='LSS', ext=1)
 
                    print(f"HDU {filename}: EXTNAME does not exist.")

if change_names:
    print('Changing names in files')

    #y1_columns_to_change = { 'lrg_mask':'LRG_MASK', 'o2c': 'O2C', 
    #                         'flux_g_dered':'FLUX_G_DERED', 'flux_r_dered':'FLUX_R_DERED',
    #                         'flux_z_dered':'FLUX_Z_DERED', 'flux_w1_dered':'FLUX_W1_DERED',
    #                         'flux_w2_dered':'FLUX_W2_DERED'
    #                         }

    y1_columns_to_change = { 'LRG_MASK':'lrg_mask', 'O2C': 'o2c', 
                            'FLUX_G_DERED':'flux_g_dered', 'FLUX_R_DERED':'flux_r_dered',
                            'FLUX_Z_DERED':'flux_z_dered', 'FLUX_W1_DERED':'flux_w1_dered',
                            'FLUX_W2_DERED':'flux_w2_dered'
                             }

    count = 0
    for root, dirs, files in os.walk(dir_):
        if 'altmtl' in root:
            continue

        for file in files:
            if file.endswith(".fits"):
                filename = os.path.join(root,file)
                hdul = fits.open(filename,mode='update')
                for i in range(len(hdul)-1):
                    hindex = i+1
                    data_to_read = hdul[hindex]
                    for y1_old in y1_columns_to_change.keys():
                        if y1_old in data_to_read.columns.names:
                            print(filename, y1_columns_to_change[y1_old])
                            count += 1
                            data_to_read.columns[y1_old].name = y1_columns_to_change[y1_old]
                hdul.close()
            else:
                pass
    if count==0:
        print('No changes applied')
    else:
        print(count,' columns changed')

# Remove_columns



def remove_columns_from_fits(input_file, output_file, columns_to_remove):
    # Open the FITS file
    with fits.open(input_file, mode='readonly') as hdul:
        # Access the first table extension (assuming it contains the data)
        hdu = hdul[1]  # Typically the first extension with a table is at index 1
        data = hdu.data
        columns = hdu.columns

        # Filter columns to keep
        columns_to_keep = [col for col in columns if col.name not in columns_to_remove]
        removed_columns = set(columns.names) - set(col.name for col in columns_to_keep)

        if removed_columns:
            print(f"Removing columns: {', '.join(removed_columns)}")
            
            # Create a new column list with the remaining columns
            new_columns = fits.ColDefs(columns_to_keep)
            
            # Create a new HDU with the updated column definitions
            new_hdu = fits.BinTableHDU.from_columns(new_columns)
            
            # Preserve the header
            new_hdu.header.extend(hdu.header, strip=True, update=True)
            
            # Write the updated FITS file
            new_hdul = fits.HDUList([hdul[0], new_hdu])
            new_hdul.writeto(output_file, overwrite=True)
            print(f"New FITS file saved as: {output_file}")
        else:
            print(f"No matching columns found to remove. No changes made.")

# Example usage
if remove_columns:
    print('Removing columns')

    columns_to_remove = [
        'sort', 'SKYMAP_MASK', 'ZDATE', 'Position', 'Z_QN_QF',
        'THRUDATE', 'elg_mask', 'ABSMAG_SDSS_G', 'ABSMAG_SDSS_R', 'COLLISION']
    count = 0

    for root, dirs, files in os.walk(dir_):
        if 'altmtl' in root:
            continue

        for file in files:
            if file.endswith(".fits"):
                filename = os.path.join(root, file)
                print(filename)
#                remove_columns_from_fits(filename, filename, columns_to_remove)
                with fits.open(filename, mode='readonly') as hdul:
                    cols_to_remove = []
                    data_to_read = hdul[1]
                    for colrem in columns_to_remove:
                        if colrem == 'COLLISION' and 'pota' in file:
                            continue
                        if colrem in data_to_read.columns.names:
                            cols_to_remove.append(colrem)

                if len(cols_to_remove) > 0:
                    count += 1
                    remove_columns_from_fits(filename, filename, cols_to_remove)
                    '''
                    print(filename, cols_to_remove)
                    count += 1
                    with fits.open(filename, mode='update') as hdul:
                        print('updating')
                        data_to_read = hdul[1]
                        colsgood = [
                            fits.Column(name=colname, array=data_to_read.data[colname], format=colfmt)
                            for colname, colfmt in zip(data_to_read.columns.names, data_to_read.columns.formats)
                            if colname not in cols_to_remove
                        ]
                        print('end selecting cols to save')
                        hdul[1] = fits.BinTableHDU.from_columns(colsgood)
                        print('end saving')
                    '''
    if count == 0:
        print('No files are affected')
    else:
        print(f'{count} files have been affected')

# Add units to columns in fits files 
if add_units:
    print('Adding units to columns')

    count=0
    for root, dirs, files in os.walk(dir_):
        print(root)
        if 'altmtl' in root:
            continue
        for file in files:
            if file.endswith(".fits"):
                filename = os.path.join(root,file)
                hdul = fits.open(filename, mode='update')
                data_to_read = hdul[1]
                for colindata in data_to_read.columns.names:
                    try:
                       
                        index_column = np.argwhere(np.array(columns_official.Name)==colindata)[0][0]
                        #index_column = np.argwhere(np.array(columns_official.Name)==colindata.upper())[0][0]
                    except IndexError:
                        print(IndexError, file, colindata, columns_official.Name)
                    if data_to_read.columns[colindata].unit != columns_official.Units[index_column] and len(columns_official.Units[index_column])!=0:
                        print(file,colindata)
                        data_to_read.columns[colindata].unit = columns_official.Units[index_column]
                        count+=1
                hdul.close() 
    print(count, ' columns affected')

# Add comments to header 
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

'''# Remove specific columns in files or directory
if remove_specific_columns:
    print('Removing specific columns from file or directory')

    isDir = True
    

    columns_to_remove = {'/global/cfs/cdirs/desi/survey/catalogs/edr_prepfor_public/LSScats/full': ['REF_EPOCH','PARALLAX','PMRA','PMDEC','OBSCONDITIONS','NUMOBS_INIT','NUMOBS_MORE','NUMOBS','ZTILEID','VERSION'], '/global/cfs/cdirs/desi/survey/catalogs/edr_prepfor_public/potential_assignments/data/':['ZTILEID', 'PARALLAX', 'PMRA', 'PMDEC', 'Z', 'REF_EPOCH', 'NUMOBS_INIT', 'NUMOBS', 'NUMOBS_MORE', 'VERSION', 'OBSCONDITIONS']}
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

'''
#Remove specidif header keys from files
if remove_specific_headers:
    print('removing headers from files')
    heads_valuable = ['XTENSION', 'BITPIX', 'NAXIS', 'NAXIS1', 'NAXIS2', 'PCOUNT', 'GCOUNT', 'TFIELDS', 'DESIDR', 'EXTNAME', 'COMMENT']

    '''
    lista = ['RRVER', 'TEMNAM00', 'TEMVER00', 'TEMNAM01', 'TEMVER01', 'TEMNAM02', 'TEMVER02', 'TEMNAM03', 'TEMVER03',
             'TEMNAM04', 'TEMVER04', 'TEMNAM05', 'TEMVER05', 'TEMNAM06', 'TEMVER06', 'TEMNAM07', 'TEMVER07', 'TEMNAM08', 'TEMVER08',
             'TEMNAM09', 'TEMVER09', 'TEMNAM10', 'TEMVER10', 'SPGRP   ', 'SURVEY  ', 'PROGRAM ', 'ZCATVER ']
    lista = ['LONGSTRN']
    lista = ['LSSMAPDIR', 'INFILES']
    lista = ['EBV_CHIANG_SFDCORR','STARDENS', 'HALPHA', 'HALPHA_ERROR','CALIB_G',
             'CALIB_R', 'CALIB_Z', 'EBV_MPF_MEAN_FW15', 'EBV_MPF_MEAN_ZPTCORR_FW15',  
             'EBV_MPF_VAR_FW15','EBV_MPF_VARCORR_FW15','EBV_MPF_MEAN_FW6P1','EBV_MPF_MEAN_ZPTCORR_FW6P1',
             'EBV_MPF_VAR_FW6P1','EBV_MPF_VARCORR_FW6P1', 'EBV_SGF14', 'BETA_ML', 'BETA_MEAN',
             'BETA_RMS', 'HI', 'KAPPA_PLANCK', 'EBV', 'PSFDEPTH_G','PSFDEPTH_R',
             'PSFDEPTH_Z', 'GALDEPTH_G', 'GALDEPTH_R', 'GALDEPTH_Z', 'PSFDEPTH_W1', 'PSFDEPTH_W2',
             'PSFSIZE_G', 'PSFSIZE_R', 'PSFSIZE_Z']
    '''
    #dar_ = os.path.join(dir_,'LSScats/v1.2/hpmaps')
    for root, dirs, files in os.walk(dir_):
        if 'altmtl' in root:# or 'healpix' in root:
            continue

        for file in files:
            if file.endswith('.fits'): #in list(headers_to_remove.keys()):
                filename = os.path.join(root,file)
                hdul = fits.open(filename, mode='update')
                header = hdul[1].header
                non_column_keys = [key for key in header.keys() if not key.startswith(('TTYPE', 'TFORM', 'TUNIT'))]
                for h in non_column_keys:
#                    if h.startswith(('TDIM','TNULL')):
#                        print(filename,h)
                    if h not in heads_valuable:
                        del hdul[1].header[h]
                        print('remove', h)
#                for header in lista:
#                    if header in list(hdul[1].header.keys()):
#                        print(filename,header)
#                        del hdul[1].header[header]

                hdul.close()


if define_extname:
    print('Setting extension name to files')

    extname_files = {'emlim_catalog':[1,'LSS'], 'rancomb':[1,'ZCATALOG'], 'zmtl_zdone':[1,'LSS'], 'tiles-':[1,'LSS']}

#    d=glob.glob('/pscratch/sd/a/acarnero/dr1_datamodel3/desi/survey/catalogs/dr1/LSS/iron/*')
#    for dd in d:
#        if dd.endswith('.fits'):
#            extname_files[dd.split('/')[-1]] = [1,'LSS']

#    d=glob.glob('/global/cfs/cdirs/desi/survey/catalogs/edr_prepfor_public/LSScats/clustering/*')
#    for dd in d:
#        if dd.endswith('.fits'):
#            extname_files[dd.split('/')[-1]] = [1,'LSS']

    for root, dirs, files in os.walk(dir_):
        for file in files:
            for ll in extname_files.keys():
                if ll in file:
                    filename = os.path.join(root,file)
                    print(filename)
                    fits.setval(filename, 'EXTNAME', value=extname_files[ll][1], ext=extname_files[ll][0])
                    print(filename, ' EXTNAME changed')


# Add DESIDR to header for EDR release 
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
                        fits.setval(filename, 'DESIDR', value='dr1', ext=i)

    print(count, 'files affected')

