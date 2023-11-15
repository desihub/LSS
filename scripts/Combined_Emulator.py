# Run this script from LSS scripts directory
import argparse
import pandas as pd
import numpy as np
import os
import subprocess
from datetime import datetime
from astropy.io import (fits, ascii)
from numpy.random import default_rng
import pickle
from scipy.io import FortranFile
from astropy.table import Table, join, vstack, unique
import LSS.common_tools as ct


startTime = datetime.now()
#python /global/homes/s/sikandar/Combined_Emulator.py -mockver ab_secondgen_cosmosim --real 1 --prog bright --base_output /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummitBGS/ --tracer BGS --emulator_dir /global/cfs/cdirs/desi/survey/catalogs/main/mocks/FAemu_preliminary/sikandar/Updated_Code_CFC/ --galcap B

parser = argparse.ArgumentParser()
parser.add_argument("--mockver", help="type of mock to use",default=None)
parser.add_argument("--mockpath", help="Location of mock file(s)",default='/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/AbacusSummit/CutSky/')
parser.add_argument("--mockfile", help="formattable name of mock file(s). e.g. cutsky_{TYPE}_{Z}_AbacusSummit_base_c000_ph{PH}.fits. TYPE will be replaced with tracer type. PH will be replaced with realization number for simulation of mock.",default='cutsky_{TYPE}_{Z}_AbacusSummit_base_c000_ph{PH}.fits')
parser.add_argument("--real", help="number for the realization",default=0)
parser.add_argument("--prog", help="dark or bright",default='dark')
parser.add_argument("--base_output", help="base directory for output",default='/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/')
parser.add_argument("--tracer", help="which tracer to do")
parser.add_argument("--emulator_dir", help="base directory of ffa emulator")
parser.add_argument("--galcap")
parser.add_argument("--prep", help = "do preprocessing if y", default='n')
parser.add_argument("--emulate", help = "do emulation if y", default='n')
args = parser.parse_args()

# if args.mockver == "ab_secondgen_cosmosim":
#     add_string = "cosmosim/"
# else:
#     add_string = ""
add_string = ""


getpota_indir = args.base_output + add_string 
getpota_outdir = args.base_output + add_string
NTILE_assign_indir = args.base_output + add_string + "SecondGenMocks/AbacusSummit/mock" + str(args.real) + "/"
DESIwemu_indir = args.base_output + add_string + "SecondGenMocks/AbacusSummit/" 
DESIwemu_outdir = args.emulator_dir + "fof_v1.0/in/"

cmd_string1 = "python /global/homes/s/sikandar/PerlmutterLSS/LSS/scripts/mock_tools/prepare_mocks_Y1.py --rbandcut 19.5 --mockver %s --prog %s --downsampling n --overwrite 1 --realmin %s --realmax %s --base_output %s" %(args.mockver, args.prog.lower(), args.real, int(args.real) + 1, args.base_output)
cmd_string2 = "python /global/homes/s/sikandar/PerlmutterLSS/LSS/scripts/getpotaY1_mock.py --base_output %s --prog %s --realization %s --base_input %s"%(getpota_indir, args.prog.upper(), args.real, getpota_outdir)
cmd_string3 = "python /global/homes/s/sikandar/PerlmutterLSS/LSS/scripts/mock_tools/NTILE_assign.py --indir %s --tracer %s --tileloc_add y --prog %s"%(NTILE_assign_indir, args.tracer, args.prog.upper())
cmd_string4 = "python /global/cfs/cdirs/desi/survey/catalogs/main/mocks/FAemu_preliminary/sikandar/Updated_Code_CFC/DESIAbacuswemu.py --mocknumber %s --tracer %s --basedir %s --prog %s --overwrite y --outdir %s --mocktype %s"%(args.real, args.tracer, DESIwemu_indir, args.prog.upper(),DESIwemu_outdir, args.mockver)
print(cmd_string1)
print(cmd_string2)
print(cmd_string3)
print(cmd_string4)

if args.prep == 'y':
    subprocess.run(cmd_string1, shell = True)
    print("Done with prepare_mocks_y1")
    subprocess.run(cmd_string2, shell = True)
    print("Done with getpotaY1_mock")
    subprocess.run(cmd_string3, shell = True)
    print("done with NTILE_assign")
    subprocess.run(cmd_string4, shell = True)
    print("all done")
else:
    print("skipped preprocessing")
print("FFA preprocessing time:")
print(datetime.now() - startTime)

# EmulatorRun

tardict = {4: 'QSO', 1: 'LRG', 34: 'ELG', 60: 'BGS'}
emulator_basedir = args.emulator_dir
fofin_basedir = emulator_basedir + "fof_v1.0/in/"
fofin_file = args.tracer + "_" + args.mockver + "_forFAemu_m" + args.real + ".dat" 
inpath = fofin_basedir + fofin_file
intable = pd.read_csv(inpath, delimiter = ' ', header = None)
new_namin = inpath
new_npart = len(intable)
new_tartype0 = np.unique(intable[4])[0]
new_target = tardict[new_tartype0]
new_galcap = args.galcap
out_dir = emulator_basedir + "fof_v1.0/out/" + args.mockver + "/"
emu_out_dir = emulator_basedir + "emulate_bitw_v1.1/out/" + args.mockver + "/"
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
    print('made ' + out_dir)
if not os.path.exists(emu_out_dir):
    os.makedirs(emu_out_dir)
    print('made ' + emu_out_dir)

fof_prefix = emulator_basedir + "fof_v1.0/"

new_namout_part = fof_prefix + 'out/' + args.mockver + '/' + args.mockver + '_forFAemu_m' + args.real +'_' + new_galcap + '_' + new_target + '_truent_part_llen0.02.dat'
new_namout_halo = fof_prefix + 'out/' + args.mockver + '/' + args.mockver + '_forFAemu_m' + args.real +'_' + new_galcap + '_' + new_target + '_truent_halo_llen0.02.dat'

# Define the file path
file_path = emulator_basedir + "fof_v1.0/INI_fof_lightcone_ang.txt"
# Read the file and modify the parameters
with open(file_path, "r") as file:
    lines = file.readlines()


# Modify the parameter values
for i in range(len(lines)):
    if lines[i].startswith("namin"):
        lines[i] = f"namin '{new_namin}'\n"
    elif lines[i].startswith("namout_part"):
        lines[i] = f"namout_part '{new_namout_part}'\n"
    elif lines[i].startswith("namout_halo"):
        lines[i] = f"namout_halo '{new_namout_halo}'\n"
    elif lines[i].startswith("npart"):
        lines[i] = f"npart {new_npart}\n"
    # elif lines[i].startswith("linklen") and args.linklen:
    #     lines[i] = f"linklen       {new_linklen}\n"
    elif lines[i].startswith("galcap"):
        lines[i] = f"galcap       {new_galcap}    !> 'N' ~~> North;   'S' ~~> South;   'B' (or anything else) ~~> Both\n"
    elif lines[i].startswith("tartype0"):
        lines[i] = f"tartype0       {new_tartype0}     !> 34 ~~> ELG;   1 ~~> LRG;   4 ~~> QSO 60 ~~> BGS\n"
    elif lines[i].startswith("fof_dir"):
        lines[i] = f"fof_dir '{emulator_basedir + 'fof_v1.0/'}'\n"

file_path_new_fof = emulator_basedir + "fof_v1.0/INI_fof_lightcone_ang_%s_m%s.txt"%(args.mockver, args.real)
# Write the modified content back to the same file
with open(file_path_new_fof, "w") as file:
    file.writelines(lines)





## emulate
if new_target == "BGS":
    selmskfile = "forFAemu_data_BGS.dat"
else:
    selmskfile = "forFAemu_data.dat"

emu_naminpart = new_namout_part
emu_naminselmsk = emulator_basedir + "fof_v1.0/in/" + selmskfile
emu_pfxselprop = emulator_basedir + "emulate_bitw_v1.1/out/" + args.mockver + "/" + "FAemu_m" + args.real + "_" + args.galcap + "_" + new_target + "_llen0.02_truent_cfc_"
emu_naminqref = emulator_basedir + "emulate_bitw_v1.1/out/TrainedEmulator_Data/FAemu_data_" + args.galcap + "_" + new_target + "_llen0.02_truent_cfc_antico_deco0.1_beta0_indselfrac_cfc_qm.dat"
emu_pfxout = emulator_basedir + "emulate_bitw_v1.1/out/" + args.mockver + "/FAemu_m" + args.real + "_" + args.galcap + "_" + new_target + "_llen0.02_truent_cfc_nbits217_"
emu_npart = len(intable)

# Define the file path
file_path = emulator_basedir + "/emulate_bitw_v1.1/INI_emulate_bitw_lightcone_ntile.txt"

# Read the file and modify the parameters
with open(file_path, "r") as file:
    lines = file.readlines()

# Modify the parameter values

#print(lines)

for i in range(len(lines)):
    if lines[i].startswith("namin_part"):
        lines[i] = f"namin_part       '{emu_naminpart}'\n"
    elif lines[i].startswith("namin_selmsk"):
        lines[i] = f"namin_selmsk     '{emu_naminselmsk}'\n"
    elif lines[i].startswith("pfx_selprop"):
        lines[i] = f"pfx_selprop      '{emu_pfxselprop}'\n"
    elif lines[i].startswith("namin_qref"):
        lines[i] = f"namin_qref       '{emu_naminqref}'\n"
    elif lines[i].startswith("pfx_out"):
        lines[i] = f"pfx_out          '{emu_pfxout}'\n"
    elif lines[i].startswith("npart"):
        lines[i] = f"npart              {emu_npart}\n"
    elif lines[i].startswith("emu_dir"):
        lines[i] = f"emu_dir     '{emulator_basedir + 'emulate_bitw_v1.1/'}'\n"

# Write the modified content back to the same file
file_path_new_emu = emulator_basedir + "/emulate_bitw_v1.1/INI_emulate_bitw_lightcone_ntile_%s_%s.txt"%(args.mockver, args.real)
with open(file_path_new_emu, "w") as file:
    file.writelines(lines)


if args.emulate == 'n':
    print("skipped emulation")

    
elif args.emulate == 'y':
    print("starting fof")
    cmd_string5 = "sh " + emulator_basedir + "fof_v1.0/COMPILE.sh"
    cmd_string6 = emulator_basedir + "fof_v1.0/fof_lightcone_ang " + file_path_new_fof
    outputfn = emulator_basedir + "fof_v1.0/log/" + args.mockver + "_fof_" + new_target + "_llen0.02_" + args.galcap + "_m" + args.real + "_datatrain.txt"
    errorfn = emulator_basedir + "fof_v1.0/log/" + "ERR_" + args.mockver + "_fof_" + new_target + "_llen0.02_" + args.galcap + "_m" + args.real + "_datatrain.txt"
    print(cmd_string5)
    print(cmd_string6)
    with open(outputfn, "w+") as output, open(errorfn, "w+") as errorf:
        subprocess.run(cmd_string5, shell=True)
        subprocess.run(cmd_string6, shell=True, stdout=output, stderr = errorf)
    print("done fof")


    print("starting emulation")
    cmd_string7 = "sh " + emulator_basedir + "emulate_bitw_v1.1/COMPILE.sh"
    cmd_string8 = emulator_basedir + "emulate_bitw_v1.1/emulate_bitw_lightcone_ntile " + file_path_new_emu
    outputfnemu = emulator_basedir + "emulate_bitw_v1.1/log/" + args.mockver + "_emu_" + new_target + "_llen0.02_" + args.galcap + "_m" + args.real + "_datatrain.txt"
    errorfnemu = emulator_basedir + "emulate_bitw_v1.1/log/" + "ERR_" + args.mockver + "_emu_" + new_target + "_llen0.02_" + args.galcap + "_m" + args.real + "_datatrain.txt"
    print(cmd_string7)
    print(cmd_string8)
    with open(outputfnemu, "w+") as output, open(errorfnemu, "w+") as errorf:
        subprocess.run(cmd_string7, shell=True)
        subprocess.run(cmd_string8, shell=True, stdout=output, stderr = errorf)
    print("done emulation")





    # WEIGHTS
    nwe = 9
    namin_fbw = emulator_basedir + "emulate_bitw_v1.1/out/" + args.mockver + "/FAemu_m" + args.real + "_" + args.galcap + "_" + new_target + "_llen0.02_truent_cfc_nbits217_" + "wemu_unformatted.dat"

    f = FortranFile(namin_fbw, 'r')
    we = f.read_ints(dtype='int32')
    npart = np.shape(we)[0]//nwe
    print('npart =', npart)
    we = we.reshape((nwe, npart)).T;
    print(np.shape(we))

    print(we)


    # -

    # # Fuctions to pack-unpack bitweights in different formats

    def pack_bitweights(array):
        """
        Creates an array of bitwise weights stored as 64-bit signed integers
        Input: a 2D boolean array of shape (Ngal, Nreal), where Ngal is the total number 
               of target galaxies, and Nreal is the number of fibre assignment realizations.
        Output: returns a 2D array of 64-bit signed integers. 
        """
        Nbits=64
        dtype=np.int64
        Ngal, Nreal = array.shape           # total number of realizations and number of target galaxies
        Nout = (Nreal + Nbits - 1) // Nbits # number of output columns
        # intermediate arrays
        bitw8 = np.zeros((Ngal, 8), dtype="i")   # array of individual bits of 8 realizations
        bitweights = np.zeros(Ngal, dtype=dtype) # array of 64-bit integers
        # array to store final output
        output_array = np.zeros((Ngal, Nout), dtype=dtype)
        idx_out = 0 # initial column in output_array
        # loop through realizations to build bitwise weights
        for i in range(Nreal):
            bitw8[array[:,i], i%8] = 1
            arr = np.array(np.packbits(bitw8[:,::-1]), dtype=dtype)
            bitweights = np.bitwise_or(bitweights, np.left_shift(arr, 8*((i%Nbits)//8)))
            #print(np.binary_repr(bitweights[0]), bitweights[0])
            if (i+1)%Nbits == 0 or i+1 == Nreal:
                output_array[:,idx_out] = bitweights
                bitweights[:] = 0
                idx_out += 1
            if (i+1)%8 == 0:
                bitw8[:] = 0
        return output_array


    def unpack_bitweights(we):
        Nbits = 64
        Ngal, Nwe = np.shape(we)
        Nreal = Nbits*Nwe
        print('Nbits, Nwe = ',Nbits,Nwe)
        print('Nreal = ',Nreal)
        print('Ngal = ',Ngal)
        true8=[np.uint8(255) for n in range(0, Ngal)]
        array_bool = np.zeros((Ngal,Nreal), dtype=bool)
        for j in range(Nwe):
            lg = np.zeros((Ngal, Nbits), dtype=bool)
            for i in range(Nbits//8):
                chunk8 = np.uint8(np.bitwise_and(np.right_shift(we[:,j],8*i), true8))
                lg[:,Nbits-8*(i+1):Nbits-i*8] = np.reshape(np.unpackbits(chunk8), (Ngal, 8))
            array_bool[:,j*Nbits:(j+1)*Nbits] = lg[:,::-1]
        return array_bool


    def pack_bitweights_BP2017(array):
        Nbits = 31
        Nwe = np.shape(array)[1]//Nbits
        Nreal = Nwe*Nbits
        Ngal = np.shape(array)[0]
        print('Nbits, Nwe = ',Nbits,Nwe)
        print('Nreal = ',Nreal)
        print('Ngal = ',Ngal)
        we = np.zeros((Ngal, Nwe), dtype=np.int32)
        for i in range(Nreal):
            k = i % Nbits
            m = i // Nbits
            we[array[:,i], m] += 2**k        
        return we


    # +
    def unpack_bitweights_BP2017(we):
        Nbits = 31
        Ngal, Nwe = np.shape(we)
        Nreal = Nbits * Nwe
        print('Nbits, Nwe = ',Nbits,Nwe)
        print('Nreal = ',Nreal)
        print('Ngal = ',Ngal)    
        array_bool = np.zeros((Ngal, Nreal), dtype=bool)

        # To speed up
        we_rs = np.reshape(we, (Ngal, Nwe, 1))

        array_bool = np.bitwise_and(np.right_shift(we_rs, np.arange(Nbits)), 1).astype(bool)
        array_bool = np.reshape(array_bool, (Ngal, Nreal))

        return array_bool



    # ## Unpacking bitweights

    bool_tot = unpack_bitweights_BP2017(we)

    emucat_outdir = emulator_basedir + "emulate_bitw_v1.1/out/" + 'emucat_out/'
    if not os.path.exists(emucat_outdir):
        os.makedirs(emucat_outdir)
        print('made ' + emucat_outdir)

    namout_bool = emucat_outdir + args.mockver + '_' + args.tracer + '_' + args.galcap + '_bool_tot_datatrain_m' + args.real + '.bin'

    # ### write bool array (just for convenience, not to have to recompute it)

    bool_file = open(namout_bool, 'wb')
    pickle.dump(bool_tot, bool_file)
    bool_file.close()

    # ### read bool array (just for convenience, in case you have already computed it)

    # +
    namin_bool = namout_bool

    bool_file = open(namin_bool, 'rb')
    bool_tot = pickle.load(bool_file)
    bool_file.close()

    #print(bool_tot)
    # -



    # +
    #bitw_test = pack_bitweights_BP2017(bool_tot)
    # -

    # ## Randomly extracting "observed" sample and alt-realisations

    rng = default_rng()
    nsmpls_out = 64*4
    nsmpls_in = 31*9 
    id_ran = rng.choice(nsmpls_in, size=nsmpls_out+1, replace=False)
    print(id_ran)

    msk_obs = bool_tot[:,id_ran[0]]
    id_msk_obs = np.where(msk_obs)
    print(msk_obs)
    print(np.sum(msk_obs))
    print(id_msk_obs)
    print(np.shape(id_msk_obs))

    bool_out = bool_tot[:,id_ran[1:]]
    print(np.shape(bool_out))

    # ## Packing bitweights with the desired format

    # +
    we_out = pack_bitweights(bool_out)
    #we_out = pack_bitweights_BP2017(bool_out)

    print(we_out)
    print(np.shape(we_out))
    # -

    # ## Computing IIP weights

    # +
    npart_out = np.shape(we_out)[0]
    print('npart_out =', npart_out)

    nbits_out = np.shape(we_out)[1] * 64
    #nbits_out = np.shape(we_out)[1] * 31
    print('nbits_out =', nbits_out)

    wiip = np.zeros(npart, dtype=float)
    #wiip[id_msk_obs] = nbits_out / np.sum(bool_out[id_msk], axis=1)
    wiip[id_msk_obs] = (nbits_out + 1) / (np.sum(bool_out[id_msk_obs], axis=1) + 1)
    # -

    print(wiip)
    #plt.hist(wiip[id_msk_obs], np.arange(-0.5,15,0.2));

    wiip[0:11]

    # ## Reading galaxy positions and other quantities of interest

    # +
    namin_ref = new_namin

    data = ascii.read(namin_ref, format='no_header')
    RA = np.array(data["col1"])
    dec = np.array(data["col2"])
    truez = np.array(data["col3"])
    rsdz = np.array(data["col4"])
    #DESI_target = np.array(data["col5"])
    #rchmsk = np.array(data["col6"])
    #selmsk = np.array(data["col7"])
    ntile = np.array(data["col8"])
    #gcap = np.array(data["col9"])
    #ntgal = np.array(data["col10"])
    #comp = np.array(data["col11"])
    targetid = np.array(data["col10"])
    # -



    # ## Writing output

    # +
    #targetid = np.arange(1, len(RA)+1, dtype=int)

    #namout_tot = '/global/cfs/cdirs/desi/survey/catalogs/main/mocks/FAemu_preliminary/sikandar/Updated_Code_CFC/emulate_bitw_v1.1/out/Abacus2Outputs_TrainedOnData/mocks_new/Catalogs_With_Weights/Ab2_corr_clustering_cat_FAemu_m' + args.real + '_B_BGS_llen0.02_truent_cfc.fits'
    namout_tot = emu_out_dir + args.mockver + '_corr_clustering_cat_FAemu_m' + args.real + '_' + args.galcap + '_' + args.tracer + '_llen0.02_truent_cfc.fits'

    col1 = fits.Column(name='TARGETID', array=targetid[id_msk_obs], format='K')
    col2 = fits.Column(name='RA', array=RA[id_msk_obs], format='E')
    col3 = fits.Column(name='DEC', array=dec[id_msk_obs], format='E')
    col4 = fits.Column(name='TRUEZ', array=truez[id_msk_obs], format='E')
    col5 = fits.Column(name='RSDZ', array=rsdz[id_msk_obs], format='E')
    col6 = fits.Column(name='NTILE', array=ntile[id_msk_obs], format='K')
    col7 = fits.Column(name='WEIGHT_IIP', array=wiip[id_msk_obs], format='E')
    col8 = fits.Column(name='BITWEIGHT', array=we_out[id_msk_obs], dim='4', format='4K')
    #col8 = fits.Column(name='DESI_TARGET', array=DESI_target[id_msk_obs], format='K')
    #col9 = fits.Column(name='TRUEZ', array=truez[id_msk_obs], format='E')
    #col10 = fits.Column(name='RCH_MASK', array=rchmsk[id_msk_obs], format='L')
    #col11 = fits.Column(name='SEL_MASK', array=selmsk[id_msk_obs], format='L')
    #col12 = fits.Column(name='GAL_CAP', array=gcap[id_msk_obs], format='A')

    tbl = fits.BinTableHDU.from_columns([col1, col2, col3, col4, col5, col6, col7, col8])

    tbl.writeto(namout_tot, overwrite=True)

    # +
    namin_test = namout_tot

    hdul = fits.open(namin_test)
    hdul.info()
    hdul[1].header
    # -

    data_test = hdul[1].data
    #print(data_test)
    print(np.shape(data_test['RA']))
    print(np.shape(data_test['BITWEIGHT']))
    print(data_test['TARGETID'])
    print(data_test['WEIGHT_IIP'])
    print(data_test['BITWEIGHT'])

    print("DONE realisation ", args.real)


    print("Time from start to emulation:")
    print(datetime.now() - startTime)


    ## Joining 

    # python /global/homes/s/sikandar/Join_Cats.py --tracer BGS --real 0 --mocktype ab_secondgen_cosmosim --galcap B --base_output /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummitBGS/ --cat_in /global/cfs/cdirs/desi/survey/catalogs/main/mocks/FAemu_preliminary/sikandar/Updated_Code_CFC/emulate_bitw_v1.1/out/Abacus2Outputs_TrainedOnData/BGS_Cutsky/ --prog bright

    # if args.mockver == "ab_secondgen_cosmosim":
    #     add_string = "cosmosim/"
    # else:
    #     add_string = ""
    add_string = ""

    emu_out_dir = args.emulator_dir + "emulate_bitw_v1.1/out/" + args.mockver + "/"

    unred_path = args.base_output + add_string + "SecondGenMocks/AbacusSummit/mock" + str(args.real) + "/"
    emu_in = Table.read(emu_out_dir + args.mockver + '_corr_clustering_cat_FAemu_m' + args.real + '_' + args.galcap + '_' + args.tracer + '_llen0.02_truent_cfc.fits')
    unreduced = Table.read(unred_path + '/pota-' + args.prog.upper() + '_joined_unreduced_%s.fits'%args.tracer)
    emu_in.remove_columns(['RA', 'DEC', 'RSDZ', 'TRUEZ', 'NTILE'])
    cat_join = join(unreduced, emu_in, keys = 'TARGETID', join_type='left')
    outdir_data = unred_path + "ffa_full_" + args.tracer + ".fits"
    join_red = cat_join[np.unique(cat_join["TARGETID"], return_index=True)[1]]
    join_red = ct.addNS(join_red)
    join_red.filled().write(outdir_data, overwrite = True)


print("script time:")
print(datetime.now() - startTime)