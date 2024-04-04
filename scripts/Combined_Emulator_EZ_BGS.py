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
 
#print('\n\n my edits are showing up\n\n')
startTime = datetime.now()
#python /global/homes/s/sikandar/Combined_Emulator.py -mockver ab_secondgen_cosmosim --real 1 --prog bright --base_output /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummitBGS/ --tracer BGS --emulator_dir /global/cfs/cdirs/desi/survey/catalogs/main/mocks/FAemu_preliminary/sikandar/Updated_Code_CFC/ --galcap B

parser = argparse.ArgumentParser()
parser.add_argument("--mockver", help="type of mock to use",default=None)
parser.add_argument("--mockpath", help="Location of mock file(s)",default='/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/AbacusSummit/CutSky/')
parser.add_argument("--mockfile", help="formattable name of mock file(s). e.g. cutsky_{TYPE}_{Z}_AbacusSummit_base_c000_ph{PH}.fits. TYPE will be replaced with tracer type. PH will be replaced with realization number for simulation of mock.",default='cutsky_{TYPE}_{Z}_AbacusSummit_base_c000_ph{PH}.fits')
parser.add_argument("--real", help="number for the realization",default=0)
parser.add_argument("--prog", help="dark or bright",default='bright')
parser.add_argument("--base_output", help="base directory for output", default = "/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/EZmock/FFA_BGS/FFA_temp/")#default='/pscratch/sd/s/sikandar/Y1/mocks/')
parser.add_argument("--tracer", help="which tracer to do", nargs = '+')
parser.add_argument("--emulator_dir", help="base directory of ffa emulator")
parser.add_argument("--galcap")
parser.add_argument("--prep", help = "do preprocessing if y", default='n')
parser.add_argument("--emulate", help = "do emulation if y", default='n')
parser.add_argument("--clear_files", default = 'n')
parser.add_argument("--copy_files", default = 'n')
parser.add_argument("--foflinklen", default = '0.02')
parser.add_argument("--emubeta0", default = '0.7')
parser.add_argument("--emubeta1", default = '0.')
parser.add_argument("--adj_qref", default = 'n')
parser.add_argument("--downsample_bgs", default = 'n')
#parser.add_argument("--deco", help = "anticorrelation parameter for emulator", default = 0.1)
args = parser.parse_args()
tracer_arr = np.array(args.tracer)
tracer_string = " ".join(tracer_arr)
# if args.mockver == "ab_secondgen_cosmosim":
#     add_string = "cosmosim/"
# else:
#     add_string = ""
add_string = ""
path2LSS='/global/homes/j/jlasker/.local/desicode/'
if (args.mockver == "EZmock") and ("BGS" not in args.tracer):
    prep_script = f"{path2LSS}/LSS/scripts/mock_tools/prepare_mocks_Y1EZ.py"
    subdir = "EZmock"
elif (args.mockver == "EZmock") and ("BGS" in args.tracer):
    prep_script = f"{path2LSS}/LSS/scripts/mock_tools/prepare_mocks_Y1.py"
    subdir = "EZmock"
elif (args.mockver.lower() == 'glam'):
    prep_script = f"{path2LSS}/LSS/scripts/mock_tools/prepare_mocks_Y1.py"
    subdir = 'GLAM'
else:
    prep_script = f"{path2LSS}/LSS/scripts/mock_tools/prepare_mocks_Y1.py"
    subdir = args.mockver #"AbacusSummit"


getpota_indir = args.base_output + add_string 
getpota_outdir = args.base_output + add_string
#getpota_tiledir = os.path.join(os.environ['SCRATCH'], 'rantiles','mock'+args.real)
getpota_tiledir = os.path.join("/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/EZmock/FFA/random_tiles/", 'mock'+args.real)
getpota_tiledir = getpota_tiledir + "/"
if not os.path.exists(getpota_tiledir):
    os.makedirs(getpota_tiledir)
    print('made ' + getpota_tiledir)
#if args.mockver.lower() == 'glam':
#    NTILE_assign_indir = args.base_output + add_string + "/mock" + str(args.real) + "/"
#    DESIwemu_indir = args.base_output + add_string + "/" 
#    DESIwemu_outdir = args.emulator_dir + "fof_v1.0/in/"
#else:
NTILE_assign_indir = args.base_output + add_string + "SecondGenMocks/" + subdir + "/mock" + str(args.real) + "/"
DESIwemu_indir = args.base_output + add_string + "SecondGenMocks/" + subdir +"/" 
DESIwemu_outdir = args.emulator_dir + "fof_v1.0/in/"
#cmd_string1 = "python /global/homes/s/sikandar/PerlmutterLSS/LSS/scripts/mock_tools/prepare_mocks_Y1.py --rbandcut 19.5 --mockver %s --prog %s --downsampling n --overwrite 1 --realmin %s --realmax %s --base_output %s" %(args.mockver, args.prog.lower(), args.real, int(args.real) + 1, args.base_output)
if args.mockver == "EZmock":
    cmd_string1 = "python %s --mockver %s --realmin %s --realmax %s --prog %s --base_output %s --downsampling n"%(prep_script, args.mockver, args.real, int(args.real) + 1, args.prog.lower(), "/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/EZmock/FFA_BGS/FFA_temp/SecondGenMocks/EZmock/forFA/")
elif args.mockver.lower() == 'glam':
    getpota_outdir = os.path.join(args.base_output, 'SecondGenMocks', 'GLAM')
    cmd_string1 = "python %s --mockpath %s --mockver %s --prog %s --downsampling n --overwrite 1 --realmin %s --realmax %s --base_output %s --mockver glam --tracer=%s " %(prep_script, args.mockpath, args.mockver, args.prog.lower(), args.real, int(args.real) + 1, getpota_outdir, tracer_string)
else:
    cmd_string1 = "python %s --rbandcut 19.5 --mockver %s --prog %s --downsampling n --overwrite 1 --realmin %s --realmax %s --base_output %s" %(prep_script, args.mockver, args.prog.lower(), args.real, int(args.real) + 1, args.base_output)
if args.mockver.lower() == 'glam':
    cmd_string2 = "python %s/LSS/scripts/getpotaY1_mock.py --base_output %s --prog %s --realization %s --base_input %s --mock %s --tile-temp-dir %s --tracer=%s"%(path2LSS,getpota_indir, args.prog.upper(), args.real, getpota_outdir, args.mockver,getpota_tiledir, tracer_string)
else:
    cmd_string2 = "python %s/LSS/scripts/getpotaY1_mock.py --base_output %s --prog %s --realization %s --base_input %s --mock %s --tile-temp-dir %s"%(path2LSS,getpota_indir, args.prog.upper(), args.real, getpota_outdir, args.mockver,getpota_tiledir)
cmd_string3 = "python %s/LSS/scripts/mock_tools/NTILE_assign.py --indir %s --tracer %s --tileloc_add y --prog %s"%(path2LSS,NTILE_assign_indir, tracer_string, args.prog.upper())

cmd_string4 = {}
for i in range(len(tracer_arr)):
    cmd_string4[tracer_arr[i]] = "python %sDESI_wemu.py --mocknumber %s --tracer %s --basedir %s --prog %s --overwrite y --outdir %s --mocktype %s"%(args.emulator_dir, args.real, tracer_arr[i], DESIwemu_indir, args.prog.upper(),DESIwemu_outdir, args.mockver)
#     if i==0:
#         cmd_string4[tracer_arr[i]] = "python %sDESIAbacuswemu.py --mocknumber %s --tracer %s --basedir %s --prog %s --overwrite y --outdir %s --mocktype %s"%(args.emulator_dir, args.real, tracer_arr[i], DESIwemu_indir, args.prog.upper(),DESIwemu_outdir, args.mockver)
#     else: 
#         cmd_string4[tracer_arr[i]] = "python %sDESIAbacuswemu.py --mocknumber %s --tracer %s --basedir %s --prog %s --overwrite y --outdir %s --mocktype %s skip_specy"%(args.emulator_dir, args.real, tracer_arr[i], DESIwemu_indir, args.prog.upper(),DESIwemu_outdir, args.mockver)
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
    for i in args.tracer:
        subprocess.run(cmd_string4[i], shell = True)
    print("all done")
else:
    print("skipped preprocessing")
print("FFA preprocessing time:")
print(datetime.now() - startTime)

# # EmulatorRun

for tr_i in tracer_arr:   
    tardict = {4: 'QSO', 1: 'LRG', 34: 'ELG', 60: 'BGS'}
    emulator_basedir = args.emulator_dir
    fofin_basedir = emulator_basedir + "fof_v1.0/in/"
    fofin_file = tr_i + "_" + args.mockver + "_forFAemu_m" + args.real + ".dat" 
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

    new_namout_part = fof_prefix + 'out/' + args.mockver + '/' + args.mockver + '_forFAemu_m' + args.real +'_' + new_galcap + '_' + new_target + '_truent_part_llen' + args.foflinklen + '.dat'
    new_namout_halo = fof_prefix + 'out/' + args.mockver + '/' + args.mockver + '_forFAemu_m' + args.real +'_' + new_galcap + '_' + new_target + '_truent_halo_llen' + args.foflinklen + '.dat'

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
        elif lines[i].startswith("linklen"):
            lines[i] = f"linklen    {args.foflinklen + 'd0'}\n"
        elif lines[i].startswith("mocknum"):
            lines[i] = f"mocknum     {args.real}\n"

    file_path_new_fof = emulator_basedir + "fof_v1.0/INI_fof_lightcone_ang_%s_m%s_ll_%s.txt"%(args.mockver, args.real, args.foflinklen)
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
    emu_pfxselprop = emulator_basedir + "emulate_bitw_v1.1/out/" + args.mockver + "/" + "FAemu_m" + args.real + "_" + args.galcap + "_" + new_target + "_llen" + args.foflinklen + "_beta" + args.emubeta1 + "_truent_cfc_"
    emu_naminqref = emulator_basedir + "emulate_bitw_v1.1/out/TrainedEmulator_Data/FAemu_data_" + args.galcap + "_" + new_target + "_llen" + "0.02" + "_truent_cfc_antico_deco0.1_beta0_indselfrac_cfc_qm.dat"
    emu_pfxout = emulator_basedir + "emulate_bitw_v1.1/out/" + args.mockver + "/FAemu_m" + args.real + "_" + args.galcap + "_" + new_target + "_llen" + args.foflinklen + "_beta" + args.emubeta1 +  "_truent_cfc_nbits217_"
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
        elif lines[i].startswith("beta0"):
            lines[i] = f"beta      {args.emubeta0 + 'd0'}\n"
        elif lines[i].startswith("beta1"):
            lines[i] = f"beta      {args.emubeta1 + 'd0'}\n"
        elif lines[i].startswith("adj_qref_lg"):
            if args.adj_qref == 'y':
                lines[i] = f"adj_qref_lg        {'T'}\n"
            else:
                lines[i] = f"adj_qref_lg        {'F'}\n"
        elif lines[i].startswith("mocknum"):
            lines[i] = f"mocknum     {args.real}\n"
        


    # Write the modified content back to the same file
    file_path_new_emu = emulator_basedir + "/emulate_bitw_v1.1/INI_emulate_bitw_lightcone_ntile_%s_%s_ll_%s_beta_%s.txt"%(args.mockver, args.real, args.foflinklen, args.emubeta1)
    with open(file_path_new_emu, "w") as file:
        file.writelines(lines)


    if args.emulate == 'n':
        print("skipped emulation")

        
    elif args.emulate == 'y':
        print("starting fof")
        cmd_string5 = "sh " + emulator_basedir + "fof_v1.0/COMPILE.sh"
        cmd_string6 = emulator_basedir + "fof_v1.0/fof_lightcone_ang " + file_path_new_fof
        outputfn = emulator_basedir + "fof_v1.0/log/" + args.mockver + "_fof_" + new_target + "_llen"+ args.foflinklen +"_" + args.galcap + "_m" + args.real + "_datatrain.txt"
        errorfn = emulator_basedir + "fof_v1.0/log/" + "ERR_" + args.mockver + "_fof_" + new_target + "_llen" + args.foflinklen + "_" + args.galcap + "_m" + args.real + "_datatrain.txt"
        print(cmd_string5)
        print(cmd_string6)
        with open(outputfn, "w+") as output, open(errorfn, "w+") as errorf:
             #subprocess.run(cmd_string5, shell=True)
             subprocess.run(cmd_string6, shell=True, stdout=output, stderr = errorf)
        print("done fof")


        print("starting emulation")
        cmd_string7 = "sh " + emulator_basedir + "emulate_bitw_v1.1/COMPILE.sh"
        cmd_string8 = emulator_basedir + "emulate_bitw_v1.1/emulate_bitw_lightcone_ntile " + file_path_new_emu
        outputfnemu = emulator_basedir + "emulate_bitw_v1.1/log/" + args.mockver + "_emu_" + new_target + "_llen"+ args.foflinklen +"_" + args.galcap + "_m" + args.real + "_datatrain.txt"
        errorfnemu = emulator_basedir + "emulate_bitw_v1.1/log/" + "ERR_" + args.mockver + "_emu_" + new_target + "_llen"+ args.foflinklen +"_" + args.galcap + "_m" + args.real + "_datatrain.txt"
        print(cmd_string7)
        print(cmd_string8)
        with open(outputfnemu, "w+") as output, open(errorfnemu, "w+") as errorf:
             #subprocess.run(cmd_string7, shell=True)
             subprocess.run(cmd_string8, shell=True, stdout=output, stderr = errorf)
        print("done emulation")





        # WEIGHTS
        nwe = 9
        namin_fbw = emulator_basedir + "emulate_bitw_v1.1/out/" + args.mockver + "/FAemu_m" + args.real + "_" + args.galcap + "_" + new_target + "_llen"+ args.foflinklen + "_beta" + args.emubeta1 + "_truent_cfc_nbits217_" + "wemu_unformatted.dat"

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

        namout_bool = emucat_outdir + args.mockver + '_' + new_target + '_' + args.galcap + '_bool_tot_datatrain_m' + args.real + '.bin'

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
        namout_tot = emu_out_dir + args.mockver + '_corr_clustering_cat_FAemu_m' + args.real + '_' + args.galcap + '_' + new_target + '_llen'+ args.foflinklen + "_beta" + args.emubeta1 + '_truent_cfc.fits'

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
        if args.mockver == "EZmock":
            pathstr = "SecondGenMocks/EZmock/"
        elif args.mockver == "ab_secondgen_cosmosim":
            pathstr = "SecondGenMocks/AbacusSummit/"
        elif args.mockver.lower() == "glam":
            pathstr = "SecondGenMocks/GLAM/"

        if new_target == "ELG":
            target_write_name = "ELG_LOP"
        else:
            target_write_name = new_target
        unred_path = args.base_output + add_string + pathstr + "mock" + str(args.real) + "/"
        emu_in = Table.read(emu_out_dir + args.mockver + '_corr_clustering_cat_FAemu_m' + args.real + '_' + args.galcap + '_' + new_target + '_llen'+ args.foflinklen + "_beta" + args.emubeta1 + '_truent_cfc.fits')
        unreduced = Table.read(unred_path + '/pota-' + args.prog.upper() + '_joined_unreduced_%s.fits'%new_target)
        emu_in.remove_columns(['RA', 'DEC', 'RSDZ', 'TRUEZ', 'NTILE'])
        cat_join = join(unreduced, emu_in, keys = 'TARGETID', join_type='left')
        #outdir_data = unred_path + "ffa_full_" + target_write_name + "_ll_" + args.foflinklen + "_beta" + args.emubeta1 + ".fits"
        outdir_data = unred_path + "ffa_full_" + target_write_name + ".fits"
        join_red = cat_join[np.unique(cat_join["TARGETID"], return_index=True)[1]]
        join_red = ct.addNS(join_red)

        if args.downsample_bgs == 'n':
            join_red.filled().write(outdir_data, overwrite = True)
        elif args.downsample_bgs == 'y':
            print("Downsampling before writing")
            bgsfull = join_red.filled()
            bgsbright = Table.read("/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/v1/BGS_BRIGHT-21.5_full_HPmapcut.dat.fits")
            sel_fraction = np.random.choice(bgsfull, size = int(1.5*len(bgsbright)), replace = False) 
            Table(sel_fraction).write(outdir_data, overwrite = True)

        print("Done with the emulation process for ", new_target)

if args.clear_files == 'y':
    cmd_rm1 = "rm %spota-DARK_withntile_LRG.fits"%(NTILE_assign_indir)
    cmd_rm2 = "rm %spota-DARK_withntile_QSO.fits"%(NTILE_assign_indir)
    cmd_rm3 = "rm %spota-DARK_withntile_ELG.fits"%(NTILE_assign_indir)
    cmd_rm4 = "rm %spota-DARK.fits"%(NTILE_assign_indir)
    print(cmd_rm1)
    print(cmd_rm2)
    print(cmd_rm3)
    print(cmd_rm4)
    subprocess.run(cmd_rm1, shell = True)
    subprocess.run(cmd_rm2, shell = True)
    subprocess.run(cmd_rm3, shell = True)
    subprocess.run(cmd_rm4, shell = True)
    print("Done removing")

if args.copy_files == 'y':
    cmd_cp1 = "cp -r -f /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/EZmock/FFA_BGS/FFA_temp/SecondGenMocks/EZmock/mock%s /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/EZmock/FFA_BGS/"%args.real
    cmd_cp2 = "cp -r -f /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/EZmock/FFA_BGS/FFA_temp/SecondGenMocks/EZmock/forFA/forFA%s.fits /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/EZmock/FFA/forFA/"%args.real
    subprocess.run(cmd_cp1, shell = True)
    subprocess.run(cmd_cp2, shell = True)
    print("Done copying")




print("script time:")
print(datetime.now() - startTime)
