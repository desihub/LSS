# Run this script from LSS scripts directory
import argparse
import subprocess
from datetime import datetime
startTime = datetime.now()
#python FFA_preprocessing.py --mockver ab_secondgen_cosmosim --realization 0 --prog bright --base_output /global/cfs/cdirs/desi/cosmosim/SecondGenMocks/AbacusSummit/CutSky/BGS/v0/z0.200/v0.1outs/ --tracer BGS --emulator_dir /global/cfs/cdirs/desi/survey/catalogs/main/mocks/FAemu_preliminary/sikandar/Updated_Code_CFC/

parser = argparse.ArgumentParser()
parser.add_argument("--mockver", help="type of mock to use",default=None)
parser.add_argument("--mockpath", help="Location of mock file(s)",default='/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/AbacusSummit/CutSky/')
parser.add_argument("--mockfile", help="formattable name of mock file(s). e.g. cutsky_{TYPE}_{Z}_AbacusSummit_base_c000_ph{PH}.fits. TYPE will be replaced with tracer type. PH will be replaced with realization number for simulation of mock.",default='cutsky_{TYPE}_{Z}_AbacusSummit_base_c000_ph{PH}.fits')
parser.add_argument("--real", help="number for the realization",default=0,type=int)
parser.add_argument("--prog", help="dark or bright",default='dark')
parser.add_argument("--base_output", help="base directory for output",default='/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/')
parser.add_argument("--tracer", help="which tracer to do")
parser.add_argument("--emulator_dir", help="base directory of ffa emulator")
args = parser.parse_args()

if args.mocktype == "ab_secondgen_cosmosim":
    add_string = "cosmosim/"
else:
    add_string = ""


getpota_indir = args.base_output + add_string 
getpota_outdir = args.base_output + add_string
NTILE_assign_indir = args.base_output + add_string + "SecondGenMocks/AbacusSummit/mock" + str(args.real) + "/"
DESIwemu_indir = args.base_output + add_string + "SecondGenMocks/AbacusSummit/" 
DESIwemu_outdir = args.emulator_dir + "fof_v1.0/in/"

cmd_string1 = "python /global/homes/s/sikandar/PerlmutterLSS/LSS/scripts/mock_tools/prepare_mocks_Y1.py --rbandcut 19.5 --mockver %s --prog %s --downsampling n --overwrite 1 --realmin %s --realmax %s --base_output %s" %(args.mockver, args.prog.lower(), args.real, args.real + 1, args.base_output)
cmd_string2 = "python /global/homes/s/sikandar/PerlmutterLSS/LSS/scripts/getpotaY1_mock.py --base_output %s --prog %s --realization %s --base_input %s"%(getpota_indir, args.prog.upper(), args.real, getpota_outdir)
cmd_string3 = "python /global/homes/s/sikandar/PerlmutterLSS/LSS/scripts/mock_tools/NTILE_assign.py --indir %s --tracer %s --tileloc_add y --prog %s"%(NTILE_assign_indir, args.tracer, args.prog.upper())
cmd_string4 = "python /global/cfs/cdirs/desi/survey/catalogs/main/mocks/FAemu_preliminary/sikandar/Updated_Code_CFC/DESIAbacuswemu.py --mocknumber %s --tracer %s --basedir %s --prog %s --overwrite y --outdir %s --mocktype %s"%(args.real, args.tracer, DESIwemu_indir, args.prog.upper(),DESIwemu_outdir, args.mockver)
print(cmd_string1)
print(cmd_string2)
print(cmd_string3)
print(cmd_string4)
subprocess.run(cmd_string1, shell = True)
print("Done with prepare_mocks_y1")
subprocess.run(cmd_string2, shell = True)
print("Done with getpotaY1_mock")
subprocess.run(cmd_string3, shell = True)
print("done with NTILE_assign")
subprocess.run(cmd_string4, shell = True)
print("all done")
print("script time:")
print(datetime.now() - startTime)