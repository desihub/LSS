from matplotlib import pyplot as plt
from astropy.io import (fits, ascii)
import numpy as np
from astropy.table import Table, join, vstack, unique
import LSS.common_tools as ct
from datetime import datetime
import argparse
startTime = datetime.now()

parser = argparse.ArgumentParser()
parser.add_argument("--emulator_dir", help = "base directory of ffa emulator")
parser.add_argument("--tracer")
parser.add_argument("--prog")
parser.add_argument("--real")
parser.add_argument("--mocktype", default = "ab_secondgen_cosmosim")
parser.add_argument("--galcap")
parser.add_argument("--cat_in", help = "base directory from where to read the catalog with weights")
parser.add_argument("--base_output")
args = parser.parse_args()

# python /global/homes/s/sikandar/Join_Cats.py --tracer BGS --real 0 --mocktype ab_secondgen_cosmosim --galcap B --base_output /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummitBGS/ --cat_in /global/cfs/cdirs/desi/survey/catalogs/main/mocks/FAemu_preliminary/sikandar/Updated_Code_CFC/emulate_bitw_v1.1/out/Abacus2Outputs_TrainedOnData/BGS_Cutsky/ --prog bright

if args.mocktype == "ab_secondgen_cosmosim":
    add_string = "cosmosim/"
else:
    add_string = ""

unred_path = args.base_output + add_string + "SecondGenMocks/AbacusSummit/mock" + str(args.real) + "/"
emu_in = Table.read(args.cat_in + args.mocktype + '_corr_clustering_cat_FAemu_m' + args.real + '_' + args.galcap + '_' + args.tracer + '_llen0.02_truent_cfc.fits')
unreduced = Table.read(unred_path + '/pota-' + args.prog.upper() + '_joined_unreduced_%s.fits'%args.tracer)
emu_in.remove_columns(['RA', 'DEC', 'RSDZ', 'TRUEZ', 'NTILE'])
cat_join = join(unreduced, emu_in, keys = 'TARGETID', join_type='left')
outdir_data = unred_path + "ffa_full_" + args.tracer + ".fits"
join_red = cat_join[np.unique(cat_join["TARGETID"], return_index=True)[1]]
join_red = ct.addNS(join_red)
join_red.filled().write(outdir_data, overwrite = True)


print("script time:")
print(datetime.now() - startTime)