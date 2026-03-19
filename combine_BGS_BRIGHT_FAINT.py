import numpy as np
from astropy.table import Table, vstack
from LSS.common_tools import write_LSS_scratchcp
from pycorr import setup_logging
import logging
import argparse
import os


setup_logging()
logger = logging.getLogger('combine_BGS_BRIGHT_FAINT')

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--ccut", help="if combining subsamples, the string that defines them (default is empty string for full catalogs)", default='')
#arguments to find input data
parser.add_argument("--basedir", help="base directory for input, note that a versioning structure is expected under this directory", default='/dvs_ro/cfs/cdirs/desi/survey/catalogs/')
parser.add_argument("--version", help="catalog version for input", default='v2')
parser.add_argument("--survey", help="e.g., Y1, DA2", default='DA2')
parser.add_argument("--verspec", help="version for redshifts", default='loa-v1')
parser.add_argument("--outdir", help="directory for output", default=os.environ['SCRATCH'])

parser.add_argument("--data", choices=['y', 'n'], help="write the data catalog?", default='y')
parser.add_argument("--random_data_ratio", choices=['unity', 'legacy'], help="how to set the random-to-data weight ratio for BRIGHT and FAINT (they must be equal). unity sets the ratio to 1; legacy is closer to the default from LSS catalogs, which takes the first ratio (in this case, BRIGHT)", default='unity')
parser.add_argument("--minr", help="minimum number for random files", default=0, type=int)
parser.add_argument("--maxr", help="maximum number for random files (plus one), 18 (0 through 17) are available (use parallel script for all)", default=18, type=int)

parser.add_argument("--par", choices=['y', 'n'], help="run different random numbers in parallel?", default='y')

args = parser.parse_args()
logger.info(f"Running with arguments: {args}")

input_dir = os.path.join(args.basedir, args.survey, 'analysis', args.verspec, 'LSScats', args.version, 'nonKP') + '/' # basedir with analysis, where we keep non-standard catalogsS
logger.info(f"Primary input directory is {input_dir}")
output_dir: str = args.outdir + '/'
logger.info(f"Output directory is {output_dir}")
input_dir_main = os.path.join(args.basedir, args.survey, 'LSS', args.verspec, 'LSScats', args.version, 'nonKP') + '/' # basedir with LSS to find the non-cut catalogs, and for fallback if the cut catalogs are not found in the analysis directory
logger.info(f"Input directory for non-cut catalogs (and fallback) is {input_dir_main}")
try_dirs = [input_dir, output_dir, input_dir_main] # order to look for input files

samples_base = ['BRIGHT', 'FAINT']
samples = [sample + args.ccut for sample in samples_base]
logger.info(f"Combining samples {samples} obtained by applying cut {args.ccut} to base samples {samples_base}")
sample_comb = '+'.join(samples_base) + args.ccut
logger.info(f"Combined sample name will be {sample_comb}")
regions = ['NGC', 'SGC']
phot_regions = ["N", "S"]


def lookup_dirs(basename: str) -> str:
    for dirname in try_dirs:
        path = os.path.join(dirname, basename)
        if os.path.isfile(path): return path
    raise FileNotFoundError(f"{basename} not found in any of the input or output directories {try_dirs}.")


def read_catalog(sample: str, reg: str, iran: int | None = None) -> Table:
    """Read the clustering catalog for a given sample and region."""
    basename = f'BGS_{sample}_{reg}' + f'_{iran}' * (iran is not None) + '_clustering.' + ('dat' if iran is None else 'ran') + '.fits'
    path = lookup_dirs(basename)
    logger.info(f"Reading clustering catalog from {path}")
    return Table.read(path)


def combine_regions(reg_cat_dict: dict[str, Table]) -> Table:
    """Combine the NGC and SGC catalogs for a given sample. Add a column to indicate the region."""
    for region, cat in reg_cat_dict.items():
        cat['REGION'] = region
    return vstack(list(reg_cat_dict.values())) # doesn't work without the list() wrapper for some reason


def get_total_weights(catalog: Table) -> np.typing.NDArray[np.float64]: return catalog['WEIGHT'] * catalog['WEIGHT_FKP']

def select_reg(catalog: Table, reg: str) -> Table: return catalog[catalog['REGION'] == reg]

def select_phot(catalog: Table, phot_region: str) -> Table: return catalog[catalog['PHOTSYS'] == phot_region]


# preparation for NX computation
nz_data = {}
comp_ntile_factors = {}
for (sample, sample_base) in zip(samples, samples_base):
    nz_data[sample] = {}
    comp_ntile_factors[sample] = {}
    for reg in regions:
        nz_name = lookup_dirs(f'BGS_{sample}_{reg}_nz.txt')
        logger.info(f"Reading nz data for sample {sample} and region {reg} from {nz_name}")
        nz_data[sample][reg] = np.loadtxt(nz_name).T
        logger.info(f"Obtaining NTILE data for sample {sample_base} and region {reg}")
        base_data = Table.read(input_dir_main + f'BGS_{sample_base}_{reg}_clustering.dat.fits') # the base sample catalog should be in the main input dir
        comp_ntl = np.bincount(base_data['NTILE']-1) / np.bincount(base_data['NTILE']-1, weights=base_data['WEIGHT_COMP']) # inverse of the mean completeness weight in data for each NTILE value (note that it is shifted down by 1)
        fttl = np.bincount(base_data['NTILE']-1, weights=base_data['FRAC_TLOBS_TILES']) / np.bincount(base_data['NTILE']-1) # mean of FRAC_TLOBS_TILES for each (positive integer) NTILE in data (randoms might be more correct for this; note that NTILE is shifted down by 1)
        comp_ntile_factors[sample][reg] = comp_ntl * fttl # indexed by NTILE-1
        del base_data # free memory


def update_NX_FKP_weight(reg_cat_dict: dict[str, Table], P0: float = 7000):
    """
    Compute the combined NX for each sample by summing contribution for each sample, separately for each region, and update the FKP weight in-place in the input catalog.
    P0=7000 is the default and probably better for full BGS (BRIGHT); P0=8400 is probably more appropriate for -21.35, but it should not make a lot of difference
    """
    for reg in regions:
        fd = reg_cat_dict[reg] # cut data catalog

        NX = np.zeros(len(fd))
        for sample in samples: # accumulate NX contribution from each sample we combine
            nzf = nz_data[sample][reg]  # nz data for this sample and region
            zmin, zmax = nzf[1][0], nzf[2][-1]
            bs = nzf[2][0]-nzf[1][0] # z step
            nzd = nzf[3]  # column with nbar values
            zl = fd['Z']
            nl = np.zeros(len(zl)) # n(z) values
            zind = ((zl - zmin) / bs).astype(int)
            valid = (zl > zmin) & (zl < zmax) & (zind >= 0) & (zind < len(nzd))
            nl[valid] = nzd[zind[valid]]

            NX += nl * comp_ntile_factors[sample][reg][fd['NTILE']-1] # apply the completeness correction factor for the corresponding NTILE value, which is indexed by NTILE-1
        
        reg_cat_dict[reg]['NX'] = NX # update the NX column for this region
        reg_cat_dict[reg]['WEIGHT_FKP'] = 1 / (1 + P0 * NX) # update the FKP weight for this region


# deal with data first, which is more straightforward, and then randoms, which is more complicated
logger.info(f"Starting data processing")
data = {sample: {reg: read_catalog(sample, reg) for reg in regions} for sample in samples}

for sample in samples:
    logger.info(f"Updating NX and FKP weight for {sample} data")
    update_NX_FKP_weight(data[sample])

os.makedirs(output_dir, exist_ok=True) # make sure the output directory exists
# now data combination is straightforward: just concatenation in each region. can write them quickly before doing the combination for randoms, which is more complicated because we want to preserve the data-to-random ratio in each region and photometric region, and also preserve the sky number density for randoms
if args.data == 'y':
    for reg in regions:
        logger.info(f"Combining {samples} data for region {reg}")
        this_data_comb = vstack([data[sample][reg] for sample in samples])
        path = os.path.join(output_dir, f'BGS_{sample_comb}_{reg}_clustering.dat.fits')
        logger.info(f"Writing {sample_comb} data for region {reg} to {path}")
        write_LSS_scratchcp(this_data_comb, path, logger=logger)
    del this_data_comb # free memory

# compute the number and sum of weights for data in each photometric region, which will be used to balance the randoms in each photometric region
n_data = {}
wsum_data = {}
for sample in samples:
    logger.info(f"Combining regions {regions} for {sample} data")
    data[sample] = combine_regions(data[sample]) # combine NGC and SGC for each sample, which is more convenient for balancing weights in photometric regions (N and S) for data
for phot_region in phot_regions:
    logger.info(f"Computing number of galaxies and sum of weights in photometric region {phot_region}")
    data_in_region = [select_phot(data[sample], phot_region) for sample in samples]
    n_data[phot_region] = np.array([len(this_data_in_region) for this_data_in_region in data_in_region])
    wsum_data[phot_region] = [get_total_weights(this_data_in_region).sum() for this_data_in_region in data_in_region]
del data_in_region, data # no longer need data
logger.info(f"Finished with data processing")

def process_random(iran: int):
    """
    Deal with randoms, one file at a time.
    This is more complicated than data because we want to preserve the data-to-random ratio in each region and photometric region, and also preserve the sky number density for randoms.
    """
    logger.info(f"Starting processing for random number {iran}")
    random = {sample: {reg: read_catalog(sample, reg, iran=iran) for reg in regions} for sample in samples}

    for sample in samples:
        logger.info(f"Updating NX and FKP weight for {sample} randoms for random number {iran}")
        update_NX_FKP_weight(random[sample])
        logger.info(f"Combining regions {regions} for {sample} randoms for random number {iran}")
        random[sample] = combine_regions(random[sample]) # combine NGC and SGC for each sample, which is more convenient for balancing weights in photometric regions (N and S)

    random[sample_comb] = {}
    np.random.seed(42+iran) # reproducible but different for each random number. otherwise, the BRIGHT and FAINT randoms will have the same TARGETID across the random numbers

    # combination
    logger.info(f"Sorting randoms by TARGETID for random number {iran}")
    for sample in samples: random[sample].sort('TARGETID') # sort randoms by TARGETID
    random_comb = []
    for phot_region in phot_regions:
        logger.info(f"Selecting randoms in photometric region {phot_region} for random number {iran}")
        # for randoms, we want to preserve the number, so we draw from two sets following the data proportion
        randoms_in_region = [select_phot(random[sample], phot_region) for sample in samples]
        for i_sample in range(1, len(samples)): assert np.array_equal(randoms_in_region[0]['TARGETID'], randoms_in_region[i_sample]['TARGETID']) # must be equal (after sorting above)
        n_randoms_tot = len(randoms_in_region[0]) # now doesn't matter which sample for randoms
        n_random_goals = np.rint(n_randoms_tot * n_data[phot_region] / n_data[phot_region].sum()).astype(int) # number of randoms to draw for each sample, following the data proportion; round to integer
        logger.info(f"Total number of randoms in photometric region {phot_region} for random number {iran} is {n_randoms_tot}, the number of randoms to draw for each sample is {n_random_goals}, and the data numbers are {n_data[phot_region]}")
        all_random_indices = np.random.permutation(n_randoms_tot) # shuffle the indices of randoms
        n_random_split = np.cumsum(n_random_goals)[:-1] # get the indices to split the shuffled randoms
        random_indices = np.split(all_random_indices, n_random_split) # split the shuffled indices according to the number of randoms to draw for each sample
        random_to_data_ratios_orig = [get_total_weights(these_randoms_in_region).sum() / wsum_data[phot_region][i] for i, these_randoms_in_region in enumerate(randoms_in_region)] # compute the random-to-data weight ratio for each sample in this photometric region before subsampling
        randoms_in_region = [these_randoms[these_indices] for these_randoms, these_indices in zip(randoms_in_region, random_indices)] # subsample/select the randoms for each sample
        random_to_data_ratios = [get_total_weights(these_randoms_in_region).sum() / wsum_data[phot_region][i] for i, these_randoms_in_region in enumerate(randoms_in_region)] # compute the random-to-data weight ratio for each sample in this photometric region after subsampling but before reweighting
        target_ratio = 1 if args.random_data_ratio == 'unity' else random_to_data_ratios_orig[0] # the target random-to-data weight ratio for all samples in this photometric region; if 'unity', set to 1; if 'legacy', set to the original ratio for the first sample (which is BRIGHT)
        logger.info(f"Reweighting in photometric region {phot_region} for random number {iran}. Original random-to-data weight ratios are {random_to_data_ratios_orig}, after subsampling they are {random_to_data_ratios}; target random-to-data weight ratio is {target_ratio}")
        for i in range(len(samples)): randoms_in_region[i]['WEIGHT'] *= target_ratio / random_to_data_ratios[i] # make random-to-data ratio target_ratio for BRIGHT and FAINT parts in each region
        logger.info(f"Stacking in photometric region {phot_region} for random number {iran}")
        random_comb.append(vstack(randoms_in_region))
        logger.info(f"Finished with photometric region {phot_region} for random number {iran}")
    del randoms_in_region, all_random_indices, n_random_split, random_indices # free memory
    logger.info(f"Stacking photometric regions {phot_regions} for random number {iran}")
    random_comb = vstack(random_comb)

    # should balance weights for z<0.5 only?
    # ok/good to manage imaging systematics separately

    for reg in regions:
        logger.info(f"Selecting {sample_comb} randoms for region {reg} and random number {iran}")
        this_random_comb = select_reg(random_comb, reg)
        path = os.path.join(output_dir, f'BGS_{sample_comb}_{reg}_{iran}_clustering.ran.fits')
        logger.info(f"Writing {sample_comb} randoms for region {reg} and random number {iran} to {path}")
        this_random_comb.remove_column('REGION') # remove the REGION column before writing, to be consistent with the original files
        write_LSS_scratchcp(this_random_comb, path, logger=logger)
    logger.info(f"Finished processing for random number {iran}")

# deal with randoms
random_inds = np.arange(args.minr, args.maxr)
logger.info(f"Starting random processing for indices {args.minr} to {args.maxr-1} " + ("in parallel" if args.par == 'y' else "sequentially"))
if args.par == 'y':
    from multiprocessing import Pool
    with Pool() as pool:
        res = pool.map(process_random, random_inds)
else:
    for ii in random_inds:#range(rm,rx):
        process_random(ii)

logger.info("Finished with random processing")
