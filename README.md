# LSS
Codes used to create LSS catalogs and randoms

Versions developed for different specific settings are being created. These are in separate folders in the py/LSS directory. Common tools go in the LSS directory.

Scripts meant to act like executables are in the bin directory. Currently, one should enter the bin directory to run everything (and then the paths work). The default mode will produce output in your CSCRATCH directory on NERSC. This code only works on NERSC.

Some examples:

Make sure to be in the DESI environment, e.g., run

source /project/projectdirs/desi/software/desi_environment.sh master

before trying anything else (master should soon switch to main).

python gatherSV_zinfo_alltiles.py --type ELG --release blanc #this gathers all of the redshift info for SV1 ELG targets in the blanc release (daily supported as release, LRG, QSO, BGS_ANY are additional supported types as should be anything in the SV1 DESIMASK)

python mkCat_singletile.py --type ELG --tile 80623 --night deep #this creates clustering catalogs for ELGs (LRG, QSO, BGS_ANY all work) on tile 80623 (see https://desi.lbl.gov/trac/wiki/SurveyValidation/SV1#ObservedTiles for the initial list) using the deep coadd (other options are a particular night or all, though all is not recommended) from the blanc release redshift data (use --release cascades to update to the more recent data). There are other options that use 'n' or 'y' to toggle whether a particular stage in the process gets run.
