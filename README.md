# LSS
Codes used to create LSS catalogs and randoms

Versions developed for different specific settings are being created. These are in separate folders in the py/LSS directory. Common tools go in the LSS directory.

Scripts meant to act like executables are in the bin directory. Currently, one should enter the bin directory to run everything (and then the paths work). Some examples:

python gatherSV_zinfo_alltiles.py --type ELG --release blanc #this gathers all of the redshift info for SV1 ELG targets in the blanc release (daily supported as release, LRG, QSO, BGS_ANY are additional supported types as should be anything in the SV1 DESIMASK)
