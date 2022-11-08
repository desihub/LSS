# LSS

Introduction
------------

This package contains codes to create LSS catalogs and randoms.

Installation
------------

First, git clone this repo:

If you want any change to the code to take place immediately, either:

1.  Add the "bin" directory to your
    ``$PATH`` environment variable and add the "py" directory to your
    ``$PYTHONPATH`` environment variable.

2.  Install (and uninstall) the current git checkout:

    $>  python setup.py develop --user

    $>  python setup.py develop --user --uninstall

You can also install a fixed version of the package:

    $>  python setup.py install --user

This will put the main scripts bin/mkCat_* in your $PATH (typically in $HOME/.local/bin).

This code only works on NERSC, in the DESI environment. Therefore, make sure to be in the DESI environment, e.g., run

    $>  source /global/common/software/desi/desi_environment.sh main

before trying anything else (master should soon switch to main). Some code also requires the cosmodesi environment to be loaded (e.g., for 2pt functions and/or reconstruction, regressis tools, etc.) After sourcing the above, source it

    $>  source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main


Structure
---------

Versions developed for different specific settings are being created. These are in separate folders in the py/LSS directory. Common tools go in the LSS directory.

Scripts meant to act like executables are in the bin directory. The default mode will produce output in your CSCRATCH directory on NERSC.

Examples
--------

python bin/mkCat_SV3_simp.py --type QSO #This re-creates the version 1 catalogs writing them into your scratch by default. One should be able to test edits to the catalogs using this (starting from some pre-defined inputs). Available types are QSO, BGS_ANY, BGS_BRIGHT, ELG, ELG_HIP, LRG, MWS_ANY

python bin/gatherSV_zinfo_alltiles.py --type ELG --release blanc #this gathers all of the redshift info for SV1 ELG targets in the blanc release (daily supported as release, LRG, QSO, BGS_ANY are additional supported types as should be anything in the SV1 DESIMASK)

python bin/mkCat_singletile.py --type ELG --tile 80623 --night deep #this creates clustering catalogs for ELGs (LRG, QSO, BGS_ANY all work) on tile 80623 (see https://desi.lbl.gov/trac/wiki/SurveyValidation/SV1#ObservedTiles for the initial list) using the deep coadd (other options are a particular night or all, though all is not recommended) from the blanc release redshift data (use --release cascades to update to the more recent data). There are other options that use 'n' or 'y' to toggle whether a particular stage in the process gets run.
