LSS
===

|Actions Status| |Coveralls Status| |Documentation Status|

.. |Actions Status| image:: https://github.com/desihub/LSS/workflows/CI/badge.svg
    :target: https://github.com/desihub/LSS/actions
    :alt: GitHub Actions CI Status

.. |Coveralls Status| image:: https://coveralls.io/repos/desihub/LSS/badge.svg
    :target: https://coveralls.io/github/desihub/LSS
    :alt: Test Coverage Status

.. |Documentation Status| image:: https://readthedocs.org/projects/lss/badge/?version=latest
    :target: https://LSS.readthedocs.io/en/latest/
    :alt: Documentation Status


Introduction
------------

This package contains codes to create LSS catalogs and randoms.

The code for k+e corrections was originally copied from Sam Moore's github (https://github.com/SgmAstro/DESI)

Installation
------------

First, git clone this repo:

If you want any change to the code to take place immediately, either:

1.  Add the "bin" directory to your
    ``$PATH`` environment variable and add the "py" directory to your
    ``$PYTHONPATH`` environment variable.
e.g., if you git cloned from your $HOME directory

.. code-block:: bash

    PYTHONPATH=$PYTHONPATH:$HOME/LSS/py
    PATH=$PATH:$HOME/LSS/bin # if needed

Or (not actually recommended as this is less common usage)

2.  Install (and uninstall) the current git checkout:

.. code-block:: bash

    pip install -e . --user

    pip uninstall LSS

You can also install a fixed version of the package:

.. code-block:: bash

    pip install -e . --user

This will put the main scripts bin/ in your $PATH (typically in $HOME/.local/bin).

This code only works on NERSC, in the DESI or cosmodesi environments. Most work in the DESI environment: e.g., run

.. code-block:: bash

    source /global/common/software/desi/desi_environment.sh main

Some code requires the cosmodesi environment to be loaded instead (e.g., for 2pt functions and/or reconstruction, regressis tools, etc.) Source it via

.. code-block:: bash

    source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main

or

.. code-block:: bash

    source /global/common/software/desi/users/adematti/cosmodesi_environment.sh test

for some of the latest capabilities.


Structure
---------

Versions developed for different specific settings are being created. These are in separate folders in the py/LSS directory. Common tools go in the LSS directory.

Scripts meant to act like executables are in the bin directory, but most of what gets run by the LSS pipeline is in the scripts directory. scripts/main/mkCat_main.py is the primary script for running LSS catalogs.

Examples
--------

The full LSS pipeline for a data version is scripts/dark_LSScat_pipe.sh

These are from early developments and might still work

python bin/mkCat_SV3_simp.py --type QSO #This re-creates the version 1 catalogs writing them into your scratch by default. One should be able to test edits to the catalogs using this (starting from some pre-defined inputs). Available types are QSO, BGS_ANY, BGS_BRIGHT, ELG, ELG_HIP, LRG, MWS_ANY

python bin/gatherSV_zinfo_alltiles.py --type ELG --release blanc #this gathers all of the redshift info for SV1 ELG targets in the blanc release (daily supported as release, LRG, QSO, BGS_ANY are additional supported types as should be anything in the SV1 DESIMASK)

python bin/mkCat_singletile.py --type ELG --tile 80623 --night deep #this creates clustering catalogs for ELGs (LRG, QSO, BGS_ANY all work) on tile 80623 (see https://desi.lbl.gov/trac/wiki/SurveyValidation/SV1#ObservedTiles for the initial list) using the deep coadd (other options are a particular night or all, though all is not recommended) from the blanc release redshift data (use --release cascades to update to the more recent data). There are other options that use 'n' or 'y' to toggle whether a particular stage in the process gets run.
