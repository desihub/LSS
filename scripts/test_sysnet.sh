source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main

PYTHONPATH=$PYTHONPATH:$LSSDIR/LSS/py #make sure to set $LSSDIR to wherever the LSS repo is installed, e.g., export LSSCODE=$HOME

run_sysnet=sysnet-app

$run_sysnet #should output something like configuration file missing
