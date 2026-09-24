# Holi pipeline for developers

## SLURM use

The pipeline employs two levels of implicit loops provided by SLURM. The first relies on the job array feature via the `#SBATCH --array` directive.
Each task in the array can start at a different time.
The second level is handled by the `srun -n xx script.sh` command, which executes the `script.sh` script in parallel `xx` times.

During the development of the pipeline, it became apparent that using nested `srun` commands complicated the exclusive use of CPU resources; consequently, the pipeline was written to use only a single level of `srun`.

By default, if one of the srun tasks fails, the others are stopped; since the simulation calculations are independent, this srun behavior must be modified using the option `--kill-on-bad-exit=0 `.

See file [sbatch_holi_pipeline.sh](./sbatch_holi_pipeline.sh)

## Namming file and good practice

It is best practice for a data processing script not to hardcode input and output filenames. Most steps in the Holi pipeline followed this principle prior to this update; to ensure consistency across all scripts, I rewrote the interfaces of the original scripts to use the same option names, aligning them with this parameter naming convention.


```python
parser = argparse.ArgumentParser()
parser.add_argument("--inputs", nargs="+")
parser.add_argument("--outputs", nargs="+")
args = parser.parse_args()
```

So, currently, the script names for the Holi pipeline steps use the original script name plus the suffix `_stdpars` (for "standard parameters").

## File parameters for script bash

To minimize the changes required to make the pipeline usable in a different context—or simply to allow a user to change the results directory—a parameter file proved necessary. The challenge lies in the fact that parameterization must be handled at the Bash script level. In fact, there is a simple solution using TOML files—which allow the various stages to be organized hierarchically—and a Python script that can be called from a Bash script as follows:

```bash
SURV=$(get_pars.py $HOLI_PARS prepare_mocks.survey)
```
where `$HOLI_PARS` is the file parameters and `prepare_mocks.survey` is the parameters `survey` in section `prepare_mocks` in TOML file


```toml
[prepare_mocks] # step 1
survey = "DA2"
```

>NOTE
>
> The python script `get_pars.py` is totaly independant of Holi pipeline and can use in different project. 

Currently, the pipeline uses only the key-value functionality of the TOML format. The use of lists and dictionaries has not been explored.

## Logs files

It proved convenient for tracking the processing status of a seed to have a separate log file for each one—something achievable with `srun`. However, the files are currently named `seed_xx_tyy.log` 

```bash
srun ...
--open-mode=append \
--output="${LOG_DIR}/logs/seed_${FIRST_ID_RANK}_t%t.log" \
--error="${LOG_DIR}/logs/seed_${FIRST_ID_RANK}_t%t.log" \
```

for the seed where `xx+yy=zz`; I haven't yet found a way to get them named `seed_zz.log`.

## Debug method

### Interactive session
Several hours can elapse between submitting the pipeline with `sbatch` and its actual start; using an interactive session for debugging is recommended.

Use interactive session with --cpus-per-task=1, like:

```bash
salloc -N 1 -C cpu --ntasks=10 --cpus-per-task=1 -t 4:00:00 --qos interactive --account desi
```

You can't used directly script `run_holi_piepeline.sh` but you can launch  `sbatch_holi_pipeline.sh <file.toml> ` and the script will execute the first rank of array job, so here 10 seeds will processed.

### Reduce number of galaxies in catalog

Step 1 takes into account the `max_gal` parameter—set in the parameter file—to limit the number of galaxies. This speeds up the pipeline up to step 7; however, I do not think it makes sense for step 8 (FA). The minimum value to avoid significantly disrupting processing is around 5,000.


## Tests

### step 7

```console
salloc -N 1 -C cpu --ntasks=10 --cpus-per-task=1 -t 4:00:00 --qos interactive --account desi


cd /global/cfs/cdirs/desi/users/colley/LSS/scripts/mock_tools/holi_pipeline
LSS_DIR=/global/cfs/cdirs/desi/users/colley/LSS
DS_DIR=/pscratch/sd/j/jcolley/test_hp
IDS=200
out6=$DS_DIR/$(printf "forFA%04d.fits" "$IDS")
echo $out6
ll $out6
source /global/common/software/desi/desi_environment.sh main
#module load desitarget/3.0.0
# use local package LSS, refresh after source env
export PYTHONPATH=$LSS_DIR/py:$PYTHONPATH
export PATH=$LSS_DIR/bin:$HOLI_DIR:$PATH
# desitarget ecsv to fits
target_dir=$(get_pars.py $HOLI_PARS TARGET_dir)
export PYTHONPATH=$target_dir/py:$PYTHONPATH
export PATH=$target_dir/bin:$PATH
altmtlxxxx=$(printf "altmtl%04d" "$IDS")
out7=$DS_DIR/$altmtlxxxx
echo $out7


cd /global/cfs/cdirs/desi/users/colley/LSS/scripts/mock_tools/holi_pipeline
source init_env_holi.sh 
6786  LSS_DIR=$(get_pars.py $HOLI_PARS LSS_dir)
 6787  DS_DIR=$(get_pars.py $HOLI_PARS mock_dir)
 6788  echo LSS_DIR
 6789  echo $LSS_DIR
 6790  target_dir=$(get_pars.py $HOLI_PARS TARGET_dir)
 6791  echo $target_dir
 6792  IDS=200
 6793  out6=$DS_DIR/$(printf "forFA%04d.fits" "$IDS")
 6794  echo $out6
 6795  DS_DIR=$(get_pars.py $HOLI_PARS mock_dir)
 6796  out6=$DS_DIR/$(printf "forFA%04d.fits" "$IDS")
 6797  echo $out6
 6798  ll /pscratch/sd/j/jcolley/holi/test_hp/forFA0200.fits
 6799  DS_DIR=$(get_pars.py $HOLI_PARS mock_dir)
 6800  out6=$DS_DIR/$(printf "forFA%04d.fits" "$IDS")
 6801  ll $out6
 6802  source /global/common/software/desi/desi_environment.sh main
 6803  #module load desitarget/3.0.0
 6804  # use local package LSS, refresh after source env
 6805  export PYTHONPATH=$LSS_DIR/py:$PYTHONPATH
 6806  export PATH=$LSS_DIR/bin:$HOLI_DIR:$PATH
 6807  # desitarget ecsv to fits
 6808  target_dir=$(get_pars.py $HOLI_PARS TARGET_dir)
 6809  export PYTHONPATH=$target_dir/py:$PYTHONPATH
 6810  export PATH=$target_dir/bin:$PATH
 6811  altmtlxxxx=$(printf "altmtl%04d" "$IDS")
 6812  out7=$DS_DIR/$altmtlxxxx
 6813  echo $out7
 6814  t6786  LSS_DIR=$(get_pars.py $HOLI_PARS LSS_dir)
 6787  DS_DIR=$(get_pars.py $HOLI_PARS mock_dir)
 6788  echo LSS_DIR
 6789  echo $LSS_DIR
 6790  target_dir=$(get_pars.py $HOLI_PARS TARGET_dir)
 6791  echo $target_dir
 6792  IDS=200
 6793  out6=$DS_DIR/$(printf "forFA%04d.fits" "$IDS")
 6794  echo $out6
 6795  DS_DIR=$(get_pars.py $HOLI_PARS mock_dir)
 6796  out6=$DS_DIR/$(printf "forFA%04d.fits" "$IDS")
 6797  echo $out6
 6798  ll /pscratch/sd/j/jcolley/holi/test_hp/forFA0200.fits
 6799  DS_DIR=$(get_pars.py $HOLI_PARS mock_dir)
 6800  out6=$DS_DIR/$(printf "forFA%04d.fits" "$IDS")
 6801  ll $out6
 6802  source /global/common/software/desi/desi_environment.sh main
 6803  #module load desitarget/3.0.0
 6804  # use local package LSS, refresh after source env
 6805  export PYTHONPATH=$LSS_DIR/py:$PYTHONPATH
 6806  export PATH=$LSS_DIR/bin:$HOLI_DIR:$PATH
 6807  # desitarget ecsv to fits
 6808  target_dir=$(get_pars.py $HOLI_PARS TARGET_dir)
 6809  export PYTHONPATH=$target_dir/py:$PYTHONPATH
 6810  export PATH=$target_dir/bin:$PATH
 6811  altmtlxxxx=$(printf "altmtl%04d" "$IDS")
 6812  out7=$DS_DIR/$altmtlxxxx
 6813  echo $out7
 6814  time ./initialize_amtl_mocks_da2_stdpars.py --inputs $out6 --outputs $out7 --obscon DARK
ime ./initialize_amtl_mocks_da2_stdpars.py --inputs $out6 --outputs $out7 --obscon DARK
```