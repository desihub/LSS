# Intro

This version of the pipeline is designed to generate at least 1,000 catalogs; it therefore uses SLURM's sbatch system rather than the interactive platform, although parts of the pipeline can still be run there. A single submission handles the entire production process, (TODO which is configurable via a parameter file.)



# Launch Holi pipeline

## 1) Create python environment

You need python >= 3.11 to read file parameters, like

```console
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main 
```
## 2) Prepare mock directory

1. Create directory, in $SCRATCH place to have best I/O performance 
2. Copy nzref files from reference

Actually the input reference is : `/global/cfs/cdirs/desi/mocks/cai/holi/webjax_v4.80`

You can configure it in file parameters, see 3.2

```console
cp /global/cfs/cdirs/desi/mocks/cai/holi/webjax_v4.80/nzref* /pscratch/sd/j/jcolley/holi/webjax_v4.80

login22:webjax_v4.80>ls /pscratch/sd/j/jcolley/holi/webjax_v4.80/
nzref_da2_elg_N.txt
nzref_da2_elg_S.txt
nzref_da2_lrg.txt
nzref_da2_qso.txt
```

note: this path is what you will set as `mock_dir` in the parameter file (see step 4)

[perlmutter scratch doc](https://docs.nersc.gov/filesystems/perlmutter-scratch/)

## 3) Install Holi pipeline

### 3.1) Clone desihub/LSS package

```console
git clone https://github.com/desihub/LSS.git
cd LSS
LSS_DIR=$PWD
git checkout fa4acm
```

note: this path is what you will set as `LSS_dir` in the parameter file (see step 4)

### 3.2) Copy/edit the parameter file

```console
cd scripts/mock_tools/holi_pipeline
cp holi_params.toml my_run_params.toml
```

Edit `my_run_params.toml` and set at least:

```toml
LSS_dir  = "<LSS_DIR, see step 3.1>"
mock_dir = "<DS_DIR, see step 2>"
first_id = <first seed ID to process>
```

and check the `[brickmask]` section (`cfitsio`, `exe_dir`, `conf_dir`)

## 4) Check parameters

Edit the `#SBATCH` directives at the top of `pipeline_step_1-8.sh`:

```console
#SBATCH --array=0-7
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=24
```

* `--array` selects how many array ranks (chunks of seeds) are submitted. The array range below only needs to start at 0: SBATCH --array=0-x. The first seed ID to process is set by the "first_id" parameter
* `--ntasks` is the number of seeds processed per array rank.
* `--cpus-per-task` selects the run mode:
  * `= 1`: **"full" mode**, steps 1 to 8 all run within the same job
    (Fiber Assignment, step 8, only uses 1 CPU per seed).
  * `> 1`: **"split" mode**, steps 1 to 7 can use the extra CPUs (e.g.
    Brickmask), but step 8 is automatically resubmitted as a separate
    job (`sbatch_step8.sh`) with fewer CPUs per task, so no CPU time is
    wasted during Fiber Assignment.

The first seed ID to process is no longer a command-line argument: it
is read from the `first_id` key of the parameter file (see step 3.2).

## 5) Launch pipeline

create a run directory for log, like
```console
run_a8_t10_c24
```
(a)rray of 8 jobs with 10 (t)asks with 24 (C)PUs.

may be copy your file parameter in this directory
and submit the pipeline, passing the path to the parameter file edited
in step 3.2:

```console
sbatch <path>/<to>/<holi_pipeline>/pipeline_step_1-8.sh <path_to>/my_run_params.toml
```

