# Intro

This version of the pipeline is designed to generate at least 1,000 catalogs; it therefore uses SLURM's sbatch system rather than the interactive platform, although parts of the pipeline can still be run there. A single submission handles the entire production process, (TODO which is configurable via a parameter file.)



# Launch Holi pipeline

## 1) Create python environment

TBD

## 2) Prepare output mock directory, dataset DS_DIR

1. Create directory, in $SCRATCH place for example
2. copy nzref files from reference 

```console
cp /global/cfs/cdirs/desi/mocks/cai/holi/webjax_v4.80/nzref* /pscratch/sd/j/jcolley/holi/webjax_v4.80

login22:webjax_v4.80>ll /pscratch/sd/j/jcolley/holi/webjax_v4.80/
total 10K
drwxrwx--- 2 jcolley jcolley 4.0K Jul 16 16:00 .
drwxrwx--- 3 jcolley jcolley 4.0K Jul 16 15:50 ..
-rw-r----- 1 jcolley jcolley  11K Jul 16 16:00 nzref_da2_elg_N.txt
-rw-r----- 1 jcolley jcolley  11K Jul 16 16:00 nzref_da2_elg_S.txt
-rw-r----- 1 jcolley jcolley 7.8K Jul 16 16:00 nzref_da2_lrg.txt
-rw-r----- 1 jcolley jcolley  18K Jul 16 16:00 nzref_da2_qso.txt
```

note the definition of DS_DIR for step 4)

[perlmutter scratch doc](https://docs.nersc.gov/filesystems/perlmutter-scratch/)

## 3) Install Holi pipeline

### 3.1) Clone desihub/LSS package

```console
git clone https://github.com/desihub/LSS.git
cd LSS
LSS_DIR=$PWD
git checkout fa4acm
```

note the definition of LSS_DIR for step 4)

### 3.2) Copy/Modify file parameters

```console
cd scripts/mock_tools/holi_pipeline
```

## 4) Check parameters
TODO

modify this parameters

```console
#SBATCH --array=0-7
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=24
```
in sbatch1_array_step1234567.sh

adapt from step 2) and 3.1)

```console
LSS_DIR=
DS_DIR=
```

## 5) Launch pipeline

create a run directory for log, like
```console
run_a8_t10_c24
```
(a)rray of 8 jobs with 10 (t)asks with 24 (C)PUs.

and submit pipeline 

```console
sbatch <path>/<to>/<holip_pipeline>/sbatch1_array_step1234567.sh <ID first seed>
```

**Time execution information**
## 7) Results

