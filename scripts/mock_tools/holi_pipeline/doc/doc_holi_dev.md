# Holi pipeline for developers

## SLURM use

The pipeline employs two levels of implicit loops provided by SLURM. The first relies on the job array feature via the `#SBATCH --array` directive.
Each task in the array can start at a different time.
The second level is handled by the `srun -n xx script.sh` command, which executes the `script.sh` script in parallel `xx` times.

During the development of the pipeline, it became apparent that using nested `srun` commands complicated the exclusive use of CPU resources; consequently, the pipeline was written to use only a single level of `srun`.

By default, if one of the srun tasks fails, the others are stopped; since the simulation calculations are independent, this srun behavior must be modified using the option `--kill-on-bad-exit=0 `.

See file [sbatch_holi_pipeline.sh](./sbatch_holi_pipeline.sh)

## DESI environment

### Versions

To ensure reproducibility, avoid the "main" tag, as the versions of the installed packages are not yet fixed. Command to get the list of available tags:
* for desi_environment.sh

```console
module use /global/common/software/desi/perlmutter/desiconda/startup/modulefiles; module avail desimodules
```

* for cosmodesi_environment.sh

```console
module use /global/common/software/desi/users/adematti/perlmutter/cosmodesiconda/startup/modulefiles; module avail cosmodesimodules
```

### Isolation of different environments

Sequencing two processes that require different environments in this manner is not recommended.

```console
#!/usr/bin/env bash

source env1.sh
processing1.sh

source env2.sh
processing2.sh
```

Defining env1 beforehand can have side effects on env2.

#### bash subshell solution with ( )

Bash has a solution for this: "subshells" using parentheses.
```console
#!/usr/bin/env bash

(
    source env1.sh
    processing1.sh
)

# no impact of "source env1.sh" here now
source env2.sh
processing2.sh
```

#### Subscripts solution

You can also use a script that sets up the environment and then launches the processing.

`do_processing1.sh`

```console
#!/usr/bin/env bash

source env1.sh
processing1.sh
```

then

```console
#!/usr/bin/env bash

do_processing1.sh

# no impact of "source env1.sh" here now
source env2.sh
processing2.sh
```

The downside is that it is possible to manage the local variables used by processing1.sh.

#### With SLURM

For high-level pipelines like Holi's, it is important to understand the rules for environment propagation with sbatch and srun. By default, SLURM generally passes on exported variables. Consequently, submitting a job using `sbatch` must be done in a neutral environment; in practice, this means avoiding having an environment sourced by default in your `.bashrc`.

There is an `--export` option that allows limiting variable transmission by providing a list, which is practically difficult to manage, while transmitting nothing via `--export=NONE` is too drastic.

If it is necessary to define an environment within the sbatch script, you must use Bash subshells (with parentheses) and include the environment initialization within the scripts launched by `srun`.


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

```