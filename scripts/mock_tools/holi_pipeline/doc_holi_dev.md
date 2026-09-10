# Holi pipeline for developers

## SLURM use

The pipeline employs two levels of implicit loops provided by SLURM. The first relies on the job array feature via the `#SBATCH --array` directive.
Each task in the array can start at a different time.
The second level is handled by the `srun -n xx script.sh` command, which executes the `script.sh` script in parallel `xx` times.

During the development of the pipeline, it became apparent that using nested `srun` commands complicated the exclusive use of CPU resources; consequently, the pipeline was written to use only a single level of `srun`.

By default, if one of the srun tasks fails, the others are stopped; since the simulation calculations are independent, this srun behavior must be modified using the option `--kill-on-bad-exit=0 `.


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
Use interactive session with --cpus-per-task=1, like:

```bash
salloc -N 1 -C cpu --ntasks=10 --cpus-per-task=1 -t 4:00:00 --qos interactive --account desi
```

You can't used directly script `run_holi_piepeline.sh` but you can launch  `sbatch_holi_pipeline.sh <file.toml> ` and the script will execute the first rank of array job, so here 10 seeds will processed.

## Tests

todo