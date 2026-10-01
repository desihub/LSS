#!/global/common/software/desi/perlmutter/desiconda/20260227-2.3.1/conda/bin/python
###!/usr/bin/env python3
"""
Check the SLURM status of the array-job tasks of a Holi pipeline run.

The Holi pipeline is launched with `run_holi_pipeline.sh`, which submits
`sbatch_holi_pipeline.sh` as a SLURM job array. Each array task writes its
own SLURM output/error log in the run's log directory, named
`holi_<JOBID>.log` (see #SBATCH --output=holi_%j.log in
sbatch_holi_pipeline.sh). JOBID here is the individual task's SLURM job id,
not the array's master job id.

This script scans a run's log directory for `holi_<JOBID>.log` files, and
for each JOBID reports:
  * if the task is still known to SLURM (squeue): its state and how long
    it has been running (or pending);
  * otherwise (task finished and only visible via accounting): its final
    state, exit code, elapsed (wall) time, and whether it was killed by
    a time limit (TIMEOUT) or an out-of-memory condition (OUT_OF_MEMORY).

When every array task has COMPLETED, been CANCELLED, or hit TIMEOUT, the
script also builds a histogram (PNG) of the number of fba*.fits files
found (recursively) under each altmtl<XXXX> directory of the run (one
directory per seed, XXXX starting at the "first_id" parameter of the run's
*.toml file, and as many directories as there are logs/seed_*.log files),
saved as fba_files_histogram.png in the run's log directory.

Usage
-----
    check_holi_jobs.py /path/to/holi_260925_23h33
    check_holi_jobs.py /path/to/holi_260925_23h33/holi_58890835.log
"""

import argparse
import glob
import os
import re
import subprocess
import sys
import tomllib

LOG_NAME_RE = re.compile(r"^holi_(\d+)\.log$")

# sacct states worth calling out explicitly, in order of priority when a
# job has several steps (batch/extern/srun steps) with different states.
STATE_PRIORITY = ["OUT_OF_MEMORY", "TIMEOUT", "NODE_FAIL", "FAILED", "CANCELLED", "COMPLETED"]


def find_job_ids(log_path):
    """Return the sorted list of job ids found from a log dir or a single log file."""
    if os.path.isdir(log_path):
        paths = sorted(glob.glob(os.path.join(log_path, "holi_*.log")))
    elif os.path.isfile(log_path):
        paths = [log_path]
    else:
        return []

    job_ids = []
    for path in paths:
        m = LOG_NAME_RE.match(os.path.basename(path))
        if m:
            job_ids.append(m.group(1))
    return job_ids


def run_cmd(cmd):
    """Run a SLURM command, returning stdout, or None if it failed/is missing."""
    try:
        res = subprocess.run(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, check=False
        )
    except FileNotFoundError:
        return None
    if res.returncode != 0 or not res.stdout.strip():
        return None
    return res.stdout


def squeue_status(job_id):
    """Return (state, elapsed) if the job is still known to squeue, else None."""
    out = run_cmd(["squeue", "-h", "-j", job_id, "-o", "%T|%M"])
    if out is None:
        return None
    state, elapsed = out.strip().splitlines()[0].split("|")
    return state, elapsed


def sacct_status(job_id):
    """
    Return a dict (state, exitcode, elapsed, start, end) for a finished job,
    combining the worst state seen across its steps (batch/extern/srun) with
    the timing of its main allocation entry. Return None if sacct has no
    record of this job.

    NOTE: `sacct -j <raw_job_id>` also matches the raw job id used as the
    array's master id, which brings back *all* array tasks, not only the
    one we asked for. We therefore filter rows on JobIDRaw (ignoring the
    ".batch"/".extern"/".N" step suffix) to keep only the requested task.
    """
    out = run_cmd(
        [
            "sacct",
            "-j", job_id,
            "--noheader",
            "--parsable2",
            "--format=JobID,JobIDRaw,State,ExitCode,Elapsed,Start,End",
        ]
    )
    if out is None:
        return None

    main_row = None
    worst_state = None
    for line in out.strip().splitlines():
        fields = line.split("|")
        if len(fields) < 7:
            continue
        jobid, jobid_raw, state, exitcode, elapsed, start, end = fields[:7]
        if jobid_raw.split(".")[0] != job_id:
            continue  # belongs to a sibling array task, not this one
        if "." not in jobid_raw:
            main_row = {
                "jobid": jobid,
                "exitcode": exitcode,
                "elapsed": elapsed,
                "start": start,
                "end": end,
            }
        base_state = state.split()[0]  # drop e.g. "CANCELLED by 12345"
        if base_state in STATE_PRIORITY:
            if worst_state not in STATE_PRIORITY or STATE_PRIORITY.index(base_state) < STATE_PRIORITY.index(worst_state):
                worst_state = base_state
        elif worst_state is None:
            worst_state = base_state

    if main_row is None:
        return None
    main_row["state"] = worst_state or "UNKNOWN"
    return main_row


def describe_state(state):
    if state == "TIMEOUT":
        return "killed: time limit exceeded"
    if state == "OUT_OF_MEMORY":
        return "killed: out of memory"
    if state == "NODE_FAIL":
        return "killed: node failure"
    if state.startswith("CANCELLED"):
        return "cancelled"
    if state == "FAILED":
        return "failed"
    if state == "COMPLETED":
        return "completed successfully"
    return state


HISTOGRAM_FILENAME = "fba_files_histogram.png"


def find_params_toml(run_dir):
    """Return the path to the pipeline parameter file copied by run_holi_pipeline.sh into run_dir."""
    matches = sorted(glob.glob(os.path.join(run_dir, "*.toml")))
    return matches[0] if matches else None


def count_seed_logs(run_dir):
    """Number of seeds processed by the run = number of logs/seed_*.log files."""
    return len(glob.glob(os.path.join(run_dir, "logs", "seed_*.log")))


def count_fba_files(mock_dir, altmtl_id):
    """Count fba*.fits files anywhere under altmtl<id>, e.g. `find altmtl0050 -name 'fba*.fits'`."""
    cdir = f"altmtl{altmtl_id:04d}"
    pattern = os.path.join(mock_dir, cdir, "**", "fba*.fits")
    count= len(glob.glob(pattern, recursive=True))
    print(f'Found {count} fba*.fits files in {cdir}')
    return count


def plot_fba_histogram(run_dir):
    """
    Build a histogram of the number of fba*.fits files found (recursively)
    in each altmtlXXXX directory of the run, and save it as PNG in run_dir.
    Returns the PNG path, or None if it could not be built.
    """
    toml_path = find_params_toml(run_dir)
    if toml_path is None:
        print(f"No *.toml parameter file found in {run_dir}, skipping histogram", file=sys.stderr)
        return None

    with open(toml_path, "rb") as f:
        pars = tomllib.load(f)
    first_id = pars["first_id"]
    mock_dir = pars["mock_dir"]

    n_seeds = count_seed_logs(run_dir)
    if n_seeds == 0:
        print(f"No logs/seed_*.log file found in {run_dir}, skipping histogram", file=sys.stderr)
        return None

    counts = [count_fba_files(mock_dir, first_id + i) for i in range(n_seeds)]

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print(
            "matplotlib is not available in this Python; "
            "e.g. run 'source /global/common/software/desi/desi_environment.sh 26.3' first",
            file=sys.stderr,
        )
        return None

    last_id = first_id + n_seeds - 1
    fig, ax = plt.subplots()
    bins = range(min(counts), max(counts) + 2)
    ax.hist(counts, align="left", rwidth=0.8)
    ax.set_xlabel("number of fba*.fits files")
    ax.set_ylabel(f"number of altmtl directories (altmtl{first_id:04d}-altmtl{last_id:04d})")
    ax.set_title("Fiber assignment realizations (fba*.fits) per altmtl directory")
    ax.grid(axis="y", alpha=0.3)

    png_path = os.path.join(run_dir, HISTOGRAM_FILENAME)
    fig.savefig(png_path)
    plt.close(fig)
    return png_path


def report(log_path):
    job_ids = find_job_ids(log_path)
    if not job_ids:
        print(f"No holi_<jobid>.log file found in {log_path}", file=sys.stderr)
        return 1

    columns = ["JOB ID", "STATUS", "ELAPSED", "EXIT CODE", "START", "END", "DETAILS"]
    widths = [12, 14, 12, 10, 20, 20, 30]
    header = "".join(c.ljust(w) for c, w in zip(columns, widths))
    print(header)
    print("-" * len(header))

    all_completed = True
    for job_id in job_ids:
        sq = squeue_status(job_id)
        if sq is not None:
            state, elapsed = sq
            row = [job_id, state, elapsed, "-", "-", "-", f"running since {elapsed}"]
            all_completed = False
        else:
            sa = sacct_status(job_id)
            if sa is None:
                row = [job_id, "UNKNOWN", "-", "-", "-", "-", "not found in squeue/sacct"]
                all_completed = False
            else:
                row = [
                    job_id,
                    sa["state"],
                    sa["elapsed"],
                    sa["exitcode"],
                    sa["start"],
                    sa["end"],
                    describe_state(sa["state"]),
                ]
                # build the histogram for COMPLETED runs, and also for CANCELLED
                # or TIMEOUT ones (fiber assignment may already have produced
                # partial fba*.fits results before being killed)
                all_completed = all_completed and (
                    sa["state"] == "COMPLETED" or sa["state"].startswith("CANCELLED") or sa["state"] == "TIMEOUT"
                )
        print("".join(str(v).ljust(w) for v, w in zip(row, widths)))

    #if all_completed and os.path.isdir(log_path):
    png_path = plot_fba_histogram(log_path)
    if png_path is not None:
        print(f"\nAll array tasks completed, cancelled or timed out, histogram saved to {png_path}")

    return 0


def main():
    parser = argparse.ArgumentParser(
        description="Report the SLURM status (running / finished, exit code, "
        "timeout/OOM, elapsed time) of the array-job tasks of a Holi pipeline "
        "run, based on its holi_<jobid>.log files."
    )
    parser.add_argument(
        "log_path",
        help="path to the Holi pipeline run log directory "
        "(e.g. .../runs/holi_260925_23h33), or to a sinreportgle holi_<jobid>.log file",
    )
    args = parser.parse_args()

    if not os.path.exists(args.log_path):
        parser.error(f"{args.log_path} does not exist")

    sys.exit(report(args.log_path))


if __name__ == "__main__":
    main()
