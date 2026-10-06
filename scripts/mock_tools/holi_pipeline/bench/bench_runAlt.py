#!/usr/bin/env python3
"""
Extract timing information from a runAlt log file.

For each matching keyword, the duration string formatted like a Python
``timedelta`` (e.g. ``0:01:22.143200``) is extracted with the same
field position as these awk commands:

    awk '/fba_run/ {print $4}' seed_50_t0.log
    awk '/runAlt/ {print $7}' seed_50_t0.log

and converted to seconds as a numpy array of floats.
"""

import argparse
import re

import numpy as np
import matplotlib.pyplot as plt

# 1-based (awk style) field index holding the duration string for each keyword
KEYWORDS = {
    "fba_run": 7,
    "runAlt": 4,
}

# matches "H:MM:SS.ffffff", optionally prefixed by "D day(s), " (str(timedelta) format)
_DURATION_RE = re.compile(r"^(?:(\d+) days?, )?(\d+):(\d{2}):(\d{2}(?:\.\d+)?)$")


def duration_to_seconds(duration_str):
    """Convert a 'H:MM:SS.ffffff' (optionally 'D days, H:MM:SS.ffffff') string to seconds."""
    m = _DURATION_RE.match(duration_str.strip())
    if not m:
        raise ValueError(f"Unrecognized duration format: {duration_str!r}")
    days, hours, minutes, seconds = m.groups()
    total = float(seconds) + int(minutes) * 60 + int(hours) * 3600
    if days:
        total += int(days) * 86400
    return total


def extract_durations(log_file, keyword):
    """Return the durations (in seconds) of all lines matching `keyword`, as a numpy array."""
    field_index = KEYWORDS[keyword] - 1  # awk fields are 1-based
    durations = []
    with open(log_file, "r") as f:
        for line in f:
            if keyword not in line:
                continue
            #print(f"Found line with keyword '{keyword}': {line.strip()}")
            fields = line.split()
            #print(f"Fields: {fields}")
            if field_index >= len(fields):
               
                continue
            # print(f"Line: {line.strip()}")
            # print(fields[field_index])
            durations.append(duration_to_seconds(fields[field_index]))
    return np.array(durations, dtype=float)

def plot_durations(fba, no_fba):
    plt.figure(figsize=(10, 6))
    plt.hist(fba, bins=30, alpha=0.5, label='fba_run', color='blue')
    plt.hist(no_fba, bins=30, alpha=0.5, label='runAlt - fba_run', color='orange')
    plt.xlabel('Duration (seconds)')
    plt.ylabel('Frequency')
    plt.title('Distribution of Durations')
    plt.legend()
    plt.grid(axis='y', alpha=0.75)
    plt.tight_layout()
   


def plot_durations2(fba, no_fba):
    """Scatter of fba vs no_fba (log-log) with marginal log-scale histograms."""
    fig = plt.figure(figsize=(8, 8))
    gs = fig.add_gridspec(
        2, 2, width_ratios=(4, 1), height_ratios=(1, 4),
        left=0.1, right=0.95, bottom=0.1, top=0.95,
        wspace=0.05, hspace=0.05,
    )
    ax = fig.add_subplot(gs[1, 0])
    ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)
    ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)

    ax.scatter(no_fba, fba, s=10, alpha=0.5, color='blue')
    # ax.set_xscale('log')
    # ax.set_yscale('log')
    ax.set_xlabel('no_fba duration (seconds)')
    ax.set_ylabel('fba_run duration (seconds)')
    ax.grid(True, which='both', alpha=0.3)

    # xbins = np.logspace(np.log10(fba.min()), np.log10(fba.max()), 30)
    # ybins = np.logspace(np.log10(no_fba.min()), np.log10(no_fba.max()), 30)


    ax_histx.hist(no_fba, color='blue', alpha=0.5)
    # #ax_histx.set_xscale('log')
    ax_histx.set_yscale('log')
    # ax_histx.tick_params(axis='x', labelbottom=False)
    ax_histx.set_ylabel('Frequency')
    ax_histx.grid(True, which='both', alpha=0.3)

    ax_histy.hist(fba, color='orange', alpha=0.5, orientation='horizontal')
    # ax_histy.set_yscale('log')
    ax_histy.set_xscale('log')
    # ax_histy.tick_params(axis='y', labelleft=False)
    ax_histy.set_xlabel('Frequency')
    ax_histy.grid(True, which='both', alpha=0.3)

    fig.suptitle(f'RunAltMTL has 2 parts: fba_run and no_fba\nUse {len(no_fba)} points', fontsize=14)
  

def plot_cumulative(fba, no_fba):
    """Plot the cumulative sum (in seconds) of fba and no_fba, in run order."""
    plt.figure(figsize=(10, 6))
    plt.plot(np.cumsum(fba), label='fba_run', color='blue')
    plt.plot(np.cumsum(no_fba), label='runAlt - fba_run', color='orange')
    plt.xlabel('fba_xxx.fits file index')
    plt.ylabel('Cumulative duration (seconds)')
    plt.title('Cumulative Duration')
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()

def main():
    parser = argparse.ArgumentParser(
        description="Extract fba_run/runAlt durations (in seconds) from a runAlt log file."
    )
    parser.add_argument("log_file", help="Path to the log file to parse")
    args = parser.parse_args()

    res = {}
    for keyword in KEYWORDS:
        durations = extract_durations(args.log_file, keyword)
        print(f"{keyword}: {len(durations)} values")
        print(durations)
        if durations.size:
            print(f"  mean={durations.mean():.3f}s  max={durations.max():.3f}s  sum={durations.sum():.3f}s")
        res[keyword] = durations
    nofba_run = res['runAlt'] - res['fba_run']
    print(f"runAlt - fba_run: {len(nofba_run)} values")
    print(f"  mean={nofba_run.mean():.3f}s median={np.median(nofba_run):.3f}s  max={nofba_run.max():.3f}s  sum={nofba_run.sum():.3f}s")
    #plot_durations(res['fba_run'], nofba_run)
    plot_durations2(res['fba_run'], nofba_run)
    plot_cumulative(res['fba_run'], nofba_run)
    
if __name__ == "__main__":
    main()
    plt.show()
