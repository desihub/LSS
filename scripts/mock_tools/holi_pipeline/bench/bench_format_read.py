#!/usr/bin/env python3
"""
Benchmark comparing READ performance between the .ecsv and .fits formats,
on files from the /pscratch/sd/j/jcolley/bench_io directory.

The .fits files are the ones already produced by bench_io_nersc.py in the
"fits_out" subdirectory (same tabular content as the .ecsv files, just a
different format).

For each file, we measure separately:
  - the raw disk read (bytes, no parsing)
  - the parsing/decoding into an astropy Table (CPU)
to see whether any gap between formats comes from disk or CPU.

Usage:
    python bench_format_read.py
    python bench_format_read.py --repeat 3 --csv results_format.csv
"""

import argparse
import glob
import io
import os
import statistics
import time

from astropy.table import Table

DIRECTORY = "/pscratch/sd/j/jcolley/bench_io"
FITS_SUBDIR = "fits_out"


def list_files(directory, pattern):
    return sorted(glob.glob(os.path.join(directory, pattern)))


def bench_format(name, files, fmt):
    read_io_times = []
    parse_cpu_times = []
    sizes = []

    for f in files:
        sizes.append(os.path.getsize(f))

        t0 = time.perf_counter()
        with open(f, "rb") as fh:
            raw_bytes = fh.read()
        t1 = time.perf_counter()

        Table.read(io.BytesIO(raw_bytes), format=fmt)
        t2 = time.perf_counter()

        read_io_times.append(t1 - t0)
        parse_cpu_times.append(t2 - t1)

    return {
        "name": name,
        "nb_files": len(files),
        "total_bytes": sum(sizes),
        "read_io_times": read_io_times,
        "parse_cpu_times": parse_cpu_times,
        "total_read_io_s": sum(read_io_times),
        "total_parse_cpu_s": sum(parse_cpu_times),
        "total_s": sum(read_io_times) + sum(parse_cpu_times),
    }


def fmt_stats(times):
    return (f"somme={sum(times):.3f}s mediane={statistics.median(times) * 1000:.1f}ms "
            f"max={max(times) * 1000:.1f}ms")


def print_result(result):
    mb = result["total_bytes"] / 1e6
    total = result["total_s"]
    print(f"\n=== {result['name']} ===")
    print(f"  nb files             : {result['nb_files']}")
    print(f"  total volume         : {mb:.1f} MB")
    print(f"  disk read (I/O)      : {fmt_stats(result['read_io_times'])}  "
          f"({mb / result['total_read_io_s']:.1f} MB/s)")
    print(f"  table parsing (CPU)  : {fmt_stats(result['parse_cpu_times'])}")
    print(f"  total time           : {total:.3f} s  ({mb / total:.1f} MB/s)")


def write_csv(results, csv_path):
    import csv

    fields = ["name", "nb_files", "total_bytes",
              "total_read_io_s", "total_parse_cpu_s", "total_s"]
    with open(csv_path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(fields)
        for r in results:
            writer.writerow([r[k] for k in fields])
    print(f"\nResults written to {csv_path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--pattern", default="*.ecsv",
                         help="glob pattern for the reference .ecsv files (default: *.ecsv)")
    parser.add_argument("--repeat", type=int, default=1,
                         help="number of times to repeat the benchmark (default: 1)")
    parser.add_argument("--csv", default=None,
                         help="path to a csv file to save the results")
    args = parser.parse_args()

    ecsv_files = list_files(DIRECTORY, args.pattern)
    fits_dir = os.path.join(DIRECTORY, FITS_SUBDIR)
    fits_files = [os.path.join(fits_dir, os.path.splitext(os.path.basename(f))[0] + ".fits")
                  for f in ecsv_files]
    missing = [f for f in fits_files if not os.path.exists(f)]
    if missing:
        print(f"Error: {len(missing)} .fits files missing in {fits_dir}.")
        print("Run bench_io_nersc.py first to generate the .fits files.")
        return

    all_results = []
    for it in range(args.repeat):
        if args.repeat > 1:
            print(f"\n########## Repetition {it + 1}/{args.repeat} ##########")
        results = [
            bench_format("ecsv", ecsv_files, "ascii.ecsv"),
            bench_format("fits", fits_files, "fits"),
        ]
        for r in results:
            print_result(r)
        all_results.extend(results)

        r_ecsv, r_fits = results
        faster = "fits" if r_fits["total_s"] < r_ecsv["total_s"] else "ecsv"
        ratio = max(r_ecsv["total_s"], r_fits["total_s"]) / min(r_ecsv["total_s"], r_fits["total_s"])
        print(f"\n>>> {faster} is faster to read (factor {ratio:.2f}x) "
              f"[ecsv: {r_ecsv['total_s']:.3f}s vs fits: {r_fits['total_s']:.3f}s]")

    if args.csv:
        write_csv(all_results, args.csv)


if __name__ == "__main__":
    main()
