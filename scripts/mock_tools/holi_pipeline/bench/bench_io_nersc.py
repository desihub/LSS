#!/usr/bin/env python3
"""
Disk read/write benchmark on the 2 NERSC filesystems (CFS and pscratch),
reading .ecsv files one by one and converting them (without modification)
into .fits files.

Usage:
    python bench_io_nersc.py
    python bench_io_nersc.py --pattern "*.ecsv" --repeat 3 --csv results.csv

The produced .fits files are written to a "fits_out" subdirectory of each
input directory, so that the write is tested on the same filesystem as the
read.
"""

import argparse
import glob
import io
import os
import time

from astropy.table import Table

# Directories to compare: (short name, path)
DEFAULT_DIRS = [
    ("pscratch", "/pscratch/sd/j/jcolley/bench_io"),
    ("cfs", "/global/cfs/cdirs/desi/users/colley/bench_io"),
]


def list_input_files(directory, pattern):
    files = sorted(glob.glob(os.path.join(directory, pattern)))
    return files


def bench_directory(name, directory, pattern, out_subdir):
    files = list_input_files(directory, pattern)
    if not files:
        print(f"[{name}] no file matching '{pattern}' in {directory}")
        return None

    out_dir = os.path.join(directory, out_subdir)
    os.makedirs(out_dir, exist_ok=True)

    # separate phases: raw disk I/O (read/write bytes) vs CPU
    # (ecsv parsing / fits serialization) done in memory on buffers.
    read_io_times = []
    parse_cpu_times = []
    serialize_cpu_times = []
    write_io_times = []
    sizes = []

    for f in files:
        sizes.append(os.path.getsize(f))

        t0 = time.perf_counter()
        with open(f, "rb") as fh:
            raw_bytes = fh.read()
        t1 = time.perf_counter()

        table = Table.read(io.BytesIO(raw_bytes), format="ascii.ecsv")
        t2 = time.perf_counter()

        buf = io.BytesIO()
        table.write(buf, format="fits", overwrite=True)
        t3 = time.perf_counter()

        out_file = os.path.join(out_dir, os.path.splitext(os.path.basename(f))[0] + ".fits")
        with open(out_file, "wb") as fh:
            fh.write(buf.getvalue())
        t4 = time.perf_counter()

        read_io_times.append(t1 - t0)
        parse_cpu_times.append(t2 - t1)
        serialize_cpu_times.append(t3 - t2)
        write_io_times.append(t4 - t3)

    total_bytes = sum(sizes)
    total_read_io = sum(read_io_times)
    total_parse_cpu = sum(parse_cpu_times)
    total_serialize_cpu = sum(serialize_cpu_times)
    total_write_io = sum(write_io_times)
    total_io = total_read_io + total_write_io
    total_cpu = total_parse_cpu + total_serialize_cpu

    result = {
        "name": name,
        "directory": directory,
        "nb_files": len(files),
        "total_bytes": total_bytes,
        "total_read_io_s": total_read_io,
        "total_parse_cpu_s": total_parse_cpu,
        "total_serialize_cpu_s": total_serialize_cpu,
        "total_write_io_s": total_write_io,
        "total_io_s": total_io,
        "total_cpu_s": total_cpu,
        "total_s": total_io + total_cpu,
    }
    return result


def print_result(result):
    mb = result["total_bytes"] / 1e6
    total = result["total_s"]
    print(f"\n=== {result['name']} ({result['directory']}) ===")
    print(f"  nb files             : {result['nb_files']}")
    print(f"  total volume         : {mb:.1f} MB")
    print(f"  disk read (I/O)      : {result['total_read_io_s']:.3f} s"
          f"  ({mb / result['total_read_io_s']:.1f} MB/s)")
    print(f"  ecsv parsing (CPU)   : {result['total_parse_cpu_s']:.3f} s")
    print(f"  fits serialization (CPU): {result['total_serialize_cpu_s']:.3f} s")
    print(f"  disk write (I/O)     : {result['total_write_io_s']:.3f} s"
          f"  ({mb / result['total_write_io_s']:.1f} MB/s)")
    print(f"  total disk I/O       : {result['total_io_s']:.3f} s"
          f"  ({100 * result['total_io_s'] / total:.0f}% of total time)")
    print(f"  total CPU (conversion): {result['total_cpu_s']:.3f} s"
          f"  ({100 * result['total_cpu_s'] / total:.0f}% of total time)")
    print(f"  total time           : {total:.3f} s"
          f"  ({mb / total:.1f} MB/s)")


def write_csv(results, csv_path):
    import csv

    fields = ["name", "directory", "nb_files", "total_bytes",
              "total_read_io_s", "total_parse_cpu_s", "total_serialize_cpu_s",
              "total_write_io_s", "total_io_s", "total_cpu_s", "total_s"]
    with open(csv_path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(fields)
        for r in results:
            writer.writerow([r[k] for k in fields])
    print(f"\nResults written to {csv_path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pattern", default="*.ecsv",
                         help="glob pattern for the input files (default: *.ecsv)")
    parser.add_argument("--out-subdir", default="fits_out",
                         help="output subdirectory for the .fits files")
    parser.add_argument("--repeat", type=int, default=1,
                         help="number of times to repeat the benchmark (default: 1)")
    parser.add_argument("--csv", default=None,
                         help="path to a csv file to save the results")
    args = parser.parse_args()

    all_results = []
    for it in range(args.repeat):
        if args.repeat > 1:
            print(f"\n########## Repetition {it + 1}/{args.repeat} ##########")
        results = []
        for name, directory in DEFAULT_DIRS:
            r = bench_directory(name, directory, args.pattern, args.out_subdir)
            if r is not None:
                print_result(r)
                results.append(r)
        all_results.extend(results)

        if len(results) == 2:
            r0, r1 = results
            faster = r0["name"] if r0["total_s"] < r1["total_s"] else r1["name"]
            ratio = max(r0["total_s"], r1["total_s"]) / min(r0["total_s"], r1["total_s"])
            print(f"\n>>> {faster} is faster (factor {ratio:.2f}x) "
                  f"[{r0['name']}: {r0['total_s']:.3f}s vs {r1['name']}: {r1['total_s']:.3f}s]")

    if args.csv:
        write_csv(all_results, args.csv)


if __name__ == "__main__":
    main()
