#!/usr/bin/env python3
"""
"Pure" disk I/O benchmark (no format conversion) on the 2 NERSC
filesystems (pscratch and cfs).

bench_io_nersc.py (ecsv read + fits conversion) turned out to be >90%
dominated by CPU (ecsv parsing, fits serialization), not by disk: the
pscratch/cfs comparison there is therefore not very relevant for judging
I/O.

This script isolates the actual filesystem performance:
  - reading the raw bytes of each file, in chunks (no parsing)
  - writing an unmodified copy of these bytes, with fsync() to force the
    actual write to disk (and not just to the kernel page cache)
  - --drop-cache option (posix_fadvise DONTNEED) after reading, to limit
    the effect of the page cache on subsequent repetitions
  - per-file statistics (median, p90, max), not just the sum, to reveal
    the "straggler" latencies typical of shared/parallel filesystems

Known limitations (NERSC, without root access):
  - impossible to flush the page cache globally (no sudo): a file already
    read during a previous run may be served again from the cache rather
    than re-read from disk. Focus mainly on the 1st pass (--repeat 1) or
    use distinct files between 2 comparisons.
  - CFS (GPFS) and pscratch (Lustre) each have their own server-side
    cache layers, outside of user control.

Usage:
    python bench_io_pure.py
    python bench_io_pure.py --repeat 3 --csv results_pure.csv
    python bench_io_pure.py --no-fsync --no-drop-cache
"""

import argparse
import glob
import os
import statistics
import time

DEFAULT_DIRS = [
    ("pscratch", "/pscratch/sd/j/jcolley/bench_io"),
    ("cfs", "/global/cfs/cdirs/desi/users/colley/bench_io"),
]

CHUNK_SIZE = 4 * 1024 * 1024  # 4 MB


def list_input_files(directory, pattern):
    return sorted(glob.glob(os.path.join(directory, pattern)))


def read_raw(path, drop_cache):
    fd = os.open(path, os.O_RDONLY)
    try:
        chunks = []
        t0 = time.perf_counter()
        while True:
            chunk = os.read(fd, CHUNK_SIZE)
            if not chunk:
                break
            chunks.append(chunk)
        elapsed = time.perf_counter() - t0
        if drop_cache:
            try:
                os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
            except (AttributeError, OSError):
                pass
    finally:
        os.close(fd)
    return elapsed, b"".join(chunks)


def write_raw(path, data, fsync):
    fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644)
    try:
        t0 = time.perf_counter()
        os.write(fd, data)
        if fsync:
            os.fsync(fd)
        elapsed = time.perf_counter() - t0
    finally:
        os.close(fd)
    return elapsed


def bench_directory(name, directory, pattern, out_subdir, drop_cache, fsync):
    files = list_input_files(directory, pattern)
    if not files:
        print(f"[{name}] no file matching '{pattern}' in {directory}")
        return None

    out_dir = os.path.join(directory, out_subdir)
    os.makedirs(out_dir, exist_ok=True)

    read_times, write_times, sizes = [], [], []
    for f in files:
        t_read, data = read_raw(f, drop_cache)
        out_file = os.path.join(out_dir, os.path.basename(f))
        t_write = write_raw(out_file, data, fsync)

        sizes.append(len(data))
        read_times.append(t_read)
        write_times.append(t_write)

    return {
        "name": name,
        "directory": directory,
        "nb_files": len(files),
        "total_bytes": sum(sizes),
        "read_times": read_times,
        "write_times": write_times,
        "total_read_s": sum(read_times),
        "total_write_s": sum(write_times),
        "total_s": sum(read_times) + sum(write_times),
    }


def fmt_stats(times):
    return (f"sum={sum(times):.3f}s median={statistics.median(times) * 1000:.1f}ms "
            f"p90={sorted(times)[int(0.9 * len(times))] * 1000:.1f}ms max={max(times) * 1000:.1f}ms")


def print_result(result):
    mb = result["total_bytes"] / 1e6
    total = result["total_s"]
    print(f"\n=== {result['name']} ({result['directory']}) ===")
    print(f"  nb files        : {result['nb_files']}")
    print(f"  total volume    : {mb:.1f} MB")
    print(f"  raw read        : {fmt_stats(result['read_times'])}  "
          f"({mb / result['total_read_s']:.1f} MB/s)")
    print(f"  raw write       : {fmt_stats(result['write_times'])}  "
          f"({mb / result['total_write_s']:.1f} MB/s)")
    print(f"  total I/O time  : {total:.3f} s  ({mb / total:.1f} MB/s)")


def write_csv(results, csv_path):
    import csv

    fields = ["name", "directory", "nb_files", "total_bytes",
              "total_read_s", "total_write_s", "total_s"]
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
                         help="glob pattern for the input files (default: *.ecsv)")
    parser.add_argument("--out-subdir", default="raw_copy",
                         help="output subdirectory for the raw copies")
    parser.add_argument("--repeat", type=int, default=1,
                         help="number of times to repeat the benchmark (default: 1)")
    parser.add_argument("--no-drop-cache", action="store_false", dest="drop_cache",
                         help="do not invalidate the page cache after reading")
    parser.add_argument("--no-fsync", action="store_false", dest="fsync",
                         help="do not force fsync() after writing")
    parser.add_argument("--csv", default=None,
                         help="path to a csv file to save the results")
    args = parser.parse_args()

    all_results = []
    for it in range(args.repeat):
        if args.repeat > 1:
            print(f"\n########## Repetition {it + 1}/{args.repeat} ##########")
        results = []
        for name, directory in DEFAULT_DIRS:
            r = bench_directory(name, directory, args.pattern, args.out_subdir,
                                 args.drop_cache, args.fsync)
            if r is not None:
                print_result(r)
                results.append(r)
        all_results.extend(results)

        if len(results) == 2:
            r0, r1 = results
            faster = r0["name"] if r0["total_s"] < r1["total_s"] else r1["name"]
            ratio = max(r0["total_s"], r1["total_s"]) / min(r0["total_s"], r1["total_s"])
            print(f"\n>>> {faster} is faster in pure I/O (factor {ratio:.2f}x) "
                  f"[{r0['name']}: {r0['total_s']:.3f}s vs {r1['name']}: {r1['total_s']:.3f}s]")

    if args.csv:
        write_csv(all_results, args.csv)


if __name__ == "__main__":
    main()
