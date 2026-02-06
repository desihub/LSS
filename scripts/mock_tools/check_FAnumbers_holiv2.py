#!/usr/bin/env python3
"""
PEP8-cleaned version of scripts/mock_tools/check_FAnumbers_holiv2.py.

Counts objects by tracer bit flags in a series of FITS files and writes
a summary file with counts per realization.
"""
from __future__ import annotations

import argparse
import os
from typing import Iterable

import fitsio
import numpy as np

def count_tracers_in_files(
    mock_dir: str,
    start: int = 0,
    end_exclusive: int = 1001,
    out_name: str = "tracer_num.txt",
    filenames: Iterable[str] | None = None,
) -> None:
    """Count LRG/ELG/QSO objects in mock FITS files and write results.

    Parameters
    ----------
    mock_dir
        Directory containing files named forFA{n}.fits (when filenames is None).
    start
        First realization index (inclusive) when filenames is None.
    end_exclusive
        Last realization index (exclusive) when filenames is None.
    out_name
        Output filename (written inside mock_dir).
    filenames
        Optional iterable of filenames (full paths) to inspect instead of the
        default forFA{n}.fits pattern.
    """
    out_path = os.path.join(mock_dir, out_name)
    header = "#realization num_LRG num_ELG num_QSO\n"

    with open(out_path, "w") as out_f:
        out_f.write(header)

        if filenames is None:
            iterable = (os.path.join(mock_dir, f"forFA{n}.fits") for n in range(start, end_exclusive))
        else:
            iterable = iter(filenames)

        for path in iterable:
            if not os.path.isfile(path):
                # skip missing files silently
                continue

            data = fitsio.read(path, columns=["DESI_TARGET"])
            # bitwise selections; ensure correct operator precedence
            sel_lrg = (data["DESI_TARGET"] & 1) > 0
            sel_elg = (data["DESI_TARGET"] & 2) > 0
            sel_qso = (data["DESI_TARGET"] & 4) > 0

            # extract realization number if possible, otherwise use filename
            base = os.path.basename(path)
            if base.startswith("forFA") and base.endswith(".fits"):
                realization = base[len("forFA") : -len(".fits")]
            else:
                realization = base

            n_lrg = int(np.sum(sel_lrg))
            n_elg = int(np.sum(sel_elg))
            n_qso = int(np.sum(sel_qso))

            line = f"{realization} {n_lrg} {n_elg} {n_qso}\n"
            out_f.write(line)
            print(line.rstrip())

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Count tracer numbers in holi_v2 mock FITS files."
    )
    parser.add_argument(
        "--mock-dir",
        default="/pscratch/sd/d/desica/DA2/mocks/holi_v2/",
        help="Directory with mock files (default: %(default)s)",
    )
    parser.add_argument(
        "--start",
        type=int,
        default=0,
        help="First realization index (inclusive).",
    )
    parser.add_argument(
        "--end",
        type=int,
        default=1001,
        help="Last realization index (exclusive).",
    )
    parser.add_argument(
        "--out-name",
        default="tracer_num.txt",
        help="Output filename to write inside mock-dir.",
    )
    return parser.parse_args()

def main() -> None:
    args = parse_args()
    count_tracers_in_files(
        mock_dir=args.mock_dir,
        start=args.start,
        end_exclusive=args.end,
        out_name=args.out_name,
    )

if __name__ == "__main__":
    main()