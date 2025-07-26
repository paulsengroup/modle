#!/usr/bin/env python3

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import math
import pathlib
import sys

import pyBigWig


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser()

    cli.add_argument("reference", type=pathlib.Path, help="Path to the reference BigWig file.")
    cli.add_argument(
        "bigwig",
        type=pathlib.Path,
        help="Path to the BigWig file to compare with the reference.",
    )
    cli.add_argument(
        "--tol",
        type=float,
        default=0,
        help="Relative tolerance to use for FP comparison.",
    )
    cli.add_argument(
        "--quiet",
        action="store_true",
        default=False,
        help="Only output to stderr in case of errors.",
    )

    return cli


def compare_chromosomes(f1, f2):
    chroms1 = f1.chroms()
    chroms2 = f2.chroms()

    if chroms1 != chroms2:
        raise RuntimeError(f"Bigwig chromosomes differ!\n{chroms1}\n{chroms2}")


def compare_intervals(f1, f2, chrom_name: str, tol: float):
    if chrom_name not in f1.chroms():
        raise RuntimeError(f"{chrom_name} not found in reference file!")
    if chrom_name not in f2.chroms():
        raise RuntimeError(f"{chrom_name} not found in bigwig file!")

    intervals1 = f1.intervals(chrom_name, 0, f1.chroms()[chrom_name])
    intervals2 = f2.intervals(chrom_name, 0, f2.chroms()[chrom_name])

    if intervals1 is None:
        intervals1 = []
    if intervals2 is None:
        intervals2 = []

    if len(intervals1) != len(intervals2):
        raise RuntimeError(f"{chrom_name}: expected {len(intervals1)} intervals, found {len(intervals2)}!")

    for i, ((start1, end1, n1), (start2, end2, n2)) in enumerate(zip(intervals1, intervals2)):
        if start1 != start2 or end1 != end2 or not math.isclose(n1, n2, rel_tol=tol):
            raise RuntimeError(
                f"{chrom_name}:{start1}-{end1}: Found difference! Expected {start1}-{end1}={n1}, found {start2}-{end2}={n2}"
            )


def main() -> int:
    args = vars(make_cli().parse_args())
    path1 = args["reference"]
    path2 = args["bigwig"]

    if not args["quiet"]:
        print(f"Comparing {path1} with {path2}...", file=sys.stderr)

    f1 = pyBigWig.open(str(path1))
    f2 = pyBigWig.open(str(path2))

    files_are_identical = True

    try:
        compare_chromosomes(f1, f2)
    except RuntimeError as e:
        files_are_identical = False
        print(e, file=sys.stderr)

    for chrom in f1.chroms():
        try:
            compare_intervals(f1, f2, chrom, args["tol"])
        except RuntimeError as e:
            files_are_identical = False
            print(e, file=sys.stderr)

    if files_are_identical:
        if not args["quiet"]:
            print("Files are identical", file=sys.stderr)
    else:
        print("Files differ", file=sys.stderr)

    return 0 if files_are_identical else 1


if __name__ == "__main__":
    ec = main()
    sys.exit(ec)
