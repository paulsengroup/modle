#!/usr/bin/env python3

# Copyright (C) 2023 Roberto Rossini (roberros@uio.no)
# SPDX-License-Identifier: MIT

import argparse
import pathlib
import shutil
import subprocess as sp
from typing import Union

import bioframe as bf
import pandas as pd
import pyBigWig


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser()
    cli.add_argument(
        "ref-genome",
        type=pathlib.Path,
        help="Path to a reference genome in FASTA format.",
    )
    cli.add_argument(
        "motif-matrix",
        type=pathlib.Path,
        help="Path to CTCF motif matrix in MEME format (e.g. https://jaspar.genereg.net/matrix/MA0139.1/).",
    )
    cli.add_argument(
        "chip-bwig",
        type=pathlib.Path,
        help="Path to a bigWig file with CTCF fold-change over control from ChIP-seq (e.g. https://www.encodeproject.org/files/ENCFF851FYN/).",
    )
    cli.add_argument(
        "--narrow-peaks",
        type=pathlib.Path,
        help="Path to a narrowPeak file (e.g. https://www.encodeproject.org/files/ENCFF330SHG/).",
    )
    cli.add_argument(
        "--output-prefix",
        type=pathlib.Path,
        required=True,
        help="Path prefix to use for output.",
    )
    cli.add_argument(
        "--coords",
        type=str,
        help="Genomic coordinates in UCSC format to include in the output files.",
    )

    cli.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Overwrite existing files (if any).",
    )

    return cli


def handle_path_collisions(*paths: pathlib.Path) -> None:
    collisions = [p for p in paths if p.exists()]

    if len(collisions) != 0:
        collisions = "\n - ".join((str(p) for p in collisions))
        raise RuntimeError(
            "Refusing to overwrite file(s):\n" f" - {collisions}\n" "Pass --force to overwrite existing file(s)."
        )


def get_apptainer() -> pathlib.Path:
    if cmd := shutil.which("apptainer"):
        return pathlib.Path(cmd)
    if cmd := shutil.which("singularity"):
        return pathlib.Path(cmd)
    raise RuntimeError("Unable to find Singularity/Apptainer in your PATH!")


def run_meme(
    path_to_fasta: pathlib.Path,
    path_to_motif: pathlib.Path,
    img="docker://memesuite/memesuite:5.5.2",
) -> pd.DataFrame:
    # TODO: it would be good to filter FASTA sequences based on --coords before running MAST
    cmd = [
        str(get_apptainer()),
        "exec",
        "-B",
        f"{path_to_fasta}:/data/ref.fa:ro",
        "-B",
        f"{path_to_motif}:/data/motif.meme",
        str(img),
        "mast",
        "-hit_list",
        "/data/motif.meme",
        "/data/ref.fa",
    ]
    with sp.Popen(" ".join(cmd), stdin=None, stderr=None, stdout=sp.PIPE, shell=True) as mast:
        df = pd.read_table(
            mast.stdout,
            delim_whitespace=True,
            comment="#",
            names=["chrom", "strand", "id", "name", "start", "end", "score", "pval"],
        )

        mast.communicate()
        if (code := mast.returncode) != 0:
            raise RuntimeError(f"{cmd} terminated with code {code}")

    df["score"] = 0.0
    df["strand"] = df["strand"].map({-1: "-", 1: "+"})
    df = df[["chrom", "start", "end", "name", "score", "strand"]]
    return df


def filter_candidate_binding_sites(path_to_narrowpeak: Union[pathlib.Path, None], df: pd.DataFrame) -> pd.DataFrame:
    if path_to_narrowpeak is None:
        return df

    cols = df.columns.tolist()
    df = bf.overlap(df, bf.read_table(path_to_narrowpeak, schema="narrowPeak"), suffixes=("", "__"))

    return df.dropna()[cols]


def generate_bwig(path_to_bigwig: pathlib.Path, output_name: pathlib.Path, coords: Union[str, None]):
    if coords is None:
        shutil.copyfile(path_to_bigwig, output_name)
        return

    chrom, start, end = bf.parse_region(coords)
    with pyBigWig.open(str(path_to_bigwig), "r") as bwin, pyBigWig.open(str(output_name), "w") as bwout:
        chrom_size = bwin.chroms()[chrom]
        bwout.addHeader([(chrom, chrom_size)])

        starts = []
        ends = []
        values = []

        for start, end, value in bwin.intervals(chrom, start, end):
            starts.append(start)
            ends.append(end)
            values.append(value)

        bwout.addEntries([chrom] * len(starts), starts, ends=ends, values=values)


def main():
    args = vars(make_cli().parse_args())
    out_prefix = args["output_prefix"]

    if not args["force"]:
        handle_path_collisions(out_prefix.with_suffix(".bed.gz"), out_prefix.with_suffix(".bw"))

    out_prefix.parent.mkdir(exist_ok=True, parents=True)

    df = run_meme(args["ref-genome"], args["motif-matrix"])
    if (coords := args["coords"]) is not None:
        df = bf.select(df, coords)

    df = filter_candidate_binding_sites(args["narrow_peaks"], df)
    df.to_csv(out_prefix.with_suffix(".bed.gz"), sep="\t", index=False, header=False)

    generate_bwig(args["chip-bwig"], out_prefix.with_suffix(".bw"), args["coords"])


if __name__ == "__main__":
    main()
