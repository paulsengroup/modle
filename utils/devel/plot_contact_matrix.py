#!/usr/bin/env python3

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as colors
from sys import exit, stderr
import argparse
from pathlib import Path
import time
import matplotlib.ticker as plticker
import re
from palettable.colorbrewer.sequential import Reds_3
from typing import Union


def parse_args() -> None:
    parser = argparse.ArgumentParser(description="Simple script to plot Hi-C distance matrices.")
    parser.add_argument("--bin-size", type=int, help="Bin size in base pairs", required=True, dest="bin_size")
    parser.add_argument("--input-dir", type=Path, required=True,
                        help="Folder containing the (compressed) CSVs to be plotted.", dest="input_dir")
    parser.add_argument("--output-dir", type=Path,
                        help="Output folder (it will be created if it doesn't already exists). Default: same as input.",
                        dest="output_dir")
    parser.add_argument("--force", action="store_true", dest="force", help="Overwrite existing files.")
    parser.add_argument("--separator", type=str, help="Field separator (default='\\t').",
                        default="\t", dest="separator")
    parser.add_argument(
        "--upper-limit-linear-scale", type=int,
        help="Upper limit for the number of contacts when plotting in linear scale.", dest="ul_linear", default=17)

    return parser.parse_args()


def bp_to_mbp(n: int, pos: int) -> str:
    n *= bin_size / 1e6
    return f"{n:.2f} Mb"


"""Generate a template for the metadata to be embedded in the heatmaps"""


def generate_metadata(input_file: Union[Path, str]) -> dict:
    buf = ""
    metadata = {"Software": "modle"}
    with open(input_file, "r") as f:
        for line in f:
            if line == "\n":
                continue
            buf = line  # Copy CLI args used to generate input data
    metadata["Description"] = buf.strip()
    metadata["Comment"] = buf.strip()
    return metadata


def make_plot(matrix: pd.DataFrame, bin_size: int, base_name: Union[Path, str], metadata: dict,
              ax_tick_fmt, extr_barriers: [pd.DataFrame, None], log_scale: bool, linear_scale_upper_limit=0) -> None:
    if log_scale:
        print("Preprocessing the contact matrix data...", file=stderr, end="")
        matrix += 1  # Can't log(0)...
        print(" DONE!", file=stderr)

    # Extract forward and reverse extr. barriers (the coords of rev. barriers are given as negative numbers)
    print("Preprocessing extrusion barrier coordinates...", file=stderr, end="")
    extr_barriers = extr_barriers.to_numpy() if extr_barriers is not None else np.empty([0, 0])
    fwd_barrier_coords = extr_barriers[extr_barriers >= 0] / bin_size
    rev_barrier_coords = extr_barriers[extr_barriers < 0] / -bin_size
    print(" DONE!", file=stderr)

    fig, ax = plt.subplots(1, 1)

    # The custom formatter prints labels in Mbp
    ax.xaxis.set_major_formatter(ax_tick_fmt)
    ax.yaxis.set_major_formatter(ax_tick_fmt)

    # Setting various properties of the plot
    plt.title(metadata["Title"] + f" bin_size={bin_size}")
    plt.xticks(rotation=45)
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)

    scale = "log" if log_scale else "linear"
    pcm = None

    # Do the actual plotting
    matrix_type = "raw" if base_name.endswith("raw") else "complete"
    print(f"Plotting {matrix_type} matrix in {scale} scale...", file=stderr, end="")
    # Plot contacts
    if log_scale:
        pcm = ax.imshow(matrix,
                        norm=colors.LogNorm(vmin=matrix.min().min(),
                                            vmax=matrix.max().max()),
                        cmap="hot", aspect="equal", interpolation=None)
    else:
        pcm = ax.imshow(matrix,
                        vmin=1, vmax=linear_scale_upper_limit,
                        cmap=Reds_3.mpl_colormap,  # This cmap seems like a good approximation
                                                   # of Juicer's colormap
                        aspect="equal", interpolation=None)

    # Plot extrusion barriers
    if matrix_type == "raw":
        ax.scatter(fwd_barrier_coords, np.zeros_like(
            fwd_barrier_coords), color="blue", s=0.5, marker=">", edgecolors='none')
        ax.scatter(rev_barrier_coords, np.zeros_like(
            rev_barrier_coords), color="black", s=0.5, marker="<", edgecolors='none')
    else:
        ax.scatter(fwd_barrier_coords, fwd_barrier_coords,
                   color="blue", s=0.5, marker=">", edgecolors='none')
        ax.scatter(rev_barrier_coords, rev_barrier_coords,
                   color="black", s=0.5, marker="<", edgecolors='none')

    fig.colorbar(pcm, ax=ax, extend="max")  # Write the plot's scale
    # Write CLI settings used to generate the input data in the plot footer
    plt.figtext(0.5, 0.01, metadata["Description"], ha="center", fontsize=3, wrap=True)
    plt.tight_layout()
    print(" DONE!", file=stderr)

    # Write plots to file
    for ext in output_formats:
        print(f" - Saving plot in {ext.upper()} format...", file=stderr, end="")
        outfile = f"{base_name}_{scale}.{ext}"
        plt.savefig(outfile, dpi=1000)
        print(f" DONE!", file=stderr)
    plt.close(None)


def group_input_files(input_dir: Path) -> dict:
    input_files = None
    # This dictionary will have the following structure:
    # {
    #   path_to_full_cmatrix: [path_to_raw_cmatrix, path_to_extr_barrier_coords]
    #   ...
    #  }
    # If the raw contact matrix or the extrusion barrier coords are not available, their path will be set to None
    input_files_grouped = {}
    if input_dir:  # Glob input dir for mydir/*.?sv*
        if not input_dir.is_dir():
            print(
                f"'{input_dir.input_dir}' does not appear to be a valid directory", file=stderr)
            exit(1)
        input_files = set(input_dir.glob("*.?sv*"))

    if len(input_files) == 0:
        print(
            f"Unable to find suitable file(s) in folder '{input_dir}'.", file=stderr)
        exit(1)

    # Exclude files corresponding to raw contacts or extr. barrier coordinates.
    # Remember to update this pattern if we change the naming scheme of any of the (compressed) csv-like files
    exclude_pattern = re.compile(
        r"^.*(_raw|\.extrusion_barriers)\.[ct]sv\.*(gz|bz2|xz)?$", re.IGNORECASE)
    # Match allowed file extensions
    ext_pattern = re.compile(r"\.[ct]sv\.*(gz|bz2|xz)?$", re.IGNORECASE)
    for file in input_files:
        f = str(file)
        if exclude_pattern.match(f):
            continue
        files = [None, None]
        ext_len = len(ext_pattern.search(f).group(0))
        # Generate file name for the raw cmatrix and extr. barriers coords.
        # based on the file name of the complete contact matrix
        raw_contacts = Path(f"{f[:-ext_len]}_raw{f[-ext_len:]}")
        extr_barriers = Path(f"{f[:-ext_len]}.extrusion_barriers{f[-ext_len:]}")

        if raw_contacts in input_files:
            files[0] = raw_contacts
        if extr_barriers in input_files:
            files[1] = extr_barriers
        input_files_grouped[file] = files

    return input_files_grouped


if __name__ == "__main__":
    # Parse argv, glob/group input files and set some variables
    args = parse_args()
    output_formats = ["png", "svg"]
    metadata_template = generate_metadata(f"{args.input_dir}/settings.log")
    grouped_input_files = group_input_files(args.input_dir)
    bin_size = args.bin_size
    ax_tick_fmt = plticker.FuncFormatter(bp_to_mbp)

    if args.output_dir:
        args.output_dir.mkdir(parents=True, exist_ok=True)

    ext_pattern = re.compile(r"\.[ct]sv\.*(gz|bz2|xz)?", flags=re.IGNORECASE)
    for cmatrix_file, (raw_cmatrix_file, extr_barr_file) in grouped_input_files.items():
        skip = False
        base_out_dir = args.output_dir if args.output_dir else cmatrix_file.parent
        base_out_name = ext_pattern.sub("", cmatrix_file.name)

        if (not args.force):  # By default we do not process files that would overwrite existing files
            for ext in output_formats:
                if Path(f"{base_out_dir}/{base_out_name}.{ext}").exists():
                    print(
                        f"File '{base_out_dir}/{base_out_name}.{ext}' already exists. SKIPPING! Pass --force to overwrite.",
                        file=stderr)
                    skip = True
                    break
        if skip:
            continue

        t0 = time.time()
        print(f"Reading data into memory...", file=stderr, end="")
        # Read contact matrices and extr. barrier coordinates
        cmatrix = pd.read_csv(cmatrix_file, sep="\t", header=None)
        raw_cmatrix = pd.read_csv(
            raw_cmatrix_file, sep="\t", header=None) if raw_cmatrix_file is not None else None
        extrusion_barriers = pd.read_csv(
            extr_barr_file, sep="\t",  header=None).T if extr_barr_file is not None else None

        # Calculate and print memory usage
        mem_usage = cmatrix.memory_usage(index=True).sum()
        if raw_cmatrix is not None:
            mem_usage += raw_cmatrix.memory_usage(index=True).sum()
        if extrusion_barriers is not None:
            mem_usage += extrusion_barriers.memory_usage(index=True).sum()
        mem_usage /= 1e6
        print(f" DONE in {(time.time() - t0):.4f}s using a total of {mem_usage:.2f} MB of RAM!", file=stderr)

        t0 = time.time()
        # We are using as title the last two levels of the path to the complete contact matrix (without extension)
        title = ext_pattern.sub("", str(cmatrix_file))[::-1].replace("/", "_", 1)[::-1].rsplit("/", 1)[-1]
        metadata = metadata_template.copy()
        metadata["Title"] = title

        # Plot the complete matrix in linear and log scale
        make_plot(cmatrix, bin_size, f"{base_out_dir}/{base_out_name}", metadata,
                  ax_tick_fmt, extrusion_barriers, False, args.ul_linear)
        make_plot(cmatrix, bin_size, f"{base_out_dir}/{base_out_name}", metadata,
                  ax_tick_fmt, extrusion_barriers, True)
        # Plot the raw matrix in linear and log scale
        if raw_cmatrix is not None:
            make_plot(raw_cmatrix, bin_size, f"{base_out_dir}/{base_out_name}_raw",
                      metadata, ax_tick_fmt, extrusion_barriers, False, args.ul_linear)
            make_plot(raw_cmatrix, bin_size, f"{base_out_dir}/{base_out_name}_raw",
                      metadata, ax_tick_fmt, extrusion_barriers, True)
