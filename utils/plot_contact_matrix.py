#!/usr/bin/env python3
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.colors as colors
from sys import exit, stderr, argv
import argparse
from pathlib import Path
import time
import matplotlib.ticker as plticker
import re


def parse_args():
    parser = argparse.ArgumentParser(description="Simple script to plot Hi-C distance matrices.")
    parser.add_argument("--bin-size", type=int, help="Bin size in base pairs",
                        required=True, dest="bin_size")
    parser.add_argument("--input-dir", type=Path,
                        help="Folder containing the (compressed) CSVs to be plotted.",
                        dest="input_dir")
    parser.add_argument("--output-dir", type=Path,
                        help="Output folder (it will be created if it doesn't already exists). Default: same as input.",
                        dest="output_dir")
    parser.add_argument("--force", action="store_true", dest="force",
                        help="Overwrite existing files.")
    parser.add_argument("--separator", type=str, help="Field separator (default='\\t').",
                        default="\t", dest="separator")
    parser.add_argument("--input-files", type=Path, nargs="+",
                        help="Path to (compressed) CSV file(s) to be plotted.", dest="input_files")
    parser.add_argument("--log-scale", action="store_true", dest="log_scale",
                        help="Plot contact counts in log-space.")

    return parser.parse_args()


def bp_to_mbp(n, pos):
    n *= bin_size / 1e6
    return f"{n:.2f} Mb"


if __name__ == "__main__":
    args = parse_args()
    input_files = []
    output_formats = ["png", "svg"]

    if args.input_dir:
        if not args.input_dir.is_dir():
            print(f"'{args.input_dir}' does not appear to be a valid directory", file=stderr)
            exit(1)
        input_files.extend(args.input_dir.glob("*.?sv*"))
    if args.input_files:
        input_files.extend([f for f in args.input_files])

    if len(input_files) == 0:
        print(f"Unable to find a suitable file among the ones provided as input.", file=stderr)
        exit(1)

    input_files.sort()

    bin_size = args.bin_size
    ax_tick_fmt = plticker.FuncFormatter(bp_to_mbp)

    if args.output_dir:
        args.output_dir.mkdir(parents=True, exist_ok=True)

    for f in input_files:
        skip = False
        for format in output_formats:
            if f.name.endswith(format):
                skip = True
                break
        base_out_dir = args.output_dir if args.output_dir else f.parent
        base_out_name = re.sub(r"\.[ct]sv\.*(gz|bz2|xz)?", "", f.name, flags=re.IGNORECASE)

        if (not args.force):
            for format in output_formats:
                if Path(f"{base_out_dir}/{base_out_name}.{format}").exists():
                    print(
                        f"File '{base_out_dir}/{base_out_name}.{format}' already exists. SKIPPING! Pass --force to overwrite.",
                        file=stderr)
                    skip = True
                    break
        if skip: continue

        t0 = time.time()
        print(f"Reading contact matrix from file '{f}'...", file=stderr, end="")
        c_matrix = pd.read_csv(f, sep="\t", header=None)
        mem_usage = c_matrix.memory_usage(index=True).sum() / 1e6
        print(
            f" DONE in {time.time() - t0}s! Read a {c_matrix.shape[0]}x{c_matrix.shape[1]} matrix using {mem_usage:.2f} MB of RAM.")
        t0 = time.time()
        c_matrix += 1
        # plt.ylim(c_matrix.shape[0] * bin_size)
        # plt.xlim(c_matrix.shape[1] * bin_size)
        fig, ax = plt.subplots(1, 1)
        ax.xaxis.set_major_formatter(ax_tick_fmt)
        ax.yaxis.set_major_formatter(ax_tick_fmt)
        plt.xticks(rotation=90)
        if args.log_scale:
            print(f" Plotting data in log-scale...", file=stderr)
            pcm = ax.imshow(c_matrix,
                            norm=colors.LogNorm(vmin=c_matrix.min().min(),
                                                vmax=c_matrix.max().max()),
                            cmap="hot", aspect="equal", interpolation=None)
        else:
            print(f" Plotting data...", file=stderr)
            pcm = ax.imshow(c_matrix, cmap="hot", aspect="equal", interpolation=None)

        fig.colorbar(pcm, ax=ax, extend="max")
        plt.tight_layout()
        for format in output_formats:
            outfile = f"{base_out_dir}/{base_out_name}.{format}"
            print(f"Saving plot {Path(outfile).name}...", file=stderr)
            plt.savefig(outfile, dpi=1000)
        plt.close(None)
        print(f" DONE in {time.time() - t0}s!", file=stderr)
