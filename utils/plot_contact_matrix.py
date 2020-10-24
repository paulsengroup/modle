#!/usr/bin/env python3
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.colors as colors
from sys import exit, stderr
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
    parser.add_argument("--upper-limit-linear-scale", type=int,
                        help="Upper limit for the number of contacts when plotting in linear scale.",
                        dest="ul_linear", default=17)

    return parser.parse_args()


def bp_to_mbp(n, pos):
    n *= bin_size / 1e6
    return f"{n:.2f} Mb"


def generate_metadata_string(input_file):
    buf = ""
    metadata = {"Software": "modle"}
    with open(input_file, "r") as f:
        for line in f:
            if line == "\n":
                continue
            buf = line
    metadata["Description"] = buf.strip()
    metadata["Comment"] = buf.strip()
    return metadata


def make_plot(matrix, bin_size, base_name, metadata, ax_tick_fmt, log_scale,
              linear_scale_upper_limit):
    fig, ax = plt.subplots(1, 1)
    ax.xaxis.set_major_formatter(ax_tick_fmt)
    ax.yaxis.set_major_formatter(ax_tick_fmt)
    plt.xticks(rotation=45)
    plt.title(metadata["Title"] + f" bin_size={bin_size}")
    scale = "log" if log_scale else "linear"
    print(f"Plotting data in {scale} scale...", file=stderr, end="")
    pcm = None
    if log_scale:
        pcm = ax.imshow(matrix,
                        norm=colors.LogNorm(vmin=matrix.min().min(),
                                            vmax=matrix.max().max()),
                        cmap="hot", aspect="equal", interpolation=None)
    else:
        pcm = ax.imshow(matrix,
                        vmin=0, vmax=linear_scale_upper_limit,
                        cmap="hot", aspect="equal", interpolation=None)

    fig.colorbar(pcm, ax=ax, extend="max")
    plt.figtext(0.5, 0.01, metadata["Description"], ha="center", fontsize=3, wrap=True)
    plt.tight_layout()
    print(" DONE!", file=stderr)
    for ext in output_formats:
        outfile = f"{base_name}_{scale}.{ext}"
        print(f"Saving plot in {ext.upper()} format...", file=stderr)
        plt.savefig(outfile, dpi=1000)
    plt.close(None)


if __name__ == "__main__":
    args = parse_args()
    input_files = []
    output_formats = ["png", "svg"]

    metadata_template = generate_metadata_string(f"{args.input_dir}/settings.log")

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
        for ext in output_formats:
            if f.name.endswith(ext):
                skip = True
                break
        base_out_dir = args.output_dir if args.output_dir else f.parent
        base_out_name = re.sub(r"\.[ct]sv\.*(gz|bz2|xz)?", "", f.name, flags=re.IGNORECASE)

        if (not args.force):
            for ext in output_formats:
                if Path(f"{base_out_dir}/{base_out_name}.{ext}").exists():
                    print(
                        f"File '{base_out_dir}/{base_out_name}.{ext}' already exists. SKIPPING! Pass --force to overwrite.",
                        file=stderr)
                    skip = True
                    break
        if skip: continue

        t0 = time.time()
        print(f"Reading contact matrix from file '{f}'...", file=stderr, end="")
        c_matrix = pd.read_csv(f, sep="\t", header=None)
        mem_usage = c_matrix.memory_usage(index=True).sum() / 1e6
        print(
            f" DONE in {time.time() - t0}s!\nRead a {c_matrix.shape[0]}x{c_matrix.shape[1]} matrix using {mem_usage:.2f} MB of RAM.",
            file=stderr)
        t0 = time.time()
        c_matrix += 1
        title = str(f).replace(".tsv.bz2", "")[::-1].replace("/", "_", 1)[::-1].rsplit("/", 1)[-1]
        metadata = metadata_template.copy()
        metadata["Title"] = title
        make_plot(c_matrix, bin_size, f"{base_out_dir}/{base_out_name}", metadata, ax_tick_fmt,
                  False, args.ul_linear)
        make_plot(c_matrix, bin_size, f"{base_out_dir}/{base_out_name}", metadata, ax_tick_fmt,
                  True, 0)

        print(f"DONE in {time.time() - t0}s!", file=stderr)
