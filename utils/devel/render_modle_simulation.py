#!/usr/bin/env python3

import argparse
import functools
import io
import multiprocessing as mp
import pathlib
import shutil
import subprocess as sp
from typing import List, Tuple, Union

import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import pandas as pd
import seaborn as sns

# Source: https://github.com/open2c/cooltools/blob/master/cooltools/lib/plotting.py
PALETTES = {
    "fall": np.array(
        (
            (255, 255, 255),
            (255, 255, 204),
            (255, 237, 160),
            (254, 217, 118),
            (254, 178, 76),
            (253, 141, 60),
            (252, 78, 42),
            (227, 26, 28),
            (189, 0, 38),
            (128, 0, 38),
            (0, 0, 0),
        )
    )
    / 255,
    "fall_mod": np.array(
        (
            (255, 255, 255),
            (255, 255, 204),
            (255, 237, 160),
            (254, 217, 118),
            (254, 178, 76),
            (253, 141, 60),
            (252, 78, 42),
            (227, 26, 28),
            (227, 26, 28),
            (227, 26, 28),
            (189, 0, 38),
        )
    )
    / 255,
    "blues": np.array(
        (
            (255, 255, 255),
            (180, 204, 225),
            (116, 169, 207),
            (54, 144, 192),
            (5, 112, 176),
            (4, 87, 135),
            (3, 65, 100),
            (2, 40, 66),
            (1, 20, 30),
            (0, 0, 0),
        )
    )
    / 255,
    "acidblues": np.array(
        (
            (255, 255, 255),
            (162, 192, 222),
            (140, 137, 187),
            (140, 87, 167),
            (140, 45, 143),
            (120, 20, 120),
            (90, 15, 90),
            (60, 10, 60),
            (30, 5, 30),
            (0, 0, 0),
        )
    )
    / 255,
    "nmeth": np.array(
        (
            (236, 250, 255),
            (148, 189, 217),
            (118, 169, 68),
            (131, 111, 43),
            (122, 47, 25),
            (41, 0, 20),
        )
    )
    / 255,
}


@functools.cache
def get_ffmpeg_args(
    ffmpeg_path: Union[None, pathlib.Path],
    encoder: str,
    fps: int,
    crf: int,
    nthreads: int,
) -> List[str]:
    if ffmpeg_path is None:
        raise RuntimeError("Unable to find ffmpeg!")
    return [
        str(ffmpeg_path),
        "-framerate",
        str(fps),
        "-f",
        "image2pipe",
        "-vcodec",
        "png",
        "-i",
        "-",
        "-c:v",
        encoder,
        "-threads",
        str(nthreads),
        "-crf",
        str(crf),
        "-pix_fmt",
        "yuv420p",
    ]


def list_to_colormap(color_list, name=None):
    color_list = np.array(color_list)
    if color_list.min() < 0:
        raise ValueError("Colors should be 0 to 1, or 0 to 255")
    if color_list.max() > 1.0:
        if color_list.max() > 255:
            raise ValueError("Colors should be 0 to 1 or 0 to 255")
        else:
            color_list = color_list / 255.0
    return mpl.colors.LinearSegmentedColormap.from_list(name, color_list, 256)


def get_cmap(name):
    is_reversed = name.endswith("_r")
    try:
        if is_reversed:
            pal = PALETTES[name[:-2]][::-1]
        else:
            pal = PALETTES[name]
    except KeyError:
        raise ValueError('Palette not found "{}"'.format(name))
    return list_to_colormap(pal)


def _register_cmaps():
    for name, pal in PALETTES.items():
        mpl.colormaps.register(name=name, cmap=list_to_colormap(pal))
        mpl.colormaps.register(name=name + "_r", cmap=list_to_colormap(pal[::-1]))


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser()

    def positive_int(arg):
        if (n := int(arg)) > 0:
            return n

        raise ValueError("Not a positive integer")

    cli.add_argument(
        "tsv",
        type=pathlib.Path,
        help="Path to the interaction file in TSV format generated by MoDLE.",
    )
    cli.add_argument(
        "resolution", type=positive_int, help="Resolution in bp to use for plotting."
    )
    cli.add_argument(
        "output-video", type=pathlib.Path, help="Path where to store the output video."
    )
    cli.add_argument(
        "--max-epochs",
        type=positive_int,
        help="Max number of epochs to render. When not specified, render all available epochs.",
    )
    cli.add_argument(
        "--step",
        type=positive_int,
        default=1,
        help="Render one frame every --step epoch.",
    )
    cli.add_argument(
        "--video-encoder",
        type=str,
        default="libx264",
        help="Video encoder to pass to ffmpeg.",
    )
    cli.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Force overwrite existing files.",
    )
    cli.add_argument(
        "--nproc",
        type=int,
        choices=range(2, mp.cpu_count() + 1),
        default=mp.cpu_count(),
        help="Maximum number of parallel processes.",
    )
    cli.add_argument(
        "--encoding-threads",
        type=int,
        choices=range(1, min(16, mp.cpu_count()) + 1),
        default=min(16, mp.cpu_count()),
        help="Maximum number of threads to use for encoding.",
    )
    cli.add_argument(
        "--fps",
        type=positive_int,
        default=240,
        help="Frames per second used for rendering.",
    )
    cli.add_argument(
        "--crf",
        type=positive_int,
        default=26,
        help="Constant rate factor (CRF) used for rendering.",
    )
    cli.add_argument(
        "--ffmpeg",
        type=pathlib.Path,
        default="ffmpeg",
        help="Path to ffmpeg binary (required if ffmpeg is not in PATH).",
    )

    return cli


def handle_path_collisions(*paths: pathlib.Path) -> None:
    collisions = [p for p in paths if p.exists()]

    if len(collisions) != 0:
        collisions = "\n - ".join((str(p) for p in collisions))
        raise RuntimeError(
            "Refusing to overwrite file(s):\n"
            f" - {collisions}\n"
            "Pass --force to overwrite existing file(s)."
        )


def compute_boundaries(df: pd.DataFrame, resolution: int) -> Tuple[int, int]:
    start_pos = (min(df["pos1"].min(), df["pos2"].min()) // resolution) * resolution
    end_pos = (
        (max(df["pos1"].max(), df["pos2"].max()) + resolution - 1) // resolution
    ) * resolution

    return start_pos, end_pos


def df_to_numpy(df, shape) -> npt.NDArray[int]:
    interactions = df.groupby(["bin1", "bin2"])["n"].sum()
    idx1 = interactions.index.get_level_values("bin1")
    idx2 = interactions.index.get_level_values("bin2")

    m = np.zeros([shape, shape], dtype=int)
    m[idx1, idx2] = interactions
    return np.triu(m, 1) + m.T


def plot_heatmap(
    m: npt.NDArray[int],
    vmin,
    vmax,
    cmap,
    ax,
    log=True,
    plot_colorbar=False,
):
    if log:
        sns.heatmap(
            m,
            norm=colors.LogNorm(vmin, vmax),
            cmap=cmap,
            ax=ax,
            cbar=plot_colorbar,
        )
    else:
        sns.heatmap(m, vmin=vmin, vmax=vmax, cmap=cmap, ax=ax, cbar=plot_colorbar)


def plot_heatmap_worker(m1, m2, i, vmax, cmap1, cmap2) -> bytes:
    m2[m2 == 0] = -1

    fig, ax = plt.subplots(1, 1)
    plot_heatmap(m1, 1, vmax, cmap1, ax=ax, log=True, plot_colorbar=True)
    plot_heatmap(m2, 0, 1, cmap2, ax=ax, log=False, plot_colorbar=False)
    ax.set(
        xticks=[],
        yticks=[],
        xticklabels=[],
        yticklabels=[],
        title=f"{i} epoch",
        aspect="equal",
    )
    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", dpi=300)
    plt.close(fig)
    buffer.seek(0)
    return buffer.read()


def import_data(
    path_to_tsv: pathlib.Path, resolution: int, max_epochs: Union[None, int]
) -> pd.DataFrame:
    df = pd.read_table(path_to_tsv, names=["epoch", "pos1", "pos2"]).sort_values(
        ["epoch", "pos1"]
    )

    if max_epochs is not None:
        df = df[df["epoch"] <= max_epochs]

    assert len(df) != 0

    df["bin1"] = df["pos1"] // resolution
    df["bin2"] = df["pos2"] // resolution
    df["n"] = 1
    return df


def compute_vmax(df: pd.DataFrame) -> int:
    return df.groupby(["bin1", "bin2"])["n"].sum().max()


def run_ffmpeg(
    args: List[str],
    outname: pathlib.Path,
    queue: mp.Queue,
):
    outname.parent.mkdir(exist_ok=True, parents=True)
    outname.unlink(missing_ok=True)
    with sp.Popen(
        args + [str(outname)],
        stdin=sp.PIPE,
        stdout=None,
        stderr=None,
    ) as ffmpeg:
        while True:
            async_res = queue.get(block=True, timeout=None)
            if async_res is None:
                ffmpeg.stdin.flush()
                ffmpeg.stdin.close()
                return

            ffmpeg.stdin.write(async_res.get())


def main():
    _register_cmaps()
    args = vars(make_cli().parse_args())

    fall_cmap = mpl.colormaps["fall"]
    grey_cmap = mpl.colormaps["Greys"]
    grey_cmap.set_under("k", alpha=0)

    resolution = args["resolution"]
    step = args["step"]
    outname = args["output-video"]
    if not args["force"]:
        handle_path_collisions(outname)

    df = import_data(args["tsv"], resolution, args["max_epochs"])

    start_pos, end_pos = compute_boundaries(df, resolution)
    shape = (end_pos - start_pos) // resolution
    first_epoch = df["epoch"].min()

    vmax = compute_vmax(df)
    ffmpeg_args = get_ffmpeg_args(
        shutil.which(args["ffmpeg"]),
        args["video_encoder"],
        args["fps"],
        args["crf"],
        args["encoding_threads"],
    )

    background_matrix = np.zeros([shape, shape], dtype=int)
    with mp.Manager() as manager:
        task_queue = manager.Queue(maxsize=args["nproc"] * 5)
        with manager.Pool(args["nproc"]) as pool:
            ffmpeg = pool.apply_async(run_ffmpeg, (ffmpeg_args, outname, task_queue))

            for epoch, df1 in df.groupby("epoch"):
                foreground_matrix = df_to_numpy(df1[["bin1", "bin2", "n"]], shape)
                background_matrix += foreground_matrix

                if (epoch - first_epoch) % step == 0:
                    task_queue.put(
                        pool.apply_async(
                            plot_heatmap_worker,
                            (
                                background_matrix,
                                foreground_matrix,
                                epoch,
                                vmax,
                                fall_cmap,
                                grey_cmap,
                            ),
                        )
                    )

            task_queue.put(None)
            ffmpeg.wait()


if __name__ == "__main__":
    main()
