#!/usr/bin/env python3

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import hashlib
import os
import sys

import cooler
import numpy as np
import pandas as pd
from PIL import Image


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser()

    # Example: https://upload.wikimedia.org/wikipedia/commons/2/29/Japanese_Squirrel_edited_version.jpg
    cli.add_argument("path-to-image", type=str)
    cli.add_argument("output-cooler", type=str)

    cli.add_argument("--bin-size", type=int, default=1000)
    cli.add_argument("--chrom-name", type=str, default="test")
    cli.add_argument("--add-noise", action="store_true", default=False)
    cli.add_argument("--img-source", type=str)

    return cli


def read_image(path_to_img: str) -> np.ndarray:
    return np.array(Image.open(path_to_img).convert("L")).astype(int)


def hash_ndarray(m: np.ndarray) -> int:
    hasher = hashlib.sha256()
    hasher.update(m.tobytes())
    digest = int(hasher.hexdigest(), 16)

    return digest % (2**32 - 1)


def add_noise_to_image(img: np.ndarray) -> np.ndarray:
    lb = np.min(img)
    ub = np.max(img)
    dynamic_range = ub - lb

    np.random.seed(hash_ndarray(img))
    noise = np.random.uniform(dynamic_range * -0.10, dynamic_range * 0.10, img.shape)

    return np.clip(np.round(img + noise, decimals=3), lb, ub)


def create_bin_table(img: np.ndarray, chrom_name: str, bin_size: int) -> pd.DataFrame:
    assert bin_size > 0
    num_bins = sum(img.shape)
    return cooler.util.binnify(pd.Series({chrom_name: num_bins * bin_size}), bin_size)


def img_to_df(img) -> pd.DataFrame:
    pixels = {"bin1_id": [], "bin2_id": [], "count": []}
    offset = img.shape[0]
    for i in range(img.shape[0]):
        for j in range(img.shape[1]):
            pixels["bin1_id"].append(i)
            pixels["bin2_id"].append(offset + j)
            pixels["count"].append(img[i, j])

    return pd.DataFrame(pixels)


def main():
    args = vars(make_cli().parse_args())

    img = read_image(args["path-to-image"])
    if args["add_noise"]:
        img = add_noise_to_image(img)

    metadata = {"script": os.path.basename(sys.argv[0])} | args

    cooler.create_cooler(
        args["output-cooler"],
        create_bin_table(img, args["chrom_name"], args["bin_size"]),
        img_to_df(img),
        dtypes={"bin1_id": int, "bin2_id": int, "count": img.dtype},
        metadata=metadata,
    )


if __name__ == "__main__":
    main()
