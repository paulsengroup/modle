#!/usr/bin/env python3

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import numpy as np
import pandas as pd
import argparse

def make_cli():
  cli = argparse.ArgumentParser()
  cli.add_argument("--input-bed",
                   type=str,
                   required=True,
                   help="Path to BED file.")
  cli.add_argument("--seed",
                   type=int,
                   default=1749789967685151308)
  cli.add_argument("--samples-per-chrom",
                   type=int,
                   default=25)

  return cli

if __name__ == "__main__":
  args = make_cli().parse_args()
  rng = np.random.default_rng(args.seed)

  bed_features = args.input_bed
  num_samples_per_chrom = args.samples_per_chrom
  chroms = ["chr1", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr2",
            "chr20", "chr21", "chr22", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chrX", "chrY"]

  # Generate ranges
  ranges = []
  for (i, chrom) in enumerate(chroms):
    for (center, range_) in zip(rng.uniform(0, 250e6, num_samples_per_chrom),
                                rng.normal(3e6, 500e3, num_samples_per_chrom)):
      center = int(center + int(rng.uniform() > 0.5))
      range_ = int(range_ + int(rng.uniform() > 0.5))
      ranges.append([chrom, max(0, center - range_), center + range_])

  records = pd.read_table(bed_features,
                          usecols=range(0, 3),
                          names=["chrom", "start", "end"])

  # Count the number of features for each range
  for (chrom, start, end) in ranges:
    num_features = ((records["chrom"] == chrom) & (records["end"] >= start) & (records["start"] < end)).sum()
    print(f"{chrom}\t{start}\t{end}\t{num_features}")
