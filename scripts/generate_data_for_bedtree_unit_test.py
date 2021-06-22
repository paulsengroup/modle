import numpy as np
import pandas as pd

seed = 1749789967685151308
rng = np.random.default_rng(seed)

bed_features = "/home/roby/github/modle/test/data/unit_tests/H1_ctcf_all_chroms.bed"
num_samples_per_chrom = 25
chroms = ["chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY"]

# Generate ranges
ranges = []
for (i, chrom) in enumerate(chroms):
    for (center, range_) in zip(rng.uniform(0, 250e6, num_samples_per_chrom),
                                rng.normal(3e6, 500e3, num_samples_per_chrom)):
        center = int(center + int(rng.uniform() > 0.5))
        range_ = int(range_ + int(rng.uniform() > 0.5))
        ranges.append([chrom, max(0, center - range_), center + range_])


records = pd.read_csv(bed_features, sep="\t", header=None)
records.columns = ["chrom", "start", "end", "name", "score", "strand", "foo"]

# Count the number of features for each range
for (chrom, start, end) in ranges:
    num_features = ((records["chrom"] == chrom) & (records["end"] >= start) & (records["start"] < end)).sum()
    print(f"{chrom}\t{start}\t{end}\t{num_features}")

