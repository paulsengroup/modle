# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import pathlib
import tempfile
from typing import Dict, Tuple

import structlog

from modle_integration_suite.tests.modle import MoDLESimulate, MoDLESimulateCli


def _run_test_001(modle: pathlib.Path, data_dir: pathlib.Path, threads: int):
    with tempfile.TemporaryDirectory() as tmpdir:
        outdir = pathlib.Path(tmpdir)
        reference_prefix = data_dir / "modle_sim_reference_001"
        out_prefix = outdir / "modle_sim_001"

        chrom_sizes = data_dir / "grch38.chrom.sizes"
        regions_bed = data_dir / "grch38_regions_of_interest.bed"
        extrusion_barriers = data_dir / "grch38_h1_extrusion_barriers.bed.xz"

        args = [
            "simulate",
            "--chrom-sizes",
            chrom_sizes,
            "--genomic-intervals",
            regions_bed,
            "--extrusion-barrier-file",
            extrusion_barriers,
            "--output-prefix",
            out_prefix,
            "--resolution",
            "20kb",
            "--verbose",
            "--target-contact-density",
            "20",
            "--ncells",
            "2",
            "--track-1d-lef-position",
            "--max-burnin-epochs",
            "5000",
            "--threads",
            str(threads),
        ]

        return MoDLESimulate(modle).run(
            args=args,
            validation_args={"test_prefix": out_prefix, "reference_prefix": reference_prefix},
        )


def run_tests(
    modle: pathlib.Path,
    data_dir: pathlib.Path,
    threads: int,
) -> Tuple[int, int, int, Dict]:

    logger = structlog.get_logger().bind()
    results = {}

    status = _run_test_001(modle, data_dir, threads)
    num_pass = status["status"] == "PASS"
    num_fail = status["status"] == "FAIL"
    num_skip = 0

    results.setdefault("modle-simulate", []).append(status)
    if status["status"] == "PASS":
        logger.bind(**status).info("")
    else:
        logger.bind(**status).warning("")

    return num_pass, num_fail, num_skip, results
