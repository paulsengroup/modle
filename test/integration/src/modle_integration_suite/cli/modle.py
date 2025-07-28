# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import pathlib
import tempfile
from typing import Dict, List, Tuple

import structlog
from immutabledict import ImmutableOrderedDict, immutabledict

from modle_integration_suite.tests.modle import MoDLESimulate, MoDLESimulateCli


def _plan_tests_cli(
    title: str = "modle-simulate-cli",
) -> List[ImmutableOrderedDict]:
    factory = {
        "title": title,
        "timeout": 5.0,
        "expect_failure": True,
        "validation_args": immutabledict(),
    }
    plans = (
        factory | {"args": tuple(("simulate",))},
        factory | {"args": tuple(("simulate", "--help")), "expect_failure": False},
    )

    plans = list(set(immutabledict(p) for p in plans))
    structlog.get_logger().debug(f"{title}: generated {len(plans)} test cases")
    return plans


def _plan_tests(
    data_dir: pathlib.Path,
    out_dir: pathlib.Path,
    threads: int,
    title: str = "modle-simulate",
) -> List[ImmutableOrderedDict]:
    reference_prefix = data_dir / "modle_sim_reference_001"
    out_prefix = out_dir / "001" / "modle_sim_001"

    chrom_sizes = data_dir / "grch38.chrom.sizes"
    regions_bed = data_dir / "grch38_regions_of_interest.bed"
    extrusion_barriers = data_dir / "grch38_h1_extrusion_barriers.bed.xz"

    factory = {
        "title": title,
        "timeout": 60.0,
        "expect_failure": False,
    }

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

    plans = (
        factory
        | {
            "args": tuple(args),
            "validation_args": immutabledict({"test_prefix": out_prefix, "reference_prefix": reference_prefix}),
        },
    )

    plans = list(set(immutabledict(p) for p in plans))
    structlog.get_logger().debug(f"{title}: generated {len(plans)} test cases")
    return plans


def run_tests(
    modle: pathlib.Path,
    data_dir: pathlib.Path,
    threads: int,
) -> Tuple[int, int, int, Dict]:

    num_pass = 0
    num_fail = 0
    num_skip = 0
    results = {}

    logger = structlog.get_logger().bind()

    with tempfile.TemporaryDirectory() as tmpdir:
        plans = _plan_tests_cli() + _plan_tests(
            data_dir=data_dir,
            out_dir=pathlib.Path(tmpdir),
            threads=threads,
        )
        for p in plans:
            title = p["title"]
            assert title.startswith("modle")
            if title.endswith("-cli"):
                test = MoDLESimulateCli(modle)
            else:
                test = MoDLESimulate(modle)

            status = test.run(**p)
            num_pass += status["status"] == "PASS"
            num_fail += status["status"] == "FAIL"
            results.setdefault(title, []).append(status)
            if status["status"] == "PASS":
                logger.bind(**status).info("")
            else:
                logger.bind(**status).warning("")

    return num_pass, num_fail, num_skip, results
