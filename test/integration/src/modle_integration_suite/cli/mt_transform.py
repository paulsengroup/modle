# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import pathlib
import tempfile
from typing import Dict, List, Tuple

import structlog
from immutabledict import ImmutableOrderedDict, immutabledict

from modle_integration_suite.tests.mt_transform import (
    MoDLEToolsTransform,
    MoDLEToolsTransformCli,
)


def _plan_tests_cli(
    title: str = "modle-tools transform-cli",
) -> List[ImmutableOrderedDict]:
    factory = {
        "title": title,
        "timeout": 5.0,
        "expect_failure": True,
        "validation_args": immutabledict(),
    }
    plans = (
        factory | {"args": tuple(("transform",))},
        factory | {"args": tuple(("transform", "--help")), "expect_failure": False},
    )

    plans = list(set(immutabledict(p) for p in plans))
    structlog.get_logger().debug(f"{title}: generated {len(plans)} test cases")
    return plans


def _plan_tests(
    data_dir: pathlib.Path,
    out_dir: pathlib.Path,
    threads: int,
    title: str = "modle-tools transform",
) -> List[ImmutableOrderedDict]:
    matrix = data_dir / "4DNFI9GMP2J8_chr20_25kbp.cool"
    reference_matrix_gauss = data_dir / "4DNFI9GMP2J8_chr20_25kbp_blurred.cool"
    reference_matrix_dog = data_dir / "4DNFI9GMP2J8_chr20_25kbp_dog.cool"

    output_matrix_gauss = out_dir / "gauss.cool"
    output_matrix_dog = out_dir / "dog.cool"

    factory = {
        "title": title,
        "timeout": 30.0,
        "expect_failure": False,
    }

    base_args = [
        "transform",
        "--input-cooler",
        matrix,
        "--diagonal-width",
        "3mbp",
        "--threads",
        str(threads),
    ]

    plans = (
        factory
        | {
            "args": tuple(base_args + ["--method", "gaussian_blur", "--output-cooler", output_matrix_gauss]),
            "reference_matrix": reference_matrix_gauss,
            "test_matrix": output_matrix_gauss,
        },
        factory
        | {
            "args": tuple(base_args + ["--method", "difference_of_gaussians", "--output-cooler", output_matrix_dog]),
            "reference_matrix": reference_matrix_dog,
            "test_matrix": output_matrix_dog,
        },
    )

    plans = list(set(immutabledict(p) for p in plans))
    structlog.get_logger().debug(f"{title}: generated {len(plans)} test cases")
    return plans


def run_tests(
    modle_tools: pathlib.Path,
    data_dir: pathlib.Path,
    threads: int,
) -> Tuple[int, int, int, Dict]:

    num_pass = 0
    num_fail = 0
    num_skip = 0
    results = {}

    logger = structlog.get_logger().bind()

    with tempfile.TemporaryDirectory() as tmpdir:
        out_dir = pathlib.Path(tmpdir)
        plans = _plan_tests_cli() + _plan_tests(data_dir=data_dir, out_dir=out_dir, threads=threads)
        for p in plans:
            title = p["title"]
            assert title.startswith("modle-tools")
            if title.endswith("-cli"):
                test = MoDLEToolsTransformCli(modle_tools)
            else:
                test = MoDLEToolsTransform(modle_tools)

            status = test.run(**p)
            num_pass += status["status"] == "PASS"
            num_fail += status["status"] == "FAIL"
            results.setdefault(title, []).append(status)
            if status["status"] == "PASS":
                logger.bind(**status).info("")
            else:
                logger.bind(**status).warning("")

    return num_pass, num_fail, num_skip, results
