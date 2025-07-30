# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import pathlib
import tempfile
from typing import Dict, List, Tuple

import structlog
from immutabledict import ImmutableOrderedDict, immutabledict

from modle_integration_suite.tests.mt_evaluate import (
    MoDLEToolsEvaluate,
    MoDLEToolsEvaluateCli,
)


def _plan_tests_cli(
    title: str = "modle-tools evaluate-cli",
) -> List[ImmutableOrderedDict]:
    factory = {
        "title": title,
        "timeout": 5.0,
        "expect_failure": True,
        "validation_args": immutabledict(),
    }
    plans = (
        factory | {"args": tuple(("evaluate",))},
        factory | {"args": tuple(("evaluate", "--help")), "expect_failure": False},
    )

    plans = list(set(immutabledict(p) for p in plans))
    structlog.get_logger().debug(f"{title}: generated {len(plans)} test cases")
    return plans


def _plan_tests(
    data_dir: pathlib.Path,
    out_dir: pathlib.Path,
    threads: int,
    title: str = "modle-tools evaluate",
) -> List[ImmutableOrderedDict]:
    matrix1 = data_dir / "4DNFI9GMP2J8_chr20_25kbp_mt_eval.cool"
    matrix2 = data_dir / "4DNFIFJH2524_chr20_25kbp_mt_eval.cool"

    reference_prefix = data_dir / "4DNFI9GMP2J8_vs_4DNFIFJH2524_mt_eval_custom_score"
    output_prefix = out_dir / "001" / "4DNFI9GMP2J8_vs_4DNFIFJH2524_mt_eval_custom_score"

    factory = {
        "title": title,
        "timeout": 5.0,
        "expect_failure": False,
    }

    args = [
        "evaluate",
        "--input-cooler",
        matrix1,
        "--reference-cooler",
        matrix2,
        "--metric",
        "custom",
        "--diagonal-width",
        "3mbp",
        "--threads",
        str(threads),
        "--output-prefix",
        output_prefix,
    ]

    plans = (
        factory
        | {
            "args": tuple(args),
            "reference_prefix": pathlib.Path(f"{reference_prefix}_custom_metric"),
            "test_prefix": pathlib.Path(f"{output_prefix}_custom_metric"),
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
                test = MoDLEToolsEvaluateCli(modle_tools)
            else:
                test = MoDLEToolsEvaluate(modle_tools)

            status = test.run(**p)
            num_pass += status["status"] == "PASS"
            num_fail += status["status"] == "FAIL"
            results.setdefault(title, []).append(status)
            if status["status"] == "PASS":
                logger.bind(**status).info("")
            else:
                logger.bind(**status).warning("")

    return num_pass, num_fail, num_skip, results
