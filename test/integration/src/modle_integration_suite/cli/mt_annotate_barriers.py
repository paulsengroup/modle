# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import pathlib
from typing import Dict, List, Tuple

import structlog
from immutabledict import ImmutableOrderedDict, immutabledict

from modle_integration_suite.tests.mt_annotate_barriers import (
    MoDLEToolsAnnotateBarriers,
    MoDLEToolsAnnotateBarriersCli,
)


def _plan_tests_cli(
    title: str = "modle-tools annotate-barriers-cli",
) -> List[ImmutableOrderedDict]:
    factory = {
        "title": title,
        "timeout": 5.0,
        "expect_failure": True,
        "validation_args": immutabledict(),
    }
    plans = (
        factory | {"args": tuple(("annotate-barriers",))},
        factory | {"args": tuple(("annotate-barriers", "--help")), "expect_failure": False},
    )

    plans = list(set(immutabledict(p) for p in plans))
    structlog.get_logger().debug(f"{title}: generated {len(plans)} test cases")
    return plans


def _plan_tests(
    data_dir: pathlib.Path,
    title: str = "modle-tools annotate-barriers",
) -> List[ImmutableOrderedDict]:
    candidate_barriers = data_dir / "ENCSR942XQI_candidate_barriers.bed.xz"
    fold_change = data_dir / "ENCSR942XQI_fc.bw"
    reference_bed = data_dir / "mt_annotate_barriers_reference_001.bed.xz"

    factory = {
        "title": title,
        "timeout": 5.0,
        "expect_failure": False,
    }

    args = [
        "annotate-barriers",
        fold_change,
        candidate_barriers,
    ]

    plans = (
        factory
        | {
            "args": tuple(args),
            "reference_bed": reference_bed,
        },
    )

    plans = list(set(immutabledict(p) for p in plans))
    structlog.get_logger().debug(f"{title}: generated {len(plans)} test cases")
    return plans


def run_tests(
    modle_tools: pathlib.Path,
    data_dir: pathlib.Path,
) -> Tuple[int, int, int, Dict]:

    num_pass = 0
    num_fail = 0
    num_skip = 0
    results = {}

    logger = structlog.get_logger().bind()

    plans = _plan_tests_cli() + _plan_tests(data_dir=data_dir)
    for p in plans:
        title = p["title"]
        assert title.startswith("modle-tools")
        if title.endswith("-cli"):
            test = MoDLEToolsAnnotateBarriersCli(modle_tools)
        else:
            test = MoDLEToolsAnnotateBarriers(modle_tools)

        status = test.run(**p)
        num_pass += status["status"] == "PASS"
        num_fail += status["status"] == "FAIL"
        results.setdefault(title, []).append(status)
        if status["status"] == "PASS":
            logger.bind(**status).info("")
        else:
            logger.bind(**status).warning("")

    return num_pass, num_fail, num_skip, results
