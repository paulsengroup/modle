#!/usr/bin/env python3

# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import datetime
import json
import multiprocessing as mp
import pathlib
import platform
import sys
import time
from typing import Dict, Tuple

import click
import structlog

import modle_integration_suite.cli.modle as modle
import modle_integration_suite.cli.mt_annotate_barriers as mt_annotate_barriers
import modle_integration_suite.cli.mt_evaluate as mt_evaluate
import modle_integration_suite.cli.mt_transform as mt_transform
from modle_integration_suite.cli.logging import setup_logger
from modle_integration_suite.runners.modle.common import version


def nproc() -> int:
    return mp.cpu_count()


def init_results(modle_bin: pathlib.Path) -> Dict:
    return {
        "platform": platform.platform(),
        "arch": platform.machine(),
        "modle-version": version(modle_bin),
        "date": datetime.datetime.now().isoformat(),
        "results": {
            "pass": 0,
            "fail": 0,
            "skip": 0,
        },
    }


def run_tests_modle(
    modle_bin: pathlib.Path | None,
    data_dir: pathlib.Path,
    threads: int,
) -> Tuple[int, int, int, Dict]:
    logger = structlog.get_logger().bind()
    if modle_bin is None:
        logger.info("skipping tests for modle as no binary was provided!")
        return 0, 0, 0, {}

    t0 = time.time()
    res = modle.run_tests(
        modle=modle_bin,
        data_dir=data_dir,
        threads=threads,
    )
    delta = time.time() - t0
    logger.info(f"running tests for modle took {delta:.2f}s")
    return res


def run_tests_mt_annotate_barriers(
    modle_tools_bin: pathlib.Path,
    data_dir: pathlib.Path,
) -> Tuple[int, int, int, Dict]:
    logger = structlog.get_logger().bind()

    t0 = time.time()
    res = mt_annotate_barriers.run_tests(
        modle_tools=modle_tools_bin,
        data_dir=data_dir,
    )
    delta = time.time() - t0
    logger.info(f"running tests for modle_tools annotate-barriers took {delta:.2f}s")
    return res


def run_tests_mt_evaluate(
    modle_tools_bin: pathlib.Path,
    data_dir: pathlib.Path,
    threads: int,
) -> Tuple[int, int, int, Dict]:
    logger = structlog.get_logger().bind()

    t0 = time.time()
    res = mt_evaluate.run_tests(
        modle_tools=modle_tools_bin,
        data_dir=data_dir,
        threads=threads,
    )
    delta = time.time() - t0
    logger.info(f"running tests for modle_tools evaluate took {delta:.2f}s")
    return res


def run_tests_mt_transform(
    modle_tools_bin: pathlib.Path,
    data_dir: pathlib.Path,
    threads: int,
) -> Tuple[int, int, int, Dict]:
    logger = structlog.get_logger().bind()

    t0 = time.time()
    res = mt_transform.run_tests(
        modle_tools=modle_tools_bin,
        data_dir=data_dir,
        threads=threads,
    )
    delta = time.time() - t0
    logger.info(f"running tests for modle_tools transform took {delta:.2f}s")
    return res


def run_tests_modle_tools(
    modle_tools_bin: pathlib.Path | None,
    data_dir: pathlib.Path,
    threads: int,
) -> Tuple[int, int, int, Dict]:
    logger = structlog.get_logger().bind()
    if modle_tools_bin is None:
        logger.info("skipping tests for modle_tools as no binary was provided!")
        return 0, 0, 0, {}
    num_pass, num_fail, num_skip, results = run_tests_mt_annotate_barriers(modle_tools_bin, data_dir)

    res = run_tests_mt_evaluate(modle_tools_bin, data_dir, threads)
    num_pass += res[0]
    num_fail += res[1]
    num_skip += res[2]
    results |= res[3]

    res = run_tests_mt_transform(modle_tools_bin, data_dir, threads)
    num_pass += res[0]
    num_fail += res[1]
    num_skip += res[2]
    results |= res[3]

    return num_pass, num_fail, num_skip, results


def run_tests(
    modle_bin: pathlib.Path | None,
    modle_tools_bin: pathlib.Path | None,
    data_dir: pathlib.Path,
    threads: int,
) -> Dict:
    assert modle_bin is not None or modle_tools_bin is not None

    num_pass = 0
    num_fail = 0
    num_skip = 0

    results = init_results(modle_bin if modle_bin is not None else modle_tools_bin)

    res = run_tests_modle(modle_bin, data_dir, threads)
    num_pass += res[0]
    num_fail += res[1]
    num_skip += res[2]
    results["results"] |= res[3]

    res = run_tests_modle_tools(modle_tools_bin, data_dir, threads)
    num_pass += res[0]
    num_fail += res[1]
    num_skip += res[2]
    results["results"] |= res[3]

    results["results"]["pass"] = num_pass
    results["results"]["fail"] = num_fail
    results["results"]["skip"] = num_skip
    results["success"] = num_fail == 0

    return results


@click.command()
@click.option(
    "--modle-bin",
    type=click.Path(
        exists=True,
        executable=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        path_type=pathlib.Path,
    ),
)
@click.option(
    "--modle-tools-bin",
    type=click.Path(
        exists=True,
        executable=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        path_type=pathlib.Path,
    ),
)
@click.option(
    "--data-dir",
    type=click.Path(exists=True, file_okay=False, dir_okay=True, path_type=pathlib.Path),
    help="Path to the folder with the test files.",
)
@click.option(
    "--threads",
    help="Specify the maximum number of CPU threads to be used.",
    type=click.IntRange(2, nproc()),
    default=2,
    show_default=True,
)
@click.option(
    "--verbosity",
    type=click.Choice(["debug", "info", "warning", "error", "critical"]),
    default="info",
    help="Set verbosity level.",
)
@click.option(
    "--result-file",
    help="Path where to write the test results.",
    type=pathlib.Path,
)
@click.option(
    "--force",
    help="Force overwrite existing output file(s).",
    default=False,
    is_flag=True,
    show_default=True,
)
def main(
    modle_bin: pathlib.Path,
    modle_tools_bin: pathlib.Path,
    data_dir: pathlib.Path,
    threads: int,
    verbosity: str,
    result_file: pathlib.Path,
    force: bool,
):
    """
    Run MoDLE and MoDLE-tools integration test suite.
    """
    if modle_bin is None and modle_tools_bin is None:
        raise click.UsageError("Please specify at least one of --modle-bin, --modle-tools-bin")

    setup_logger(verbosity)

    logger = structlog.get_logger()

    if result_file and result_file.exists():
        if force:
            result_file.unlink()
        else:
            raise RuntimeError(f'refusing to overwrite file "{result_file}"')

    t0 = time.time()
    results = run_tests(
        modle_bin,
        modle_tools_bin,
        data_dir,
        threads,
    )
    t1 = time.time()

    unexpected_exit_codes = {}
    num_unexpected_exit_code = 0
    for k, runs in results["results"].items():
        if not k.startswith("modle"):
            continue
        for res in runs:
            ec = res["exit-code"]
            if ec not in {0, 1}:
                if ec in unexpected_exit_codes:
                    unexpected_exit_codes[ec].append(res)
                else:
                    unexpected_exit_codes[ec] = [res]
                num_unexpected_exit_code += 1

    num_pass = results["results"]["pass"]
    num_fail = results["results"]["fail"]
    num_skip = results["results"]["skip"]

    delta = t1 - t0
    logger.info(f"running {num_pass + num_fail} tests took {delta:.2f}s")

    if result_file is not None:
        with open(result_file, "w") as f:
            f.write(json.dumps(results, indent=2))

    if num_unexpected_exit_code != 0:
        logger.warn(
            "some of the tests returned non-zero exit codes with unexpected values: "
            f"{', '.join(str(x) for x in sorted(unexpected_exit_codes.keys()))}. "
            "Please carefully review the test report."
        )
        for results in unexpected_exit_codes.values():
            for res in results:
                logger.warn(f"FAIL: {res}")

    print("", file=sys.stderr)
    print(f"# PASS: {num_pass}", file=sys.stderr)
    print(f"# SKIP: {num_skip}", file=sys.stderr)
    print(f"# FAIL: {num_fail}", file=sys.stderr)
    if num_unexpected_exit_code != 0:
        print(f"# UNEXPECTED EXIT CODE: {num_unexpected_exit_code}", file=sys.stderr)

    sys.exit(num_fail != 0 or num_unexpected_exit_code != 0)


if __name__ == "__main__":
    main()
