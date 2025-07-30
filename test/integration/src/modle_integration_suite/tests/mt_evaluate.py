# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import os
import pathlib
from timeit import default_timer as timer
from typing import Any, Dict, List

import pandas as pd

from modle_integration_suite import validators
from modle_integration_suite.runners.modle import MoDLETestHarness

from .cli import MoDLECli


class MoDLEToolsEvaluateCli(MoDLECli):
    def __repr__(self) -> str:
        return "modle-tools-evaluate-cli"


class MoDLEToolsEvaluate(MoDLETestHarness):
    def __repr__(self) -> str:
        return "modle-tools-evaluate"

    def _validate(
        self,
        reference_prefix: pathlib.Path,
        test_prefix: pathlib.Path,
        expect_failure: bool,
    ):  # noqa
        if expect_failure:
            self._handle_expected_failure()
            return

        if self.returncode != 0:
            self._failures["unexpected return code"] = f"expected zero, found {self.returncode}"
            self._failures["stderr"] = self.stderr(500)
            return

        for suffix in ("horizontal", "vertical"):
            self._failures |= validators.compare_tables(
                pd.read_table(f"{reference_prefix}_{suffix}.tsv.gz"),
                pd.read_table(f"{test_prefix}_{suffix}.tsv.gz"),
            )
            self._failures |= validators.compare_bigwigs(
                pathlib.Path(f"{reference_prefix}_{suffix}.bw"),
                pathlib.Path(f"{test_prefix}_{suffix}.bw"),
            )

    def run(
        self,
        args: List[str],
        reference_prefix: pathlib.Path,
        test_prefix: pathlib.Path,
        timeout: int = 3600,
        env_variables: Dict[str, str] | None = None,
        expect_failure: bool = False,
        max_attempts: int = 1,
        title: str | None = None,
        id: str | None = None,  # noqa
    ) -> Dict[str, Any]:
        if title is None:
            title = str(self)

        self.clear()
        self._id = id
        self._title = title
        self._args = args
        self._expect_failure = expect_failure

        if env_variables is None:
            env_variables = os.environ.copy()
        else:
            env_variables = dict(env_variables.copy())

        if "LLVM_PROFILE_FILE" in env_variables:
            env_variables["LLVM_PROFILE_FILE"] = env_variables["LLVM_PROFILE_FILE"].replace("%id", str(id))

        t0 = timer()
        self._run_modle_tools(args, timeout=timeout, env_variables=env_variables, max_attempts=max_attempts)

        t1 = timer()
        self._validate(reference_prefix=reference_prefix, test_prefix=test_prefix, expect_failure=expect_failure)
        t2 = timer()

        self._modle_duration = t1 - t0
        self._validation_duration = t2 - t1
        self._duration = t2 - t0

        return self.status()
