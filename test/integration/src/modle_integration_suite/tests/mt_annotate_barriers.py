# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import os
import pathlib
from io import StringIO
from timeit import default_timer as timer
from typing import Any, Dict, List

import pandas as pd

from modle_integration_suite import validators
from modle_integration_suite.runners.modle import MoDLETestHarness

from .cli import MoDLECli


class MoDLEToolsAnnotateBarriersCli(MoDLECli):
    def __repr__(self) -> str:
        return "modle-tools-annotate-barriers-cli"


class MoDLEToolsAnnotateBarriers(MoDLETestHarness):
    def __repr__(self) -> str:
        return "modle-tools-annotate-barriers"

    def _validate(
        self,
        reference_bed: pd.DataFrame,
        test_bed: pd.DataFrame,
        expect_failure: bool,
    ):  # noqa
        if expect_failure:
            self._handle_expected_failure()
            return

        if self.returncode != 0:
            self._failures["unexpected return code"] = f"expected zero, found {self.returncode}"
            self._failures["stderr"] = self.stderr(500)
            return

        self._failures = validators.compare_bed6(reference_bed, test_bed)

    def run(
        self,
        args: List[str],
        reference_bed: pathlib.Path,
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
        self._run_modle(args, timeout=timeout, env_variables=env_variables, max_attempts=max_attempts)

        t1 = timer()
        cols = ["chrom", "start", "end", "name", "score", "strand"]
        df1 = pd.read_table(reference_bed, names=cols)
        df2 = pd.read_table(StringIO(self.stdout()), names=cols)
        self._validate(reference_bed=df1, test_bed=df2, expect_failure=expect_failure)
        t2 = timer()

        self._modle_duration = t1 - t0
        self._validation_duration = t2 - t1
        self._duration = t2 - t0

        return self.status()
