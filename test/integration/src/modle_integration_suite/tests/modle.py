# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import pathlib

from modle_integration_suite import validators
from modle_integration_suite.runners.modle import MoDLETestHarness

from .cli import MoDLECli


class MoDLESimulateCli(MoDLECli):
    def __repr__(self) -> str:
        return "modle-simulate-cli"


class MoDLESimulate(MoDLETestHarness):
    def __repr__(self) -> str:
        return "modle-simulate"

    def _validate(
        self,
        test_prefix: pathlib.Path,
        reference_prefix: pathlib.Path,
        expect_failure: bool,
    ):  # noqa
        if expect_failure:
            self._handle_expected_failure()
            return

        if len(self.stdout()) != 0:
            self._failures["unexpected output on stdout"] = self.stdout(500).strip()
        if self.returncode != 0:
            self._failures["unexpected return code"] = f"expected zero, found {self.returncode}"
            self._failures["stderr"] = self.stderr(500)
            return

        self._failures |= validators.compare_coolers(
            reference_prefix.with_suffix(".cool"),
            test_prefix.with_suffix(".cool"),
        )
        self._failures |= validators.compare_bigwigs(
            reference_prefix.with_suffix(".bw"),
            pathlib.Path(f"{test_prefix}_lef_1d_occupancy.bw"),
        )

        if not test_prefix.with_suffix(".log").is_file():
            self._failures["log file is missing"] = ""

        if not pathlib.Path(f"{test_prefix}_config.toml").is_file():
            self._failures["config file is missing"] = ""
