# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

from modle_integration_suite.runners.modle import MoDLETestHarness


class MoDLECli(MoDLETestHarness):
    def __repr__(self) -> str:
        return "modle-cli"

    def _validate(self, expect_failure: bool):
        if expect_failure:
            if len(self.stderr()) == 0:
                self._failures["missing help message"] = ""
            if len(self.stdout()) != 0:
                self._failures["unexpected output on stdout"] = self.stdout(500).strip()
            if self.returncode == 0:
                self._failures["unexpected return code"] = f"expected non-zero, found {self.returncode}"
            return

        if len(self.stdout()) == 0:
            self._failures["missing help message"] = ""
        if len(self.stderr()) != 0:
            self._failures["unexpected output on stderr"] = self.stderr(500).strip()
        elif not expect_failure and self.returncode != 0:
            self._failures["unexpected return code"] = f"expected zero, found {self.returncode}"


class MoDLEToolsCli(MoDLECli):
    def __repr__(self) -> str:
        return "modle_tools-cli"
