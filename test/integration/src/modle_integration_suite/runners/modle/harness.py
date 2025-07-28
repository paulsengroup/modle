# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import os
import pathlib
from datetime import timedelta
from timeit import default_timer as timer
from typing import Any, Dict, List

import pandas as pd

from modle_integration_suite.runners import Runner


class MoDLETestHarness:
    def __init__(
        self,
        modle_exec: pathlib.Path | str,
        cwd: pathlib.Path | str | None = None,
        tmpdir: pathlib.Path | str | None = None,
    ):
        self._exec = pathlib.Path(modle_exec) if modle_exec else None
        self._cwd = pathlib.Path(cwd) if cwd else None
        self._tmpdir = pathlib.Path(tmpdir) if tmpdir else None

        self._id = None
        self._args = []
        self._expect_failure = None

        self._returncode = None
        self._stdout = None
        self._stderr = None

        self._failures = {}
        self._title = None

        self._modle_duration = None
        self._validation_duration = None
        self._duration = None

    def _get_modle_keyword_option(self, option_name: str, default=None) -> Any:
        try:
            i = self._args.index(option_name)
        except ValueError:
            return default

        if i + 1 == len(self._args):
            return default

        return self._args[i + 1]

    def _get_modle_flag_value(self, flag_name: str) -> bool:
        return flag_name in self._args

    def _run_modle(
        self,
        args_: List[str],
        timeout: int = 1,
        env_variables: Dict[str, str] | None = None,
        max_attempts: int = 1,
    ):
        with Runner(self._exec, args_, cwd=self._cwd, tmpdir=self._tmpdir) as runner:
            self._args = runner.args
            self._returncode, self._stdout, self._stderr = runner.run(
                timeout=timeout,
                env_variables=env_variables,
                max_attempts=max_attempts,
            )

    def _run_modle_tools(self, *args, **kwargs):
        self._run_modle(*args, **kwargs)

    def _handle_expected_failure(self):
        if len(self.stderr()) == 0:
            self._failures["missing error message"] = ""
        if len(self.stdout()) != 0:
            self._failures["unexpected output on stdout"] = self.stdout(500).strip()
        if self.returncode == 0:
            self._failures["unexpected return code"] = f"expected non-zero, found {self.returncode}"

    def _validate(self, **kwargs):
        raise NotImplementedError

    def run(
        self,
        args: List[str],
        validation_args: Dict[str, Any],
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
        self._validate(expect_failure=expect_failure, **validation_args)
        t2 = timer()

        self._modle_duration = t1 - t0
        self._validation_duration = t2 - t1
        self._duration = t2 - t0

        return self.status()

    def ok(self) -> bool:
        return self._returncode == 0 and len(self._failures) == 0

    @property
    def returncode(self) -> int:
        return self._returncode

    @property
    def duration(self) -> float:
        return self._duration

    def stderr(self, max_length: int | None = None) -> str:
        payload = "".join(self._stderr)
        if max_length is None:
            return payload

        if len(payload) > max_length:
            return f"{payload[-max_length:]}\n -- truncated"
        return payload

    def stdout(self, max_length: int | None = None) -> str:
        if isinstance(self._stdout, pd.DataFrame):
            if len(self._stdout) == 0:
                payload = ""
            else:
                payload = str(self._stdout)
        else:
            payload = "".join(self._stdout)

        if max_length is None:
            return payload

        if len(payload) > max_length:
            return f"{payload[-max_length:]}\n -- truncated"
        return payload

    @property
    def args(self) -> List[str]:
        return self._args

    @property
    def failures(self) -> Dict[str, str]:
        return self._failures

    def clear(self):
        self._id = None
        self._args = []
        self._expect_failure = None
        self._returncode = None
        self._stdout = None
        self._stderr = None
        self._failures = {}
        self._title = None
        self._modle_duration = None
        self._validation_duration = None
        self._duration = None

    def status(self) -> Dict[str, Any]:
        s = {
            "id": str(self._id),
            "title": str(self._title),
            "args": self.args[1:],
            "modle-runtime": str(timedelta(seconds=self._modle_duration)),
            "validation-runtime": str(timedelta(seconds=self._validation_duration)),
            "elapsed-time": str(timedelta(seconds=self._duration)),
            "exit-code": self._returncode,
            "expect-failure": self._expect_failure,
            "errors": [],
            "status": "PASS" if len(self._failures) == 0 else "FAIL",
        }

        if s["status"] == "PASS":
            return s

        for k, v in self._failures.items():
            if len(v) == 0:
                s["errors"].append(k)
            else:
                s["errors"].append(f"{k}: {v}")

        return s


class MoDLEToolsTestHarness(MoDLETestHarness):
    def status(self) -> Dict[str, Any]:
        return {k.replace("modle", "modle_tools"): v for k, v in super().status()}
