# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import pathlib
import subprocess as sp


def version(modle: pathlib.Path) -> str:
    ver = sp.check_output([modle, "--version"], encoding="utf-8").strip()
    if ver.startswith("MoDLE-tools"):
        return ver.removeprefix("MoDLE-tools-")

    if ver.startswith("MoDLE"):
        return ver.removeprefix("MoDLE-")

    raise RuntimeError(f'"{modle}" does not seem to be a modle or modle_tools binary')
