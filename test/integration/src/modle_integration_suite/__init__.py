# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

from importlib.metadata import version

from . import cli, runners, tests, validators
from .main import main

__version__ = version("modle_integration_suite")
