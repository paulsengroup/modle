<!--
Copyright (C) 2022 Roberto Rossini <roberros@uio.no>

SPDX-License-Identifier: MIT
-->

# Using git hooks for this repo

- Install clang-format and cmake-format
- Copy the content of this folder to `.git/hooks/` (path is relative to the root of this repo)
- Make sure scripts inside `.git/hooks/*` are executable

Next time you commit something to this repo, git will run clang-format and cmake-format
on the files you are about to commit.
