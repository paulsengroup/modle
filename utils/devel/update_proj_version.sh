#!/usr/bin/env bash

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -u

### Helper script to bump MoDLE's version prior to making a new release

argc=$#

if [ $argc -ne 5 ]; then
  >&2 echo "Usage:   $0 major minor patch suffix path/to/modle/repo/"
  >&2 echo "Example: $0 1 0 0 rc.1 modle/ -> modle-v1.0.0-rc.1"
  >&2 echo "Example: $0 1 0 0 \"\" modle/   -> modle-v1.0.0"

  exit 1
fi

function check_is_integer() {
  n="$1"
  re='^[0-9]+$'
  if ! [[ $n =~ $re ]] ; then
    >&2 echo "error: $n does not look like a valid integer"
    exit 1
  fi
}

major="$1"
minor="$2"
patch="$3"
suffix="$4"

check_is_integer "$major"
check_is_integer "$minor"
check_is_integer "$patch"

version="$major.$minor.$patch"
if [ -n "$suffix" ]; then
  version+="-$suffix"
else
  suffix='""'
fi

repo="$5"

# Update Versioning.cmake
sed -i "0,/set(MODLE_PROJECT_VERSION_MAJOR.*/s//set(MODLE_PROJECT_VERSION_MAJOR $major)/" "$repo/cmake/Versioning.cmake"
sed -i "0,/set(MODLE_PROJECT_VERSION_MINOR.*/s//set(MODLE_PROJECT_VERSION_MINOR $minor)/" "$repo/cmake/Versioning.cmake"
sed -i "0,/set(MODLE_PROJECT_VERSION_PATCH.*/s//set(MODLE_PROJECT_VERSION_PATCH $patch)/" "$repo/cmake/Versioning.cmake"
sed -i "0,/set(MODLE_PROJECT_VERSION_SUFFIX.*/s//set(MODLE_PROJECT_VERSION_SUFFIX $suffix)/" "$repo/cmake/Versioning.cmake"

# Update conanfile.py
# sed -i "0,/  version.*/s//  version = \"$version\"/" "$repo/conanfile.py"
