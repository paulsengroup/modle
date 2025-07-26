#!/usr/bin/env bash

# Copyright (c) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -eu
set -o pipefail

if [ $# -ne 1 ]; then
  1>&2 echo "Usage: $0 modle:latest"
  exit 1
fi

IMG="$1"

tmpdir="$(mktemp -d)"
# shellcheck disable=SC2064
trap "rm -rf '$tmpdir'" EXIT

cat > "$tmpdir/runme.sh" <<- 'EOM'

set -eu

whereis -b modle
whereis -b modle_tools

modle --version
modle_tools --version
EOM

chmod 755 "$tmpdir/runme.sh"

sudo docker run --rm --entrypoint=/bin/bash \
  --volume "$tmpdir/runme.sh:/tmp/runme.sh:ro" \
  "$IMG" \
  /tmp/runme.sh
