<!--
Copyright (C) 2025 Roberto Rossini <roberros@uio.no>

SPDX-License-Identifier: MIT
-->

# README

Installing the pipeline:

```bash
python3 -m venv venv
venv/bin/pip install . -v
```

Running the tests:

```console
user@dev:/tmp$ venv/bin/modle_integration_suite --modle-bin /tmp/modle/build/src/modle/modle --modle-tools-bin /tmp/modle/build/src/modle_tools/modle_tools --data-dir ../data/ --threads="$(nproc)"

2025-07-28 16:11:24.456902 [info     ] [simulate (CLI)          ] [id=None    ] [runtime=31.989ms  ] PASS
2025-07-28 16:11:24.489230 [info     ] [simulate (CLI)          ] [id=None    ] [runtime=31.948ms  ] PASS
2025-07-28 16:12:26.396835 [info     ] [simulate                ] [id=None    ] [runtime=1m:1.907s ] PASS
2025-07-28 16:12:26.397082 [info     ] running tests for modle took 61.97s
2025-07-28 16:12:26.429557 [info     ] [tools annotate-barriers (CLI)] [id=None    ] [runtime=32.291ms  ] PASS
2025-07-28 16:12:26.461938 [info     ] [tools annotate-barriers (CLI)] [id=None    ] [runtime=32.169ms  ] PASS
2025-07-28 16:12:26.529217 [info     ] [tools annotate-barriers ] [id=None    ] [runtime=67.120ms  ] PASS
2025-07-28 16:12:26.529279 [info     ] running tests for modle_tools annotate-barriers took 0.13s
2025-07-28 16:12:26.561381 [info     ] [tools evaluate (CLI)    ] [id=None    ] [runtime=31.901ms  ] PASS
2025-07-28 16:12:26.593410 [info     ] [tools evaluate (CLI)    ] [id=None    ] [runtime=31.884ms  ] PASS
2025-07-28 16:12:27.065743 [info     ] [tools evaluate          ] [id=None    ] [runtime=472.197ms ] PASS
2025-07-28 16:12:27.065889 [info     ] running tests for modle_tools evaluate took 0.54s
2025-07-28 16:12:27.098014 [info     ] [tools transform (CLI)   ] [id=None    ] [runtime=31.947ms  ] PASS
2025-07-28 16:12:27.130176 [info     ] [tools transform (CLI)   ] [id=None    ] [runtime=32.011ms  ] PASS
2025-07-28 16:12:31.987124 [info     ] [tools transform         ] [id=None    ] [runtime=4.857s    ] PASS
2025-07-28 16:12:35.719145 [info     ] [tools transform         ] [id=None    ] [runtime=3.732s    ] PASS
2025-07-28 16:12:35.719707 [info     ] running tests for modle_tools transform took 8.65s
2025-07-28 16:12:35.719748 [info     ] running 13 tests took 71.32s

# PASS: 13
# SKIP: 0
# FAIL: 0
```
