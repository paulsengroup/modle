#!/usr/bin/env bash

set -e

cd "$(readlink -f "$(dirname "$0")")"
cp -r "dummy_data/" /tmp || :
tmpfile=$(mktemp /tmp/jl.XXXXXX)

cat > "$tmpfile" <<-EOF
using Pkg
Pkg.instantiate()
Pkg.precompile()
cd("..")

using PackageCompiler
create_app("compute_statistics_for_bwigs", "compute_statistics_for_bwigs_app", force=true, filter_stdlibs=true, precompile_execution_file="compute_statistics_for_bwigs/precompile_app.jl")
exit()
EOF

julia --project=. -q "$tmpfile"

rm -f "$tmpfile"
