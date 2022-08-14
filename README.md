<!--
Copyright (C) 2022 Roberto Rossini <roberros@uio.no>

SPDX-License-Identifier: MIT
-->

# MoDLE

[![Unit tests Ubuntu](https://github.com/paulsengroup/modle/actions/workflows/unit-tests-ubuntu.yml/badge.svg?branch=main)](https://github.com/paulsengroup/modle/actions/workflows/unit-tests-ubuntu.yml)
[![Unit tests macOS](https://github.com/paulsengroup/modle/actions/workflows/unit-tests-macos.yml/badge.svg)](https://github.com/paulsengroup/modle/actions/workflows/unit-tests-macos.yml)
[![Build Docker image](https://github.com/paulsengroup/modle/actions/workflows/build-docker-image.yml/badge.svg)](https://github.com/paulsengroup/modle/actions/workflows/build-docker-image.yml)

## Using MoDLE

The recommended way to run MoDLE is using the Docker images hosted
on [ghcr.io](https://github.com/paulsengroup/modle/pkgs/container/modle)
and [dockerhub](https://hub.docker.com/repository/docker/paulsengroup/modle).

```bash
# Using Docker
sudo docker run ghcr.io/paulsengroup/modle:1.0.0-rc.7 --help

# Using Singularity/Apptainer
singularity run docker://ghcr.io/paulsengroup/modle:1.0.0-rc.7 --help
```

## Building MoDLE

### Pre-requisites

Building MoDLE requires a compiler toolchain supporting C++17, such as:

- GCC 8 and newer
- Clang 8 and newer
- Apple-Clang 10.0 and newer

MoDLE is developed on a Linux machine and should run on most UNIX-like OSes (including macOS).

It is in theory possible to compile and run MoDLE on Windows, but this is not officially supported.
If you really need to run MoDLE on Windows, please consider using the Docker images hosted
on [ghcr.io](https://github.com/paulsengroup/modle/pkgs/container/modle)
and [dockerhub](https://hub.docker.com/repository/docker/paulsengroup/modle).

In addition to a C++17 compiler, building MoDLE requires the following tools:

- CMake >= 3.18
- Conan >= 1.50

#### Installing Conan

Conan is a package manager for C and C++ applications, and can be installed using pip or Homebrew:

- `pip3 install "conan>=1.50"`
- `brew install conan`

### Getting MoDLE source code

We highly recommend users to download MoDLE's source code for the latest stable release from
the [Release](https://github.com/paulsengroup/modle/releases) page.

Cloning the repository is only recommended if you intend to test features/bugfixes that have not yet landed in a release.

### Compiling MoDLE

Run the following command from inside the folder where MoDLE's source code was extracted.

Here we assume the machine where MoDLE will be compiled has 8 CPU cores.
Feel free to adjust `-j 8` to match the number of CPU cores available on your machine to improve compilation speed.

```bash
mkdir build/
cd build/

# Configure project
cmake -DCMAKE_BUILD_TYPE=Release ..

# Compile project
cmake --build . -j 8
```

<details>
<summary>Notes</summary>

By default, running the commands listed in
section [Installing MoDLE](https://github.com/paulsengroup/modle#installing-modle) will install MoDLE
under `/usr/local/` (i.e. the actual binary will be located at `/usr/local/bin/modle`).

Pass `-DCMAKE_INSTALL_PREFIX="$HOME/.local/"` to the first CMake command (before `..`) to install MoDLE for your user only. In this case MoDLE binary will be located at `~/.local/bin/modle`

The path passed to CMake through `-DCMAKE_INSTALL_PREFIX` can be in principle any path where your user has write permissions.
</details>

<details>
<summary>Troubleshooting common errors</summary>

#### Incorrect or incomplete Conan profile

This will cause CMake to exit with an error during project configuration.

When this is the case, the error message should look similar to the following:

```
ERROR: libBigWig/0.4.6: 'settings.compiler' value not defined
CMake Error at build/conan.cmake:631 (message):
  Conan install failed='1'
```

This issue is usually fixed by forcing Conan to re-detect compiler information:

```bash
# Backup old profile
mv ~/.conan/profiles/default ~/.conan/profiles/default.bak

# Write the new profile
conan profile new ~/.conan/profiles/default --detect
```

If after running the previous command you see a warning mentioning `GCC OLD ABI COMPATIBILITY`, run:

```bash
conan profile update settings.compiler.libcxx=libstdc++11 default
```

On a Linux x86_64 machine with GCC 11 installed, the default profile should look similar to:

```
[settings]
os=Linux
os_build=Linux
arch=x86_64
arch_build=x86_64
compiler=gcc
compiler.version=11
compiler.libcxx=libstdc++11
build_type=Release
[options]
[build_requires]
[env]
```

On a Mac with Intel CPU the profile will be like:

```
[settings]
os=Macos
os_build=Macos
arch=x86_64
arch_build=x86_64
compiler=apple-clang
compiler.version=11.0
compiler.libcxx=libc++
build_type=Release
[options]
[build_requires]
[env]
```

Now remove the content of the build folder with e.g. `rm -r build/*` and re-run the steps listed in the [Compiling MoDLE](https://github.com/paulsengroup/modle#compiling-modle) section.

#### Need more help?
If the above troubleshooting steps did not help, feel free to get in touch by starting a new [discussion](https://github.com/paulsengroup/modle/discussions/new).
</details>


### Running automated tests

To ensure that the compiled code works as intended we highly recommend running MoDLE's automated test.

The integration test depends on `h5diff`, `shasum` and `xz`, which can be installed as follows:

```bash
# Instructions for Ubuntu
apt-get install -y hdf5-tools libdigest-sha-perl xz-utils

# Instructions for macOS
brew install hdf5 xz
```

Once the test dependencies are satisfied, the test suite can be started by running the following commands from the repository root:

```bash
cd build/

# Run unit tests
ctest -j 8                 \
      --test-dir .         \
      --schedule-random    \
      --output-on-failure  \
      --no-tests=error     \
      -E '(SciPy)|(wCorr)'

# Run integration test
../test/scripts/modle_integration_test.sh src/modle/modle

../test/scripts/modle_tools_transform_integration_test.sh src/modle_tools/modle_tools
../test/scripts/modle_tools_eval_integration_test.sh src/modle_tools/modle_tools
```

<details>
<summary>Example output</summary>

The first command should produce an output similar to the following:
```
101/110 Test #110: Generate LEF moves 001 - LONG ....................................................................   Passed    6.13 sec
        Start  69: Detect LEF-BAR collisions 001 - wo soft collisions fwd CTCFs - SHORT
102/110 Test  #69: Detect LEF-BAR collisions 001 - wo soft collisions fwd CTCFs - SHORT .............................   Passed    0.02 sec
        Start  96: Variance - SHORT
103/110 Test  #96: Variance - SHORT .................................................................................   Passed    0.01 sec
104/110 Test  #31: Writer lzma - SHORT ..............................................................................   Passed    9.40 sec
105/110 Test  #24: Reader lzma - SHORT ..............................................................................   Passed    7.46 sec
106/110 Test  #19: Reader plain - SHORT .............................................................................   Passed   14.62 sec
107/110 Test  #23: Reader lz4 - SHORT ...............................................................................   Passed    7.17 sec
108/110 Test  #30: Writer bzip2 - SHORT .............................................................................   Passed    7.67 sec
109/110 Test  #28: Writer plain - SHORT .............................................................................   Passed    6.86 sec
110/110 Test  #20: Reader plain sv - SHORT ..........................................................................   Passed   14.04 sec

100% tests passed, 0 tests failed out of 110

Total Test time (real) =  18.45 sec
```

While the output of the second command should look something like this.
```
[2022-06-15 13:28:02.649] [info]: Simulation of "chr2" successfully completed.
[2022-06-15 13:28:02.869] [info]: Writing contacts for "chr2" to file "/tmp/ci-OdNlvn6LME/out.cool"...
[2022-06-15 13:28:02.909] [info]: Written 1816500 contacts for "chr2" across 0.21M out of 1.82M pixels to file "/tmp/ci-OdNlvn6LME/out.cool".
[2022-06-15 13:28:02.909] [info]: Writing contacts for "chr20" to file "/tmp/ci-OdNlvn6LME/out.cool"...
[2022-06-15 13:28:02.909] [info]: Written 483450 contacts for "chr20" across 0.05M out of 0.48M pixels to file "/tmp/ci-OdNlvn6LME/out.cool".
[2022-06-15 13:28:02.909] [info]: Writing contacts for "chr21" to file "/tmp/ci-OdNlvn6LME/out.cool"...
[2022-06-15 13:28:02.909] [info]: Written 350400 contacts for "chr21" across 0.04M out of 0.35M pixels to file "/tmp/ci-OdNlvn6LME/out.cool".
[2022-06-15 13:28:02.909] [info]: Writing contacts for "chr22" to file "/tmp/ci-OdNlvn6LME/out.cool"...
[2022-06-15 13:28:02.909] [info]: Written 381150 contacts for "chr22" across 0.04M out of 0.38M pixels to file "/tmp/ci-OdNlvn6LME/out.cool".
[2022-06-15 13:28:03.279] [info]: Simulation terminated without errors in 4.259878566s!

Bye.
Comparing /tmp/modle-6n3WSvOXxQ/out.cool with /tmp/modle/test/data/integration_tests/reference_001.cool...

### PASS ###
/tmp/modle-6n3WSvOXxQ/out_lef_1d_occupancy.bw: OK
```

If the second test reports one or more differences between `out.cool` and `reference_001.cool`, then the test failed.

Test failure example:
```
Comparing /tmp/modle-6n3WSvOXxQ/out.cool with /home/roby/github/modle/test/data/integration_tests/reference_001.cool...

dataset: </indexes/bin1_offset> and </indexes/bin1_offset>
20154 differences found
Not comparable: </pixels/bin1_id> has rank 1, dimensions [355352], max dimensions [18446744073709551615]
and </pixels/bin1_id> has rank 1, dimensions [356001], max dimensions [18446744073709551615]
Not comparable: </pixels/bin2_id> has rank 1, dimensions [355352], max dimensions [18446744073709551615]
and </pixels/bin2_id> has rank 1, dimensions [356001], max dimensions [18446744073709551615]
Not comparable: </pixels/count> has rank 1, dimensions [355352], max dimensions [18446744073709551615]
and </pixels/count> has rank 1, dimensions [356001], max dimensions [18446744073709551615]

### FAIL ###
```

</details>

<details>
<summary>For developers</summary>
To run the full test suite, remove `-E '(SciPy)|(wCorr)` from the above snippet.

Some of MoDLE's unit tests depend on the following libraries:

- [SciPy](https://scipy.org/)
- [wCorr](https://cran.r-project.org/web/packages/wCorr/index.html)

These libraries can be installed as follows:

```bash
python3 -m pip install scipy
Rscript --no-save -e 'install.packages("wCorr", dependencies=c("Depends", "Imports", "LinkingTo"), repos="https://cloud.r-project.org")'
```

</details>

### Installing MoDLE

The following command will install MoDLE files under the prefix specified through `-DCMAKE_INSTALL_PREFIX` (`/usr/local`
by default).

```bash
# Run from the repository root
cmake --install build/
```

### Running MoDLE

Commands in this section assume you are running MoDLE using Singularity/Apptainer from the root of this repository.

Test datasets are hosted on Zenodo [10.5281/zenodo.6625788](https://doi.org/10.5281/zenodo.6625788) and can be
downloaded as follows:

```bash
# IMPORTANT! You should be in the repository root when running the following command (otherwise test files will be extracted in the wrong place)
curl -L 'https://zenodo.org/record/6638906/files/modle_test_data.tar.gz?download=1' | tar -xzf -
```

Datasets are automatically downloaded by CMake when running steps from inside folder `test/data/integration_tests`.

<details>
<summary>If you are not using containers...</summary>
If you are building MoDLE and have followed the <a href="https://github.com/paulsengroup/modle#compiling-modle">instructions</a> for compiling MoDLE, then test datasets have already been downloaded and extracted by CMake, so you can skip the above step.
</details>

#### Required input files

Running a simulation with default settings only requires two input files:

- A chrom.sizes with the list of chromosome to be simulated
- A BED file with the list of extrusion barriers to use in the simulation

The extrusion barrier BED file should have at least the first 6 columns defined (i.e. chrom, chromStart, chromEnd, name,
score and strand).

The ___name___ field is ignored, and the ___score___ field is optional (and should be set to 0 when not used).

When ___score___ is non-zero, its value will be used to set the occupancy for the extrusion barrier defined by the
current line.

The ___strand___ field is required and is used to define the extrusion barrier direction.
As of `v1.0.0-rc.7`, this field should be populated with the direction of the corresponding CTCF binding site.
Barriers without strand information (i.e. with strand '.') will be ignored.

Sample chrom.sizes and BED file(s) are available inside folder `test/data/integration_tests`.

#### Running a simulation with default settings

```bash
singularity run docker://ghcr.io/paulsengroup/modle:1.0.0-rc.7   \
    simulate \
    --chrom-sizes test/data/integration_tests/grch38.chrom.sizes \
    --extrusion-barrier-file test/data/integration_tests/grch38_h1_extrusion_barriers.bed.xz \
    --output-prefix path/to/output/prefix
```

This will create folder `path/to/output` (if it doesn't already exist), and write the following files inside it:

```
path/to/output
├── prefix_config.toml
├── prefix.cool
└── prefix.log
```

Contacts are stored in `path/to/output/prefix.cool`.

In case any of the output files already exist, MoDLE will refuse to run and print an error message listing the file name
collisions.

Passing the `--force` flag overrides this behavior and will cause MoDLE to overwrite existing files.

#### Running a simulation using config files

File `prefix_config.toml` from the previous section is a config file that can be used to re-run a simulation using the
same parameters.

```bash
# Run a simulation with the same parameter as the previous example
singularity run docker://ghcr.io/paulsengroup/modle:1.0.0-rc.7 \
    --config path/to/output/prefix_config.toml
```

Parameters read from a config file have lower priority than parameter specified through the CLI,
meaning that CLI options can be used to override parameters read from the config file.

This can be useful when running a batch of simulations using a large number of optional parameters, and where only a
handful of parameters need to change across simulation runs.

```bash
# The following command will run a simulation using parameters from the previous example as starting point,
# but using a custom lef density and overriding the output prefix specified by the config file.
singularity run docker://ghcr.io/paulsengroup/modle:1.0.0-rc.7 \
    simulate \
    --config path/to/output/prefix_config.toml \
    --lef-density 15 \
    --output-prefix path/to/output/prefix2
```

The config file is a text file in TOML format.

Adding a line like `my-option=my_value` to the config file is equivalent to passing `--my-option=my_value` on the CLI.

For an up-to-date list of supported CLI options, please refer to MoDLE's help message:

```bash
singularity run docker://ghcr.io/paulsengroup/modle:1.0.0-rc.7 simulate --help
```

<details>
<summary>Tips and tricks</summary>
<b>Compressed input files</b>

MoDLE automatically detects and handles compressed input files.

As of `v1.0.0-rc.6`, the following compression algorithms are supported:

- bzip2
- gzip
- LZ4
- LZO
- XZ/LZMA
- ZSTD

<b>Visualizing simulation result</b>

To quickly visualize .cool files we recommend using [cooler](https://github.com/open2c/cooler) show.

Example:

```bash
# Visualize a region from chr1 (10-15Mbp)
cooler show my_cooler.cool chr1:10000000-15000000

# Save heatmap as .png
cooler show -o my_matrix.png my_cooler.cool chr1:10000000-15000000

# Save high resolution heatmap as .png
cooler show -o my_matrix.png --dpi 600 my_cooler.cool chr1:10000000-15000000
```

For a better visualization experience we recommend using [HiGlass](https://github.com/higlass/higlass), in particular the containerized version of HiGlass which is installed and managed through [higlass-manage](https://github.com/higlass/higlass-manage).
