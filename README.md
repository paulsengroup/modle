<!--
Copyright (C) 2022 Roberto Rossini <roberros@uio.no>

SPDX-License-Identifier: MIT
-->

# MoDLE

[![Unit tests Ubuntu](https://github.com/paulsengroup/modle/actions/workflows/unit-tests-ubuntu.yml/badge.svg?branch=main)](https://github.com/paulsengroup/modle/actions/workflows/unit-tests-ubuntu.yml)
[![Build Docker image](https://github.com/paulsengroup/modle/actions/workflows/build-docker-image.yml/badge.svg)](https://github.com/paulsengroup/modle/actions/workflows/build-docker-image.yml)

## Using MoDLE

The recommended way to run MoDLE is using the Docker images hosted
on [ghcr.io](https://github.com/paulsengroup/modle/pkgs/container/modle)
and [dockerhub](https://hub.docker.com/repository/docker/paulsengroup/modle).

```bash
# Using Docker
sudo docker run ghcr.io/paulsengroup/modle:1.0.0-rc.3 --help

# Using Singularity/Apptainer
singularity run docker://ghcr.io/paulsengroup/modle:1.0.0-rc.3 --help
```

## Building MoDLE

### Pre-requisites

Building MoDLE requires a compiler toolchain supporting C++17, such as:

- GCC 7 and newer
- Clang 7 and newer
- Apple-Clang 10.0 and newer

MoDLE is developed on a Linux machine and should run on most UNIX-like OSes (including macOS).

It is in theory possible to compile and run MoDLE on Windows, but this is not officially supported.
If you really need to run MoDLE on Windows, please consider using the Docker images hosted
on [ghcr.io](https://github.com/paulsengroup/modle/pkgs/container/modle)
and [dockerhub](https://hub.docker.com/repository/docker/paulsengroup/modle).

In addition to a C++17 compiler, building MoDLE requires the following tools:

- CMake >= 3.16
- Conan >= 1.43

#### Installing Conan

Conan is a package manager for C and C++ applications, and can be installed using pip or Homebrew:

- `pip3 install "conan>=1.43"`
- `brew install conan`

### Getting MoDLE source code

We highly recommend users to download MoDLE's source code for the latest stable release from
the [Release](https://github.com/paulsengroup/modle/releases) page.

Using a simple `git clone` is only recommended if you intend to test features/bugfixes that have not yet landed in a
release.

### Compiling MoDLE

Run the following command from inside MoDLE's source tree.
Here we assume the machine where MoDLE will be compiled has 8 CPU cores.
Feel free to adjust `-j 8` to match the number of CPU cores available on your machine to improve compilation speed.

```bash
mkdir build/
cd build/

# Configure project
cmake ..

# Compile project
cmake --build . -j 8
```

#### Notes

By default, running the commands listed in
section [Installing MoDLE](https://github.com/paulsengroup/modle#installing-modle) will install MoDLE
under `/usr/local/` (i.e. the actual binary will be located at `/usr/local/bin/modle`).

Replace `cmake ..` with `cmake -DCMAKE_INSTALL_PREFIX="$HOME/.local/" ..` to install MoDLE for your user only.

The path passed to CMake through `-DCMAKE_INSTALL_PREFIX` can be in principle any path where your user has write
permissions.

### Testing MoDLE

Some of MoDLE's unit test depend on [SciPy](https://scipy.org/)
and [wCorr](https://cran.r-project.org/web/packages/wCorr/index.html), so make sure to have both packages installed
before running `ctest`.

Alternatively, pass `-E '(Corr\.)|(Binom) test` to `ctest` to exclude unit tests with external dependencies from the
test suite.

To launch the test suite, run the following commands from the repository root:

```bash
cd build/

ctest -j 8                 \
      --test-dir .         \
      --schedule-random    \
      --output-on-failure  \
      --no-tests=error
```

To launch the integration tests, run the following command from the repository root:

```bash
test/scripts/modle_integration_test_simple.sh build/src/modle/modle
```

### Installing MoDLE

The following command will install MoDLE files under the prefix specified through `-DCMAKE_INSTALL_PREFIX` (`/usr/local`
by default).

```bash
# Run from the repository root
cmake --install build/
```

### Troubleshooting build errors

#### Incorrect or incomplete Conan profile

In this case the build process will fail during project configuration (i.e. when running the first CMake command
in [this](https://github.com/paulsengroup/modle#compiling-modle) section) with an error message similar to the
following:

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

If after running the previous command you see a warning about `GCC OLD ABI COMPATIBILITY` run:

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

### Running MoDLE

Commands in this section assume you are running MoDLE using Singularity/Apptainer from the root of this repository.

Test datasets are located under `test/data/integration_tests`.

If you are running MoDLE without using containers, replace `singularity run docker://ghcr.io/paulsengroup/modle:1.0.0-rc.3` with `modle` in the coming examples.

#### Required input files

Running a simulation with default settings only requires two input files:
- A chrom.sizes with the list of chromosome to be simulated
- A BED file with the list of extrusion barriers to use in the simulation

The extrusion barrier BED file should have at least the first 6 columns defined (i.e. chrom, chromStart, chromEnd, name, score and strand).

The ___name___ field is ignored, while the ___score___ field is optional.

When ___score___ is non-zero, its value will be used to set the occupancy for the extrusion barrier defined by the current line.

The ___strand___ field is required and is used to define the extrusion barrier direction.
As of `v1.0.0-rc.3`, this field should be set to the direction of the CTCF motif.
Barriers without strand information (i.e. with strand '.') will be ignored.

Sample chrom.sizes and BED file(s) are available inside folder `test/data/integration_tests`.

#### Running a simulation with default settings
```bash
singularity run docker://ghcr.io/paulsengroup/modle:1.0.0-rc.3   \
    simulate \
    --chrom-sizes test/data/integration_tests/grch38.chrom.sizes \
    --extrusion-barrier-file test/data/integration_tests/grch38_h1_extrusion_barriers.bed.xz \
    --output-prefix path/to/ouput/prefix
```

This will create folder `path/to/output` (if it doesn't already exists), and write the following files inside it:
```
path/to/ouput
├── prefix_config.toml
├── prefix.cool
└── prefix.log
```

In case any of the output files already exist, MoDLE will refuse to run and print an error message listing the file name collisions.

Passing the `--force` flag overrides this behavior and will cause MoDLE to override existing files.

#### Running a simulation using config files

File `prefix_config.toml` from the previous section is a config file that can be used to re-run a simulation using the same parameters.

```bash
# Run a simulation with the same parameter as the previous example
singularity run docker://ghcr.io/paulsengroup/modle:1.0.0-rc.3 \
    --config path/to/output/prefix_config.toml
```

Parameters read from a config file have lower priority than parameter specified through the CLI,
meaning that CLI options can be used to override parameters read from the config file.

This can be useful when running a batch of simulations using a large number of optional parameters, and where only a handful of parameters need to change across simulation runs.

```bash
# The following command will run a simulation using parameters from the previous example as starting point,
# but using a custom lef density and overriding the output prefix specified by the config file.
singularity run docker://ghcr.io/paulsengroup/modle:1.0.0-rc.3 \
    simulate \
    --config path/to/output/prefix_config.toml \
    --lef-density 15 \
    --output-prefix path/to/ouput/prefix2
```

The config file is a text file in TOML format.

Adding a line like `my-option=my_value` to the config file is equivalent to passing `--my-option=my_value` on the CLI.

For an up-to-date list of supported CLI options, please refer to MoDLE's help message:
```bash
singularity run docker://ghcr.io/paulsengroup/modle:1.0.0-rc.3 simulate --help
```

#### Closing remarks
MoDLE automatically detects and handles compressed input files.

As of `v1.0.0-rc.3`, the following compression algorithms are supported:

- bzip2
- gzip
- LZ4
- LZO
- XZ/LZMA
<!-- - ZSTD -->
