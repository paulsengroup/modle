<!--
Copyright (C) 2022 Roberto Rossini <roberros@uio.no>

SPDX-License-Identifier: MIT
-->

# MoDLE

[![License](https://img.shields.io/badge/license-MIT-green)](./LICENSE)
[![Ubuntu CI](https://github.com/paulsengroup/modle/actions/workflows/ubuntu-ci.yml/badge.svg)](https://github.com/paulsengroup/modle/actions/workflows/ubuntu-ci.yml)
[![MacOS CI](https://github.com/paulsengroup/modle/actions/workflows/macos-ci.yml/badge.svg)](https://github.com/paulsengroup/modle/actions/workflows/macos-ci.yml)
[![Build Docker image](https://github.com/paulsengroup/modle/actions/workflows/build-docker-image.yml/badge.svg)](https://github.com/paulsengroup/modle/actions/workflows/build-docker-image.yml)
[![Download from Bioconda](https://img.shields.io/conda/vn/bioconda/modle?label=bioconda&logo=Anaconda)](https://anaconda.org/bioconda/modle)

[![MoDLE paper - Genome Biology 2022](https://img.shields.io/badge/CITE-Genome%20Biology%20(2022)-blue
)](https://doi.org/10.1186/s13059-022-02815-7)
[![Zenodo DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6424697.svg)](https://doi.org/10.5281/zenodo.6424697)
---

MoDLE is a computational tool for fast, stochastic modeling of molecular contacts from DNA loop extrusion capable of simulating realistic contact patterns genome wide in a few minutes.

![Figure 1 from MoDLE's paper](https://media.springernature.com/full/springer-static/image/art%3A10.1186%2Fs13059-022-02815-7/MediaObjects/13059_2022_2815_Fig1_HTML.png)

## Using MoDLE

MoDLE is developed on Linux and tested on Linux and MacOS.
Running MoDLE on Windows is only possible using Docker containers (see below).

MoDLE can be installed or compiled with one of the following methods.

### Docker or Singularity/Apptainer

First, make sure you follow the instructions on how to install Docker or Singularity/Apptainer on your OS.

<details>
<summary>Installing Docker</summary>

The following instructions assume you have root/admin permissions.

- [Linux](https://docs.docker.com/desktop/install/linux-install/#generic-installation-steps/)
- [MacOS](https://docs.docker.com/desktop/install/mac-install/)
- [Windows](https://docs.docker.com/desktop/install/windows-install/)

On some Linux distributions just installing Docker is not enough.
You also need to start (and optionally enable) the appropriate service(s).
This is usually done with one of the following:

```bash
sudo systemctl start docker
sudo systemctl start docker.service
```

Refer to [Docker](https://docs.docker.com/engine/install/) or your distribution documentation for more details.

</details>

<details>
<summary>Installing Singularity/Apptainer</summary>

The following instructions assume you have root/admin permissions.

Apptainer can be easily installed using your system package manager.

[Here](https://apptainer.org/docs/admin/main/installation.html#install-from-pre-built-packages) you can find instructions for common Linux distributions such as Ubuntu.

Even if your distribution is not listed in the above documentation, your system package manager likely includes a package for Singularity or Apptainer. If this is not the case, then you must install Apptainer from source (instructions available [here](https://github.com/apptainer/apptainer/blob/release-1.1/INSTALL.md)).

</details>

#### Pulling MoDLE Docker image

MoDLE Docker images are available on [ghcr.io](https://github.com/paulsengroup/modle/pkgs/container/modle)
and [dockerhub](https://hub.docker.com/repository/docker/paulsengroup/modle).

Downloading and running the latest stable release can be done as follows:

```console
# Using Docker, may require sudo
user@dev:/tmp$ docker run ghcr.io/paulsengroup/modle:1.0.0 --help

# Using Singularity/Apptainer
user@dev:/tmp$ singularity run ghcr.io/paulsengroup/modle:1.0.0 --help

High-performance stochastic modeling of DNA loop extrusion interactions.
Usage: /usr/local/bin/modle [OPTIONS] SUBCOMMAND

Options:
-h,--help                   Print this help message and exit
-V,--version                Display program version information and exit
--config                    Path to MoDLE's config file (optional).

Subcommands:
simulate, sim               Simulate loop extrusion and write resulting molecular contacts in a .cool file.
```

The above will print MoDLE's help message, and is equivalent to running `modle --help` on the command line (assuming MoDLE is available on your machine).

<details>
<summary>Troubleshooting</summary>

<!-- TODO: Should this section be moved elsewhere? -->

__Errors regarding missing input file(s)__

Assuming you have placed the required input files inside a folder named `mydata`

```console
user@dev:/tmp$ ls -lah mydata
total 248K
drwx------  2 user group   80 Jan  5 12:00 .
drwxrwxrwt 20 user group  540 Jan  5 12.00 ..
-rw-------  1 user group  365 Jan  5 12:00 grch38.chrom.sizes
-rw-------  1 user group 241K Jan  5 12:00 grch38_h1_extrusion_barriers.bed.xz
```

Running MoDLE using Docker as shown below fails due to missing input file(s)
```bash
docker run ghcr.io/paulsengroup/modle:1.0.0 simulate \
  --chrom-sizes mydata/grch38.chrom.sizes \
  --extrusion-barrier-file mydata/grch38_h1_extrusion_barriers.bed.xz \
  --output-prefix mydata/output

# Error message
# --chrom-sizes: File does not exist: mydata/grch38.chrom.sizes
# Run with --help for more information.
```

This issue is caused by Docker not mounting/binding the host filesystem (e.g. the filesystem where `mydata/` resides) inside the container.

To make `mydata/` accessible from inside the container, use one of the following methods

```bash
# Method 1: mydata/ folder on the host is mapped to /mydata (notice the /) inside the container
docker run -v "$(pwd -P)/mydata/:/data/" \
  ghcr.io/paulsengroup/modle:1.0.0 simulate \
  --chrom-sizes /data/grch38.chrom.sizes \
  --extrusion-barrier-file /data/grch38_h1_extrusion_barriers.bed.xz \
  --output-prefix /data/output

# Mehtod 2: use two different folders for input/output. Input folder is mounted in read-only mode
docker run -v "$(pwd -P)/mydata/:/input_data/:ro" \
           -v "$(pwd -P)/output/:/output_data/" \
  ghcr.io/paulsengroup/modle:1.0.0 simulate \
  --chrom-sizes /input_data/grch38.chrom.sizes \
  --extrusion-barrier-file /input_data/grch38_h1_extrusion_barriers.bed.xz \
  --output-prefix /output_data/output

# Methods for Singularity/Apptainer (note that manually mounting paths is usually not required when using Singularity/Apptainer)
singularity run -B "$(pwd -P)/mydata/:/data/" \
  docker://ghcr.io/paulsengroup/modle:1.0.0 simulate \
  --chrom-sizes /data/grch38.chrom.sizes \
  --extrusion-barrier-file /data/grch38_h1_extrusion_barriers.bed.xz \
  --output-prefix /data/output

singularity run -B "$(pwd -P)/mydata/:/input_data/:ro" \
                -B "$(pwd -P)/output/:/output_data/" \
  docker://ghcr.io/paulsengroup/modle:1.0.0 simulate \
  --chrom-sizes /input_data/grch38.chrom.sizes \
  --extrusion-barrier-file /input_data/grch38_h1_extrusion_barriers.bed.xz \
  --output-prefix /output_data/output
```

Output produced by Method 1:

```console
user@dev:/tmp$ ls -lah mydata
total 57M
drwx------  2 user group  140 Jan  5 12:20 .
drwxrwxrwt 20 user group  560 Jan  5 12:20 ..
-rw-------  1 user group  365 Jan  5 12:00 grch38.chrom.sizes
-rw-------  1 user group 241K Jan  5 12:00 grch38_h1_extrusion_barriers.bed.xz
-rw-r--r--  1 root root  3.7K Jan  5 12:00 output_config.toml
-rw-r--r--  1 root root   57M Jan  5 12:05 output.cool
-rw-r--r--  1 root root   14K Jan  5 12:05 output.log
```

Output produced by Method 2:
```console
user@dev:/tmp$ ls -lah mydata/
total 57M
drwx------  2 user group  140 Jan  5 12:20 .
drwxrwxrwt 20 user group  560 Jan  5 12:20 ..
-rw-------  1 user group  365 Jan  5 12:00 grch38.chrom.sizes
-rw-------  1 user group 241K Jan  5 12:00 grch38_h1_extrusion_barriers.bed.xz

user@dev:/tmp$ ls -lah output/
drwx------  2 user group  140 Jan  5 12:20 .
drwxrwxrwt 20 user group  560 Jan  5 12:20 ..
-rw-r--r--  1 root root  3.7K Jan  5 12:00 output_config.toml
-rw-r--r--  1 root root   57M Jan  5 12:05 output.cool
-rw-r--r--  1 root root   14K Jan  5 12:05 output.log
```

__Troubleshooting problems inside the container is a bit tedious. Any tips?__

Yes! For troubleshooting we recommend using an interactive shell running inside the container.

```bash
docker run -v "$(pwd -P)/mydata/:/input_data/:ro" \
           -v "$(pwd -P)/output/:/output_data/" \
           -it --rm --entrypoint=/bin/bash \
  ghcr.io/paulsengroup/modle:1.0.0
```

This will drop you in an interactive bash shell.
Press `Ctrl+C`, `Ctrl+D` or type `exit` to leave the shell and stop the container.

```console
root@25e8efe74294:/data# ls -lah
total 4.0K
drwxr-xr-x  2 root root  2 Jan  5 12:00 .
drwxr-xr-x 20 root root 27 Jan  5 12:00 ..
root@25e8efe74294:/data# cd /input_data/
root@25e8efe74294:/input_data# ls -lah
total 250K
drwx------  2 1000 1000   80 Jan  5 12:00 .
drwxr-xr-x 20 root root   27 Jan  5 12:00 ..
-rw-------  1 1000 1000  365 Jan  5 12:00 grch38.chrom.sizes
-rw-------  1 1000 1000 241K Jan  5 12:00 grch38_h1_extrusion_barriers.bed.xz
root@25e8efe74294:/input_data# touch foo
touch: cannot touch 'foo': Read-only file system
root@25e8efe74294:/input_data# cd /output_data/
root@25e8efe74294:/output_data# ls -lah
total 2.0K
drwxr-xr-x  2 1000 1000 40 Jan  5 12:00 .
drwxr-xr-x 20 root root 27 Jan  5 12:00 ..
root@25e8efe74294:/output_data# touch foo
root@25e8efe74294:/output_data# ls -lah
total 2.0K
drwxr-xr-x  2 1000 1000 60 Jan  5 12:00 .
drwxr-xr-x 20 root root 27 Jan  5 12:00 ..
-rw-r--r--  1 root root  0 Jan  5 12:01 foo
root@25e8efe74294:/output_data# whereis modle
modle: /usr/local/bin/modle
root@25e8efe74294:/output_data# modle --help
High-performance stochastic modeling of DNA loop extrusion interactions.
Usage: modle [OPTIONS] SUBCOMMAND

Options:
  -h,--help                   Print this help message and exit
  -V,--version                Display program version information and exit
  --config                    Path to MoDLE's config file (optional).

Subcommands:
  simulate, sim               Simulate loop extrusion and write resulting molecular contacts in a .cool file.

root@25e8efe74294:/output_data# exit
exit
```

For Singularity/Apptainer this is even easier:

```bash
singularity shell docker://ghcr.io/paulsengroup/modle:1.0.0
```

</details>


### Conda (bioconda)

MoDLE package for Linux and MacOS is available on [bioconda](https://anaconda.org/bioconda/modle) and can be installed as follows:

```console
user@dev:/tmp$ conda create -n modle -c conda-forge -c bioconda modle

(modle) user@dev:/tmp$ conda activate modle

(modle) user@dev:/tmp$ whereis modle
modle: /home/user/.miniconda3/envs/modle/bin/modle

(modle) user@dev:/tmp$ modle --version
MoDLE-v1.0.0
```
<!-- TODO: add conda-specific version suffix -->

### Pre-built Linux binaries (not recommended)

Download file `modle-1.0.0-x86_64-linux.tar.xz` from the [latest release](https://github.com/paulsengroup/modle/releases/latest).

Extract the archive and copy `modle` and/or `modle_tools` to a location in your `PATH`

<details>
<summary> click to expand </summary>

<!-- TODO make sure the link is correct once the release goes live -->
```console
user@dev:/tmp$ curl -LO 'https://github.com/paulsengroup/modle/releases/download/v1.0.0/modle-1.0.0-x86_64-linux.tar.xz'

user@dev:/tmp$ tar -xvf modle-1.0.0-x86_64-linux.tar.xz
modle-1.0.0-x86_64-linux/
modle-1.0.0-x86_64-linux/share/
modle-1.0.0-x86_64-linux/share/licenses/
modle-1.0.0-x86_64-linux/share/licenses/modle_tools/
modle-1.0.0-x86_64-linux/share/licenses/modle_tools/LICENSE
modle-1.0.0-x86_64-linux/share/licenses/modle/
modle-1.0.0-x86_64-linux/share/licenses/modle/LICENSE
modle-1.0.0-x86_64-linux/bin/
modle-1.0.0-x86_64-linux/bin/modle_tools
modle-1.0.0-x86_64-linux/bin/modle

user@dev:/tmp modle-1.0.0-x86_64-linux/bin/modle --version
MoDLE-v1.0.0

# Optional: add modle and modle_tools to your PATH
user@dev:/tmp$ install -Dm0755 modle-1.0.0-x86_64-linux/bin/modle "$HOME/.local/bin/"
user@dev:/tmp$ install -Dm0755 modle-1.0.0-x86_64-linux/bin/modle_tools "$HOME/.local/bin/"

user@dev:/tmp$ whereis modle
modle: /home/user/.local/bin/modle

user@dev:/tmp$ modle --version
MoDLE-v1.0.0-bioconda
```

</details>

### Installing from source

MoDLE can be compiled on most UNIX-like systems, including many Linux distributions and MacOS (10.15+).

<details>
<summary> <b>Build instructions</b> </summary>

##### Build requirements

Compiling MoDLE requires a compiler toolchain supporting C++17, such as:

- GCC 8+
- Clang 8+
- Apple-Clang 10.0+

Furthermore, the following tools are required:
- CMake 3.20+
- Conan 2+
- git 2.7+
- make or ninja
- Python3.6+ (including `pip`, required to install Conan)


We recommend to install CMake and Conan in a Python [virtualenv](https://virtualenvwrapper.readthedocs.io/en/stable/), but you are of course free to install the build dependencies in any way you want.

```bash
python3 -m venv /tmp/venv
/tmp/venv/bin/python3 -m pip install pip setuptools --upgrade
/tmp/venv/bin/python3 -m pip  install 'cmake>=3.20' 'conan>=2' ninja

# NOTE: It's important to activate the venv after installing CMake
. /tmp/venv/bin/activate

whereis cmake  # cmake: /tmp/venv/bin/cmake
whereis conan  # conan: /tmp/venv/bin/conan
whereis ninja  # ninja: /tmp/venv/bin/ninja

cmake --version
conan --version

# Detect compiler toolchain. It is usually a good idea to explicitly set CC and CXX
CC=gcc CXX=g++ conan profile detect --force
```

Compiling MoDLE on Windows using MSVC is currently not possible.
Compiling with MinGW may be possible but has never been tested.
If you really want to compile MoDLE on Windows, then you can try to do so using the [Windows Subsystem for Linux (WSL)](https://learn.microsoft.com/en-us/windows/wsl/install).

##### Getting the source code

Download from the [Release](https://github.com/paulsengroup/modle/releases) page (recommended).
```bash
mkdir /tmp/modle
curl -L 'https://github.com/paulsengroup/modle/archive/refs/tags/v1.0.0.tar.gz' | tar --strip-components=1 -C /tmp/modle -xzf -
```

Using git.
```bash
git clone https://github.com/paulsengroup/modle.git /tmp/modle

cd /tmp/modle
git checkout v1.0.0  # Skip this step if you want to build the latest commit from main
```

##### Compiling MoDLE

```bash
# Activate venv
. /tmp/venv/bin/activate

# Set these variables to the number of CPU cores available on your machine
# You can check this with e.g.
# python -c 'import multiprocessing as mp; print(mp.cpu_count())')
export CONAN_CPU_COUNT=8
export CMAKE_BUILD_PARALLEL_LEVEL=8

# Install/build dependencies with Conan
conan install --build=missing \
              --build=cascade \
              --update \
              -pr default \
              -s build_type=Release \
              -s compiler.cppstd=17 \
              --output-folder=./build/ \
              .

# This may take a while, as CMake will run Conan to build MoDLE dependencies.
# Do not pass -G Ninja if you want CMake to use make instead of ninja
cmake -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_PREFIX_PATH="$PWD/build" \
      -DMODLE_ENABLE_TESTING=ON  \
      -G Ninja \
      -S /tmp/modle \
      -B /tmp/modle/build

cmake --build /tmp/modle/build
```

To override the default compiler used by CMake, pass the following arguments to the first CMake command: `-DCMAKE_C_COMPILER=path/to/cc -DCMAKE_CXX_COMPILER=path/to/c++`

We highly recommend using the same compiler when running Conan and CMake.

</details>

<details>
<summary> <b> Running automated tests </b> </summary>

Steps outlined in this section are optional but highly recommended.

##### Unit tests

```bash
# Activate venv
. /tmp/venv/bin/activate

cd /tmp/modle
ctest --test-dir build/   \
      --schedule-random   \
      --output-on-failure \
      --no-tests=error    \
      --timeout 600       \
      -j8  # Change this to the number of available CPU cores
```

A successful run of the test suite will produce an output like the following:
```console
user@dev:/tmp/modle$ ctest --test-dir build/ ...
...
124/127 Test #120: ContactMatrix internal: encode/decode roundtrip - LONG ...........................................   Passed   35.05 sec
125/127 Test #117: ContactMatrixSparse serialization - LONG .........................................................   Passed   42.19 sec
126/127 Test #119: ContactMatrixSparse serialization (FP) - LONG ....................................................   Passed   44.71 sec
127/127 Test #118: ContactMatrixSparse serialization (destructive) - LONG ...........................................   Passed   32.53 sec

100% tests passed, 0 tests failed out of 127

Total Test time (real) =  45.83 sec
```

__All tests are expected to pass. Do not ignore test failures!__

<details>
<summary> Troubleshooting test failures </summary>
If one or more test fail, try the following troubleshooting steps before reaching out for help.

1. Make sure you are running `ctest` from the root of the source tree (`/tmp/modle` if you are following the instructions).
2. Make sure you are passing the correct build folder to `--test-dir`. Pass the absolute path if necessary (i.e. `--test-dir=/tmp/modle/build/` if you are following the instructions).
3. Re-run `ctest` with `-j1`. This can be necessary on machines with very little memory (e.g. less than 2GB).
4. Before running `ctest`, create a temporary folder where your user has read-write permissions and where there are at least 100-200MB of space available.
   Then set variable `TMPDIR` to that folder and re-run `ctest`.
5. Checksum the test dataset located under `test/data/` by running `sha256sum -c checksums.sha256`.
   If the checksumming fails or the folder doesn't exist, download and extract the `.tar.xz` file listed in file `cmake/FetchTestDataset.cmake`. Make sure you run `tar -xf` from the root of the repository (`/tmp/modle` if you are following the instructions).

Example:
```bash
# Activate venv
. /tmp/venv/bin/activate

cd /tmp/modle

# Make sure this is the URL listed  in file cmake/FetchTestDataset.cmake
curl -L 'https://zenodo.org/record/7506960/files/modle_test_data.tar.xz?download=1' | tar -xJf -

# This should print "OK" if the check is successful
(cd test/data && sha256sum --quiet -c checksums.sha256 && 2>&1 echo OK)

mkdir ~/modle-test-dir  # Remember to delete this folder

TMPDIR="$HOME/modle-test-dir"      \
ctest --test-dir=/tmp/modle/build/ \
      --schedule-random            \
      --output-on-failure          \
      --no-tests=error             \
      --timeout 600                \
      -j1

# rm -r ~/modle-test-dir
```

If after trying the above steps the tests are still failing, feel free to start [discussion](https://github.com/paulsengroup/modle/discussions) asking for help.

</details>


##### Integration tests (Linux only)

The integration test scripts depend on the following tools:

- cooler>=0.8.11
- xz
- common UNIX shell commands such as (namely `cmp`, `diff`, `grep` and `mktemp`)

cooler can be installed using pip:
```bash
# See this PR for why we need numpy<1.24
# https://github.com/open2c/cooler/pull/298

/tmp/venv/bin/pip3 install cython 'numpy<1.24'
/tmp/venv/bin/pip3 install 'cooler>=0.8.11'
```

If not already installed, `xz` can usually be installed with your system package manager (on some Linux distributions the relevant package is called `xz-utils`).

```bash
# Activate venv
. /tmp/venv/bin/activate

# Test modle
test/scripts/modle_integration_test.sh build/src/modle/modle

# Test modle_tools
test/scripts/modle_tools_transform_integration_test.sh build/src/modle_tools/modle_tools
test/scripts/modle_tools_eval_integration_test.sh build/src/modle_tools/modle_tools
```

<details>
<summary> Sample outputs of successful test runs </summary>

```console
user@dev:/tmp$ test/scripts/modle_integration_test.sh build/src/modle/modle
[2023-01-14 22:13:25.505] [info]: Running MoDLE v1.0.0
[2023-01-14 22:13:25.506] [info]: Complete log will be written to file "/tmp/modle-QCJPRtcPHp/out.log"
[2023-01-14 22:13:25.506] [info]: Writing simulation parameters to config file "/tmp/modle-QCJPRtcPHp/out_config.toml"
[2023-01-14 22:13:25.506] [info]: Command: build/src/modle/modle sim -c /dev/fd/63 -b /dev/fd/62 -o /tmp/modle-QCJPRtcPHp/out -r 20kb --target-contact-density 20 --ncells 2 --track-1d-lef-position --max-burnin-epochs 5000
[2023-01-14 22:13:25.506] [info]: Simulation will use up to 16 out of 16 available CPU cores.
[2023-01-14 22:13:25.506] [info]: Using --target-contact-density=20.00 as stopping criterion.
[2023-01-14 22:13:25.506] [info]: Contact sampling strategy: tad-plus-loop-with-noise.
[2023-01-14 22:13:25.506] [info]: Contact matrix resolution: 20000bp
[2023-01-14 22:13:25.506] [info]: Importing chromosomes from file "/dev/fd/63"...
[2023-01-14 22:13:25.506] [info]: Imported 2 chromosomes in 420.684us.
[2023-01-14 22:13:25.506] [info]: Importing extrusion barriers from file "/dev/fd/62"...
[2023-01-14 22:13:25.523] [info]: Imported 3160 barriers in 17.105856ms.
[2023-01-14 22:13:25.526] [info]: Spawning simulation thread 0...
[2023-01-14 22:13:25.526] [info]: Spawning simulation thread 2...
[2023-01-14 22:13:25.526] [info]: Spawning simulation thread 5...
[2023-01-14 22:13:25.526] [info]: Spawning simulation thread 3...
[2023-01-14 22:13:25.526] [info]: Spawning simulation thread 11...
[2023-01-14 22:13:25.526] [info]: Spawning simulation thread 13...
[2023-01-14 22:13:25.526] [info]: Spawning simulation thread 1...
[2023-01-14 22:13:25.526] [info]: Spawning simulation thread 6...
[2023-01-14 22:13:25.526] [info]: Spawning simulation thread 15...
[2023-01-14 22:13:25.526] [info]: Spawning simulation thread 9...
[2023-01-14 22:13:25.526] [info]: Spawning simulation thread 10...
[2023-01-14 22:13:25.526] [info]: Spawning simulation thread 12...
[2023-01-14 22:13:25.526] [info]: Spawning simulation thread 4...
[2023-01-14 22:13:25.526] [info]: Spawning simulation thread 14...
[2023-01-14 22:13:25.526] [info]: Spawning simulation thread 7...
[2023-01-14 22:13:25.526] [info]: Spawning simulation thread 8...
[2023-01-14 22:13:25.527] [info]: Begin processing "chr19": simulating ~11724 epochs across 2 cells using 1172 LEFs and 1441 barriers (~5862 epochs per cell)...
[2023-01-14 22:13:25.527] [info]: Begin processing "chr17": simulating ~11716 epochs across 2 cells using 1665 LEFs and 1719 barriers (~5858 epochs per cell)...
[2023-01-14 22:13:27.870] [info]: Simulation of "chr19" successfully completed.
[2023-01-14 22:13:28.791] [info]: Simulation of "chr17" successfully completed.
[2023-01-14 22:13:28.851] [info]: Writing contacts for "chr17" to file "/tmp/modle-QCJPRtcPHp/out.cool"...
[2023-01-14 22:13:28.899] [info]: Written 12489000 contacts for "chr17" to file "/tmp/modle-QCJPRtcPHp/out.cool" (0.15M nnz out of 0.62M pixels).
[2023-01-14 22:13:28.899] [info]: Writing contacts for "chr19" to file "/tmp/modle-QCJPRtcPHp/out.cool"...
[2023-01-14 22:13:28.922] [info]: Written 8793000 contacts for "chr19" to file "/tmp/modle-QCJPRtcPHp/out.cool" (0.10M nnz out of 0.44M pixels).
[2023-01-14 22:13:29.094] [info]: Simulation terminated without errors in 3.587667161s!

Bye.
Comparing /tmp/modle-QCJPRtcPHp/out.cool with test/data/integration_tests/reference_001.cool...
Files are identical
Comparing /tmp/modle-QCJPRtcPHp/out_lef_1d_occupancy.bw with test/data/integration_tests/reference_001.bw...
Files are identical

### PASS ###


user@dev:/tmp$ test/scripts/modle_tools_eval_integration_test.sh buils/src/modle_tools/modle_tools
[2023-01-14 22:18:03.209] [info]: Computing metric(s) for the following 24 intervals:
 - chr1:0-248956422
 - chr2:0-242193529
 - chr3:0-198295559
 - chr4:0-190214555
 - chr5:0-181538259
 - chr6:0-170805979
 - chr7:0-159345973
 - chr8:0-145138636
 - chr9:0-138394717
 - chr10:0-133797422
 - chr11:0-135086622
 - chr12:0-133275309
 - chr13:0-114364328
 - chr14:0-107043718
 - chr15:0-101991189
 - chr16:0-90338345
 - chr17:0-83257441
 - chr18:0-80373285
 - chr19:0-58617616
 - chr20:0-64444167
 - chr21:0-46709983
 - chr22:0-50818468
 - chrX:0-156040895
 - chrY:0-57227415
[2023-01-14 22:18:03.214] [warning]: Read 0 contacts for chr1:248956422. SKIPPING!
[2023-01-14 22:18:03.215] [warning]: Read 0 contacts for chr2:242193529. SKIPPING!
[2023-01-14 22:18:03.216] [warning]: Read 0 contacts for chr3:198295559. SKIPPING!
[2023-01-14 22:18:03.217] [warning]: Read 0 contacts for chr4:190214555. SKIPPING!
[2023-01-14 22:18:03.218] [warning]: Read 0 contacts for chr5:181538259. SKIPPING!
[2023-01-14 22:18:03.218] [warning]: Read 0 contacts for chr6:170805979. SKIPPING!
[2023-01-14 22:18:03.219] [warning]: Read 0 contacts for chr7:159345973. SKIPPING!
[2023-01-14 22:18:03.220] [warning]: Read 0 contacts for chr8:145138636. SKIPPING!
[2023-01-14 22:18:03.221] [warning]: Read 0 contacts for chr9:138394717. SKIPPING!
[2023-01-14 22:18:03.222] [warning]: Read 0 contacts for chr10:133797422. SKIPPING!
[2023-01-14 22:18:03.222] [warning]: Read 0 contacts for chr11:135086622. SKIPPING!
[2023-01-14 22:18:03.223] [warning]: Read 0 contacts for chr12:133275309. SKIPPING!
[2023-01-14 22:18:03.224] [warning]: Read 0 contacts for chr13:114364328. SKIPPING!
[2023-01-14 22:18:03.225] [warning]: Read 0 contacts for chr14:107043718. SKIPPING!
[2023-01-14 22:18:03.225] [warning]: Read 0 contacts for chr15:101991189. SKIPPING!
[2023-01-14 22:18:03.226] [warning]: Read 0 contacts for chr16:90338345. SKIPPING!
[2023-01-14 22:18:03.227] [warning]: Read 0 contacts for chr17:83257441. SKIPPING!
[2023-01-14 22:18:03.228] [warning]: Read 0 contacts for chr18:80373285. SKIPPING!
[2023-01-14 22:18:03.229] [warning]: Read 0 contacts for chr19:58617616. SKIPPING!
[2023-01-14 22:18:03.229] [info]: Reading contacts for chr20:64444167...
[2023-01-14 22:18:03.235] [info]: Read 56772 contacts for chr20:64444167 in 5.904802ms
[2023-01-14 22:18:03.235] [info]: Custom metric for vertical stripes from interval chr20:0-64444167 computed in 419.907us.
[2023-01-14 22:18:03.236] [info]: Custom metric for horizontal stripes from interval chr20:0-64444167 computed in 1.522523ms.
[2023-01-14 22:18:03.256] [info]: 2578 values have been written to files "out_custom_score_custom_metric_vertical.{tsv.gz,bw}" in 20.776087ms.
[2023-01-14 22:18:03.258] [info]: 2578 values have been written to files "out_custom_score_custom_metric_horizontal.{tsv.gz,bw}" in 21.645827ms.
[2023-01-14 22:18:03.258] [warning]: Read 0 contacts for chr21:46709983. SKIPPING!
[2023-01-14 22:18:03.258] [warning]: Read 0 contacts for chr22:50818468. SKIPPING!
[2023-01-14 22:18:03.258] [warning]: Read 0 contacts for chrX:156040895. SKIPPING!
[2023-01-14 22:18:03.258] [warning]: Read 0 contacts for chrY:57227415. SKIPPING!
[2023-01-14 22:18:03.258] [info]: DONE in 77.414949ms!
Comparing test/data/integration_tests/4DNFI9GMP2J8_vs_4DNFIFJH2524_mt_eval_custom_score_custom_metric_horizontal.bw with /tmp/modle-tools-XxmbEEIUgY/out_custom_score_custom_metric_horizontal.bw...
Files are identical
Comparing test/data/integration_tests/4DNFI9GMP2J8_vs_4DNFIFJH2524_mt_eval_custom_score_custom_metric_vertical.bw with /tmp/modle-tools-XxmbEEIUgY/out_custom_score_custom_metric_vertical.bw...
Files are identical
Comparing test/data/integration_tests/4DNFI9GMP2J8_vs_4DNFIFJH2524_mt_eval_custom_score_custom_metric_horizontal.tsv.gz with /tmp/modle-tools-XxmbEEIUgY/out_custom_score_custom_metric_horizontal.tsv.gz...
Files are identical
Comparing test/data/integration_tests/4DNFI9GMP2J8_vs_4DNFIFJH2524_mt_eval_custom_score_custom_metric_vertical.tsv.gz with /tmp/modle-tools-XxmbEEIUgY/out_custom_score_custom_metric_vertical.tsv.gz...
Files are identical

### PASS ###


user@dev:/tmp$  test/scripts/modle_tools_transform_integration_test.sh build/src/modle_tools/modle_tools
[2023-01-14 22:20:46.906] [info]: Transforming contacts from Cooler at URI "test/data/integration_tests/4DNFI9GMP2J8_chr20_25kbp.cool"...
[2023-01-14 22:20:46.911] [warning]: Read 0 contacts for chr1. SKIPPING!
[2023-01-14 22:20:46.911] [warning]: Read 0 contacts for chr2. SKIPPING!
[2023-01-14 22:20:46.912] [warning]: Read 0 contacts for chr3. SKIPPING!
[2023-01-14 22:20:46.913] [warning]: Read 0 contacts for chr4. SKIPPING!
[2023-01-14 22:20:46.914] [warning]: Read 0 contacts for chr5. SKIPPING!
[2023-01-14 22:20:46.914] [warning]: Read 0 contacts for chr6. SKIPPING!
[2023-01-14 22:20:46.915] [warning]: Read 0 contacts for chr7. SKIPPING!
[2023-01-14 22:20:46.916] [warning]: Read 0 contacts for chr8. SKIPPING!
[2023-01-14 22:20:46.917] [warning]: Read 0 contacts for chr9. SKIPPING!
[2023-01-14 22:20:46.918] [warning]: Read 0 contacts for chr10. SKIPPING!
[2023-01-14 22:20:46.918] [warning]: Read 0 contacts for chr11. SKIPPING!
[2023-01-14 22:20:46.919] [warning]: Read 0 contacts for chr12. SKIPPING!
[2023-01-14 22:20:46.920] [warning]: Read 0 contacts for chr13. SKIPPING!
[2023-01-14 22:20:46.921] [warning]: Read 0 contacts for chr14. SKIPPING!
[2023-01-14 22:20:46.921] [warning]: Read 0 contacts for chr15. SKIPPING!
[2023-01-14 22:20:46.922] [warning]: Read 0 contacts for chr16. SKIPPING!
[2023-01-14 22:20:46.923] [warning]: Read 0 contacts for chr17. SKIPPING!
[2023-01-14 22:20:46.924] [warning]: Read 0 contacts for chr18. SKIPPING!
[2023-01-14 22:20:46.925] [warning]: Read 0 contacts for chr19. SKIPPING!
[2023-01-14 22:20:46.925] [info]: Processing contacts for chr20...
[2023-01-14 22:20:47.090] [info]: Applying Gaussian blur with sigma=1 to contacts for chr20...
[2023-01-14 22:20:47.104] [info]: chr20 processing took 178.637745ms
[2023-01-14 22:20:47.172] [warning]: Read 0 contacts for chr21. SKIPPING!
[2023-01-14 22:20:47.172] [warning]: Read 0 contacts for chr22. SKIPPING!
[2023-01-14 22:20:47.172] [warning]: Read 0 contacts for chrX. SKIPPING!
[2023-01-14 22:20:47.172] [warning]: Read 0 contacts for chrY. SKIPPING!
[2023-01-14 22:20:47.172] [info]: DONE! Processed 24 chromosomes in 266.321673ms!
[2023-01-14 22:20:47.172] [info]: Transformed contacts have been saved to file "/tmp/modle-tools-kA5f4J1HMA/out_blurred.cool"
[2023-01-14 22:20:47.380] [info]: Transforming contacts from Cooler at URI "test/data/integration_tests/4DNFI9GMP2J8_chr20_25kbp.cool"...
[2023-01-14 22:20:47.385] [warning]: Read 0 contacts for chr1. SKIPPING!
[2023-01-14 22:20:47.386] [warning]: Read 0 contacts for chr2. SKIPPING!
[2023-01-14 22:20:47.386] [warning]: Read 0 contacts for chr3. SKIPPING!
[2023-01-14 22:20:47.387] [warning]: Read 0 contacts for chr4. SKIPPING!
[2023-01-14 22:20:47.388] [warning]: Read 0 contacts for chr5. SKIPPING!
[2023-01-14 22:20:47.389] [warning]: Read 0 contacts for chr6. SKIPPING!
[2023-01-14 22:20:47.389] [warning]: Read 0 contacts for chr7. SKIPPING!
[2023-01-14 22:20:47.390] [warning]: Read 0 contacts for chr8. SKIPPING!
[2023-01-14 22:20:47.391] [warning]: Read 0 contacts for chr9. SKIPPING!
[2023-01-14 22:20:47.392] [warning]: Read 0 contacts for chr10. SKIPPING!
[2023-01-14 22:20:47.392] [warning]: Read 0 contacts for chr11. SKIPPING!
[2023-01-14 22:20:47.393] [warning]: Read 0 contacts for chr12. SKIPPING!
[2023-01-14 22:20:47.394] [warning]: Read 0 contacts for chr13. SKIPPING!
[2023-01-14 22:20:47.394] [warning]: Read 0 contacts for chr14. SKIPPING!
[2023-01-14 22:20:47.395] [warning]: Read 0 contacts for chr15. SKIPPING!
[2023-01-14 22:20:47.396] [warning]: Read 0 contacts for chr16. SKIPPING!
[2023-01-14 22:20:47.397] [warning]: Read 0 contacts for chr17. SKIPPING!
[2023-01-14 22:20:47.397] [warning]: Read 0 contacts for chr18. SKIPPING!
[2023-01-14 22:20:47.398] [warning]: Read 0 contacts for chr19. SKIPPING!
[2023-01-14 22:20:47.399] [info]: Processing contacts for chr20...
[2023-01-14 22:20:47.565] [info]: Computing the difference of Gaussians for chr20 (sigma1=1; sigma2=1.6)...
[2023-01-14 22:20:47.597] [info]: chr20 processing took 197.779235ms
[2023-01-14 22:20:47.665] [warning]: Read 0 contacts for chr21. SKIPPING!
[2023-01-14 22:20:47.665] [warning]: Read 0 contacts for chr22. SKIPPING!
[2023-01-14 22:20:47.665] [warning]: Read 0 contacts for chrX. SKIPPING!
[2023-01-14 22:20:47.665] [warning]: Read 0 contacts for chrY. SKIPPING!
[2023-01-14 22:20:47.665] [info]: DONE! Processed 24 chromosomes in 284.595243ms!
[2023-01-14 22:20:47.665] [info]: Transformed contacts have been saved to file "/tmp/modle-tools-kA5f4J1HMA/out_dog.cool"
Comparing /tmp/modle-tools-kA5f4J1HMA/out_blurred.cool with test/data/integration_tests/4DNFI9GMP2J8_chr20_25kbp_blurred.cool...
Files are identical
Comparing /tmp/modle-tools-kA5f4J1HMA/out_dog.cool with test/data/integration_tests/4DNFI9GMP2J8_chr20_25kbp_dog.cool...
Files are identical

### PASS ###
```

</details>

</details>

<details>
<summary> <b> Installation </b> </summary>

Once all tests have passed, `modle` and `modle_tools` can be installed as follows:

```console
# Activate venv
user@dev:/tmp$ . /tmp/venv/bin/activate

# Install system-wide (requires root/admin rights)
user@dev:/tmp$ cmake --install /tmp/modle/build
-- Install configuration: "Release"
-- Installing: /usr/local/bin/modle_tools
-- Installing: /usr/local/share/licenses/modle_tools/LICENSE
-- Installing: /usr/local/bin/modle
-- Installing: /usr/local/share/licenses/modle/LICENSE

# Alternatively, install to custom path
user@dev:/tmp$ cmake --install /tmp/modle/build --prefix "$HOME/.local/"
-- Install configuration: "Release"
-- Installing: /home/user/.local/bin/modle_tools
-- Installing: /home/user/.local/share/licenses/modle_tools/LICENSE
-- Installing: /home/user/.local/bin/modle
-- Installing: /home/user/.local/share/licenses/modle/LICENSE
```

</details>

<details>
<summary> <b> Cleaning build artifacts </b> </summary>

After successfully compiling MoDLE the following folders safely be removed:
- Python virtualenv: `/tmp/venv`
- MoDLE source tree: `/tmp/modle`

If you are not using Conan in any other project feel free to also delete Conan's folder `~/.conan2/`

</details>

### Running MoDLE

Here we assume MoDLE was correctly installed and is available in your PATH (i.e. running `whereis modle` prints something like `modle: /usr/local/bin/modle`).

Folder `examples/data/` contains the input files mentioned in this section.

<details>
<summary> Data provenance </summary>

- `examples/data/hg38.chrom.sizes` was generated from the [hg38.chrom.sizes](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes) released by UCSC.
Chromosmes were filtered and sorted using `grep -E 'chr[0-9XY]+[[:space:]]' hg38.chrom.sizes | sort -V`
- `examples/data/hg38_extrusion_barriers.bed.xz` was generated as described the Methods of MoDLE's [paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02815-7) (third paragraph of section "MoDLE simulations").
  Refer to code hosted on [github.com/paulsengroup/2021-modle-paper-001-data-analysis](https://github.com/paulsengroup/2021-modle-paper-001-data-analysis) for more details (useful links: [link1](https://github.com/paulsengroup/2021-modle-paper-001-data-analysis/blob/main/workflows/preprocess_data.nf#L41-L100), [link2](https://github.com/paulsengroup/2021-modle-paper-001-data-analysis/blob/main/configs/preprocess_data.config), [link3](https://github.com/paulsengroup/2021-modle-paper-001-data-analysis/blob/main/data/download_list.txt), [link4](https://github.com/paulsengroup/2021-modle-paper-001-data-analysis/blob/main/scripts/convert_chip_signal_to_occupancy.py). Feel free to get in touch if you are struggling to navigate code hosted in the data analysis repository).

</details>

#### Required input files

Running a simulation with default settings only requires two input files:

- A [chrom.sizes file](https://software.broadinstitute.org/software/igv/chromSizes) with the list of chromosomes to be simulated
- A [BED file]((https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bed-format)) with the list of extrusion barriers to include as part of the simulation

Input files using common compression formats (e.g. Bzip2, Gzip, LZ4, xz/LZMA or zstd) are supported.

The extrusion barrier BED file should have at least the first 6 columns defined (i.e. chrom, start, end, name,
score and strand. See [bedtools docs](https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bed-format) for more details).

Some remarks:
- ___chrom___ should contain the name of one of the chromosome defined in the `.chrom.sizes`
- The ___start___ and ___end___ fields are used to position barriers along chromosomes. MoDLE models extrusion barriers as 1bp entities. Thus, when $end - start \gt 1$ the extrusion barrier will be placed in the middle of interval $[start, end)$.
- The ___name___ field is required but ignored, and can be safely set to `.`
- ___score___ is optional and should be set to 0 when not used.
  When non-zero, the ___score___ value will be used to set the average occupancy of the extrusion barrier site defined by the
current line. Refer to [MoDLE's paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02815-7) for more details regarding the extrusion barrier occupancy.
- ___strand___ is required and is used to define the extrusion barrier direction (should be one of `-`, `+` or `.`).
  As of MoDLE v1.0.0, extrusion barriers are modeled after CTCF barriers.
  Thus, the strand field should be populated with the direction of the CTCF binding site that is being defined.
  Barriers without strand information (i.e. with strand `.`) are ignored.

#### Running a simulation with default settings

```console
user@dev:/tmp/modle$ modle simulate \
    --chrom-sizes examples/data/hg38.chrom.sizes \
    --extrusion-barrier-file examples/data/hg38_extrusion_barriers.bed.xz \
    --output-prefix examples/out/hg38_default

[2023-01-06 23:01:24.551] [info]: Running MoDLE v1.0.0
[2023-01-06 23:01:24.552] [info]: Complete log will be written to file "examples/out/hg38_default.log"
[2023-01-06 23:01:24.552] [info]: Writing simulation parameters to config file "examples/out/hg38_default_config.toml"
[2023-01-06 23:01:24.552] [info]: Command: modle simulate --chrom-sizes examples/data/hg38.chrom.sizes --extrusion-barrier-file examples/data/hg38_extrusion_barriers.bed.xz --output-prefix examples/out/hg38_default
[2023-01-06 23:01:24.552] [info]: Simulation will use up to 16 out of 16 available CPU cores.
[2023-01-06 23:01:24.552] [info]: Using --target-contact-density=1.00 as stopping criterion.
[2023-01-06 23:01:24.552] [info]: Contact sampling strategy: tad-plus-loop-with-noise.
[2023-01-06 23:01:24.552] [info]: Contact matrix resolution: 5000bp
[2023-01-06 23:01:24.552] [info]: Importing chromosomes from file "examples/data/hg38.chrom.sizes"...
[2023-01-06 23:01:24.554] [info]: Imported 24 chromosomes in 1.650957ms.
[2023-01-06 23:01:24.554] [info]: Importing extrusion barriers from file "examples/data/hg38_extrusion_barriers.bed.xz"...
[2023-01-06 23:01:24.603] [info]: Imported 38815 barriers in 48.915543ms.
[2023-01-06 23:01:24.606] [info]: Spawning simulation thread 0...
[2023-01-06 23:01:24.606] [info]: Spawning simulation thread 2...
[2023-01-06 23:01:24.606] [info]: Spawning simulation thread 3...
...
[2023-01-06 23:01:24.620] [info]: Begin processing "chr1": simulating ~37485 epochs across 512 cells using 4979 LEFs and 3518 barriers (~73 epochs per cell)...
[2023-01-06 23:01:44.824] [info]: Begin processing "chr2": simulating ~37501 epochs across 512 cells using 4844 LEFs and 2974 barriers (~73 epochs per cell)...
[2023-01-06 23:01:51.086] [info]: Simulation of "chr1" successfully completed.
[2023-01-06 23:01:51.499] [info]: Writing contacts for "chr1" to file "examples/out/hg38_default.cool"...
[2023-01-06 23:01:52.990] [info]: Written 29875200 contacts for "chr1" to file "examples/out/hg38_default.cool" (3.47M nnz out of 29.88M pixels).
[2023-01-06 23:02:04.607] [info]: Begin processing "chr3": simulating ~37474 epochs across 512 cells using 3966 LEFs and 2469 barriers (~73 epochs per cell)...
...
[2023-01-06 23:05:58.431] [info]: Writing contacts for "chrY" to file "examples/out/hg38_default.cool"...
[2023-01-06 23:05:58.810] [info]: Written 6867600 contacts for "chrY" to file "examples/out/hg38_default.cool" (0.92M nnz out of 6.87M pixels).
[2023-01-06 23:06:10.770] [info]: Simulation terminated without errors in 4m46.217517852s!

Bye.
```

The above command will create folder `examples/out/` (if it doesn't already exist), and write the following files inside it:

```
examples/out
├── hg38_default_config.toml
├── hg38_default.cool
└── hg38_default.log
```

Simulated interactions are written to file `hg38_default.cool`. More information about the .cool format is available at [open2c/cooler](https://github.com/open2c/cooler).

In case any of the output files already exist, MoDLE will refuse to run and print an error message listing the file name
collisions.

Passing the `--force` flag overrides this behavior and will cause MoDLE to overwrite existing files.

To visualize `.cool` files we recommend using [cooler show](https://open2c.github.io/cooler/) or [higlass](https://docs.higlass.io/) (more specifically [higlass-docker](https://github.com/higlass/higlass-docker) or [higlass-manage](https://github.com/higlass/higlass-manage)).

<details>
<summary> <b> Tips for visualization </b> </summary>

##### Quickly visualizing interactions with Cooler

```bash
cooler show examples/out/hg38_default.cool chr1:5000000-8000000
```

![chr1 5-8Mbp heatmap](examples/images/hg38_default_001.avif)

##### Visualizing interactions with HiGlass

Before visualizing heatmaps with HiGlass, the single-resolution `.cool` file produced by MoDLE needs to be converted to a multi-resolution Cooler (`.mcool`) file.
This is easily done with `cooler zoomify`:

```console
user@dev:/tmp/modle$ cooler zoomify examples/out/hg38_default.cool

INFO:cooler.cli.zoomify:Recursively aggregating "examples/out/hg38_default.cool"
INFO:cooler.cli.zoomify:Writing to "examples/out/hg38_default.mcool"
INFO:cooler.reduce:Copying base matrices and producing 12 new zoom levels.
INFO:cooler.reduce:Bin size: 5000
INFO:cooler.reduce:Aggregating from 5000 to 10000.
INFO:cooler.create:Creating cooler at "examples/out/hg38_default.mcool::/resolutions/10000"
INFO:cooler.create:Writing chroms
INFO:cooler.create:Writing bins
INFO:cooler.create:Writing pixels
INFO:cooler.reduce:0 10000066
INFO:cooler.reduce:10000066 20000079
INFO:cooler.reduce:20000079 30000056
INFO:cooler.reduce:30000056 40000018
INFO:cooler.reduce:40000018 43586992
INFO:cooler.create:Writing indexes
INFO:cooler.create:Writing info
INFO:cooler.reduce:Aggregating from 10000 to 20000.
...
```

Next, you need to start HiGlass server.

If you are using `higlass-manage`, you may have to prefix commands with `sudo` or `sudo -E`.
```console
# Using higlass-docker
user@dev:/tmp/modle$ sudo docker pull higlass/higlass-docker

Using default tag: latest
latest: Pulling from higlass/higlass-docker
Digest: sha256:44086069ee7d4d3f6f3f0012569789ec138f42b84aa44357826c0b6753eb28de
Status: Image is up to date for higlass/higlass-docker:latest
docker.io/higlass/higlass-docker:latest

user@dev:/tmp/modle$ sudo docker run \
   --detach                 \
   --publish 8989:80        \
   --volume ~/hg-data:/data \
   --volume ~/hg-tmp:/tmp   \
   --name higlass-modle     \
   higlass/higlass-docker

eb3fe07ea23b57eb844a590e7f1e939cf82a11e5d25db87c4b25024a41c9a546

# Using higlass-manage
user@dev:/tmp/modle$ sudo higlass-manage start --use-redis --no-public-data

Stopping previously running container
Attempting to remove existing Docker network instance
Stopping previously running Redis container
Pulling redis:5.0.3-alpine
done
Pulling latest image...
done
Data directory: /home/user/hg-data
Temp directory: ()
Starting... default 8989
Docker started: higlass-manage-container-default
sending request 1
Waiting to start (tilesets)...
sending request 2
...
public_data: False
ret: b'{"uid": "default_local"}'
ret: ExecResult(exit_code=0, output=b'')
Replaced js file
Started
```

Now if you visit [localhost:8989/app](http://localhost:8989/app) with your browser you should see an empty HiGlass page like the following:

![higlass empty page](examples/images/higlass_001.avif)

Next, we need to ingest our `.mcool` into HiGlass:

```console
user@dev:/tmp/modle$ mkdir ~/hg-data/media/modle
user@dev:/tmp/modle$ cp examples/out/hg38_default.mcool ~/hg-data/media/modle

# With higlass-docker
user@dev:/tmp/modle$ sudo docker exec higlass-modle \
  python higlass-server/manage.py ingest_tileset \
    --project-name "modle" \
    --filename modle/hg38_default.mcool \
    --filetype cooler \
    --datatype matrix \
    --no-upload

# With higlass-manage
user@dev:/tmp/modle$ sudo higlass-manage ingest --project-name=modle modle/hg38_default.mcool --no-upload
state True
Inferred filetype: cooler
Inferred datatype: matrix
state True
name_text:
hg_name: default
command: python higlass-server/manage.py ingest_tileset --filename modle/hg38_default.mcool --filetype cooler --datatype matrix  --project-name "modle"  --no-upload --uid emFfbA9tTCiBTt6b5vTG8Q
```

Finally, after refreshing page [localhost:8989/app](http://localhost:8989/app) you should be able to add tileset `hg38_default.mcool` to the center panel:

![higlass default settings](examples/images/higlass_002.avif)

In case you are struggling to visualize MoDLE's output feel free to reach out by starting a new [discussion](https://github.com/paulsengroup/modle/discussions).

</details>

<details>
<summary> <b> Tips for programmatic access to Cooler files </b> </summary>

Programmatic access to contacts stored in Cooler files is possible using the `cooler` [CLI interface](https://cooler.readthedocs.io/en/latest/cli.html) or Python [API](https://cooler.readthedocs.io/en/latest/api.html).

##### CLI interface

```console
# Dumping interactions in BEDPE format
user@dev:/tmp/modle$ cooler dump --table pixels --join examples/out/hg38_default.cool | head

chr1	0	5000	chr1	0	5000	9
chr1	0	5000	chr1	5000	10000	5
chr1	0	5000	chr1	10000	15000	2
chr1	0	5000	chr1	15000	20000	2
chr1	0	5000	chr1	20000	25000	1
chr1	0	5000	chr1	25000	30000	2
chr1	0	5000	chr1	30000	35000	3
chr1	0	5000	chr1	35000	40000	1
chr1	0	5000	chr1	40000	45000	2
chr1	0	5000	chr1	45000	50000	1

user@dev:/tmp/modle$ cooler dump --table pixels --join examples/out/hg38_default.mcool::/resolutions/200000 | head

chr1	0	20000	chr1	0	20000	114
chr1	0	20000	chr1	20000	40000	120
chr1	0	20000	chr1	40000	60000	69
chr1	0	20000	chr1	60000	80000	59
chr1	0	20000	chr1	80000	100000	34
chr1	0	20000	chr1	100000	120000	28
chr1	0	20000	chr1	120000	140000	25
chr1	0	20000	chr1	140000	160000	18
chr1	0	20000	chr1	160000	180000	20
chr1	0	20000	chr1	180000	200000	12
```

##### Python API

```python-console
In [1]: import cooler
In [2]: c = cooler.Cooler("examples/out/hg38_default.cool")

In [3]: selector = c.matrix(balance=False)
In [4]: selector.fetch("chr1:0-50000")
Out[4]:
array([[ 9,  5,  2,  2,  1,  2,  3,  1,  2,  1],
       [ 5,  9,  9, 15,  3,  3,  8,  1,  8,  1],
       [ 2,  9, 12, 26, 11,  7,  3,  8,  5,  2],
       [ 2, 15, 26, 25, 27, 16, 12, 14, 10,  7],
       [ 1,  3, 11, 27, 20, 42, 30, 17, 15, 10],
       [ 2,  3,  7, 16, 42, 21, 27, 17, 15, 12],
       [ 3,  8,  3, 12, 30, 27, 26, 27, 26, 17],
       [ 1,  1,  8, 14, 17, 17, 27, 23, 44, 23],
       [ 2,  8,  5, 10, 15, 15, 26, 44, 31, 35],
       [ 1,  1,  2,  7, 10, 12, 17, 23, 35, 23]], dtype=int32)

In [5]: selector = c.matrix(balance=False, as_pixels=True, join=True)
In [6]: selector.fetch("chr1:0-50000")
Out[6]:
   chrom1  start1   end1 chrom2  start2   end2  count
0    chr1       0   5000   chr1       0   5000      9
1    chr1       0   5000   chr1    5000  10000      5
2    chr1       0   5000   chr1   10000  15000      2
3    chr1       0   5000   chr1   15000  20000      2
4    chr1       0   5000   chr1   20000  25000      1
5    chr1       0   5000   chr1   25000  30000      2
6    chr1       0   5000   chr1   30000  35000      3
7    chr1       0   5000   chr1   35000  40000      1
8    chr1       0   5000   chr1   40000  45000      2
9    chr1       0   5000   chr1   45000  50000      1
...
```

Refer to Cooler [documentation](https://cooler.readthedocs.io/en/latest/quickstart.html) for more details.

</details>

#### Running a simulation using config files

The `hg38_default_config.toml` file produced by MoDLE in the [previous section](#running-a-simulation-with-default-settings) is a config file that can be used to repeat simulations using the same set of parameters and input files.
Information found in config files is also embedded in the `.cool` file as a JSON string, and can be extracted using `cooler info`:
<!-- TODO: Update once https://github.com/paulsengroup/modle/issues/75 is implemented -->
```console
user@dev:/tmp/modle$ cooler info --metadata examples/out/hg38_default.cool
{
    "simulate": {
        "avg-lef-processivity": 300000,
        "burnin-extr-speed-coefficient": 1,
        "burnin-history-length": 100,
        "burnin-smoothing-window-size": 5,
        "burnin-target-epochs-for-lef-activation": 320,
        "chrom-sizes": "examples/data/hg38.chrom.sizes",
...
```

The following will run a simulation with the same settings as the simulation from the [previous example](#Running-a-simulation-with-default-settings)

```bash
modle simulate --config examples/out/hg38_default_config.toml --force
```

Parameters read from a config file have lower priority than parameter specified through the CLI,
meaning that CLI options can be used to override parameters read from the config file.

This can be useful when running a batch of simulations using a large number of optional parameters, and where only a
handful of parameters need to change across simulation runs.

The following command will run a simulation using parameters from the [previous example](#Running-a-simulation-with-default-settings) as starting point, but using a different number of cells and target contact density.

```bash
modle simulate \
    --config examples/out/hg38_default_config.toml \
    --ncells 16 \
    --target-contact-density 10 \
    --output-prefix examples/out/hg38_deep
```

![Simulation w/ high contact density](examples/images/higlass_003.avif)

MoDLE's config files are in [TOML format](https://toml.io/en/).

Adding a line like `my-option=my_value` to the config file has the same effect as passing `--my-option=my_value` on the CLI.

For an up-to-date list of supported CLI options, please refer to MoDLE's help message:

```console
user@dev:/tmp/modle$ modle simulate --help

Simulate loop extrusion and write resulting molecular contacts in a .cool file.
Usage: /usr/local/modle simulate [OPTIONS]

Options:
  -h,--help                   Print this help message and exit
[Option Group: IO]
  Options controlling MoDLE input, output and logs.
  Options:
    -c,--chrom-sizes REQUIRED   Path to file with chromosome sizes in chrom.sizes format.
    -b,--extrusion-barrier-file REQUIRED
                                Path to a file in BED6+ format with the genomic coordinates of extrusion barriers to be
                                simulated. The score field in a BED record should be a number between 0 and 1 and is
                                interpreted as the extrusion barrier occupancy for the extrusion barrier described
                                by the record.
                                Barriers mapping on chromosomes not listed in the chrom.sizes file passed through
                                the --chrom-sizes option are ignored.
    -f,--force                  Overwrite existing files (if any).
    -o,--output-prefix REQUIRED Output prefix.
                                Can be an absolute or relative path including the file name but without the extension.
                                Example: running modle sim -o /tmp/my_simulation ... yields the following files:
                                         - /tmp/my_simulation.cool
                                         - /tmp/my_simulation.log
                                         - /tmp/my_simulation_config.toml
    --assembly-name=unknown     Name of the genome assembly to be simulated.
                                This is only used to populate the "assembly" attribute in the output .cool file.
...
```

## Citing

If you use MoDLE in your research, please cite the following paper:

Rossini, R., Kumar, V., Mathelier, A. _et al._ MoDLE: high-performance stochastic modeling of DNA loop extrusion interactions. _Genome Biol_ __23__, 247 (2022). [https://doi.org/10.1186/s13059-022-02815-7](https://doi.org/10.1186/s13059-022-02815-7)


<details>
<summary>BibTex</summary>

```bibtex
@article{rossini_2022,
author = {Rossini, Roberto and Kumar, Vipin and Mathelier, Anthony and Rognes, Torbjørn and Paulsen, Jonas},
doi = {10.1186/s13059-022-02815-7},
journal = {Genome Biology},
month = {11},
title = {{MoDLE: high-performance stochastic modeling of DNA loop extrusion interactions}},
url = {https://doi.org/10.1186/s13059-022-02815-7},
volume = {23},
year = {2022}
}
```

</details>

In addition to the above reference, you can also cite specific version(s) of MoDLE using one of the DOIs from Zenodo: [10.5281/zenodo.6424697](https://doi.org/10.5281/zenodo.6424697)
