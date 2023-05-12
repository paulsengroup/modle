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
user@dev:/tmp$ docker run ghcr.io/paulsengroup/modle:1.0.1 --help

# Using Singularity/Apptainer
user@dev:/tmp$ singularity run ghcr.io/paulsengroup/modle:1.0.1 --help

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
docker run ghcr.io/paulsengroup/modle:1.0.1 simulate \
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
  ghcr.io/paulsengroup/modle:1.0.1 simulate \
  --chrom-sizes /data/grch38.chrom.sizes \
  --extrusion-barrier-file /data/grch38_h1_extrusion_barriers.bed.xz \
  --output-prefix /data/output

# Mehtod 2: use two different folders for input/output. Input folder is mounted in read-only mode
docker run -v "$(pwd -P)/mydata/:/input_data/:ro" \
           -v "$(pwd -P)/output/:/output_data/" \
  ghcr.io/paulsengroup/modle:1.0.1 simulate \
  --chrom-sizes /input_data/grch38.chrom.sizes \
  --extrusion-barrier-file /input_data/grch38_h1_extrusion_barriers.bed.xz \
  --output-prefix /output_data/output

# Methods for Singularity/Apptainer (note that manually mounting paths is usually not required when using Singularity/Apptainer)
singularity run -B "$(pwd -P)/mydata/:/data/" \
  docker://ghcr.io/paulsengroup/modle:1.0.1 simulate \
  --chrom-sizes /data/grch38.chrom.sizes \
  --extrusion-barrier-file /data/grch38_h1_extrusion_barriers.bed.xz \
  --output-prefix /data/output

singularity run -B "$(pwd -P)/mydata/:/input_data/:ro" \
                -B "$(pwd -P)/output/:/output_data/" \
  docker://ghcr.io/paulsengroup/modle:1.0.1 simulate \
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
-rw-r--r--  1 root root  2.4M Jan  5 12:05 output_lef_1d_occupancy.bw
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
-rw-r--r--  1 root root  2.4M Jan  5 12:05 output_lef_1d_occupancy.bw
-rw-r--r--  1 root root   57M Jan  5 12:05 output.cool
-rw-r--r--  1 root root   14K Jan  5 12:05 output.log
```

__Troubleshooting problems inside the container is a bit tedious. Any tips?__

Yes! For troubleshooting we recommend using an interactive shell running inside the container.

```bash
docker run -v "$(pwd -P)/mydata/:/input_data/:ro" \
           -v "$(pwd -P)/output/:/output_data/" \
           -it --rm --entrypoint=/bin/bash \
  ghcr.io/paulsengroup/modle:1.0.1
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
singularity shell docker://ghcr.io/paulsengroup/modle:1.0.1
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
MoDLE-v1.0.1-bioconda
```

### Pre-built Linux binaries (not recommended)

Download file `modle-1.0.1-x86_64-linux.tar.xz` from the [latest release](https://github.com/paulsengroup/modle/releases/latest).

Extract the archive and copy `modle` and/or `modle_tools` to a location in your `PATH`

<details>
<summary> click to expand </summary>

<!-- TODO make sure the link is correct once the release goes live -->
```console
user@dev:/tmp$ curl -LO 'https://github.com/paulsengroup/modle/releases/download/v1.0.1/modle-1.0.1-x86_64-linux.tar.xz'

user@dev:/tmp$ tar -xvf modle-1.0.1-x86_64-linux.tar.xz
modle-1.0.1-x86_64-linux/
modle-1.0.1-x86_64-linux/share/
modle-1.0.1-x86_64-linux/share/licenses/
modle-1.0.1-x86_64-linux/share/licenses/modle_tools/
modle-1.0.1-x86_64-linux/share/licenses/modle_tools/LICENSE
modle-1.0.1-x86_64-linux/share/licenses/modle/
modle-1.0.1-x86_64-linux/share/licenses/modle/LICENSE
modle-1.0.1-x86_64-linux/bin/
modle-1.0.1-x86_64-linux/bin/modle_tools
modle-1.0.1-x86_64-linux/bin/modle

user@dev:/tmp modle-1.0.1-x86_64-linux/bin/modle --version
MoDLE-v1.0.1

# Optional: add modle and modle_tools to your PATH (please ensure $HOME/.local/bin/ is in your PATH)
user@dev:/tmp$ install -Dm0755 modle-1.0.1-x86_64-linux/bin/modle "$HOME/.local/bin/"
user@dev:/tmp$ install -Dm0755 modle-1.0.1-x86_64-linux/bin/modle_tools "$HOME/.local/bin/"

user@dev:/tmp$ whereis modle
modle: /home/user/.local/bin/modle

user@dev:/tmp$ modle --version
MoDLE-v1.0.1
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
curl -L 'https://github.com/paulsengroup/modle/archive/refs/tags/v1.0.1.tar.gz' | tar --strip-components=1 -C /tmp/modle -xzf -
```

Using git.
```bash
git clone https://github.com/paulsengroup/modle.git /tmp/modle

cd /tmp/modle
git checkout v1.0.1  # Skip this step if you want to build the latest commit from main
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
125/128 Test #120: ContactMatrix internal: encode/decode roundtrip - LONG ...........................................   Passed   35.05 sec
126/128 Test #117: ContactMatrixSparse serialization - LONG .........................................................   Passed   42.19 sec
127/128 Test #119: ContactMatrixSparse serialization (FP) - LONG ....................................................   Passed   44.71 sec
128/128 Test #118: ContactMatrixSparse serialization (destructive) - LONG ...........................................   Passed   32.53 sec

100% tests passed, 0 tests failed out of 128

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

- cooler>=0.9
- xz
- common UNIX shell commands such as (namely `cmp`, `diff`, `grep` and `mktemp`)

cooler can be installed using pip:
```bash
/tmp/venv/bin/pip3 install 'cooler>=0.9'
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
[2023-05-11 11:42:03.315] [info]: running MoDLE v1.0.1
[2023-05-11 11:42:03.315] [info]: complete log will be written to file "/tmp/modle-ci-SpASy09Lp3/out.log"
[2023-05-11 11:42:03.315] [info]: writing simulation parameters to config file "/tmp/modle-ci-SpASy09Lp3/out_config.toml"
[2023-05-11 11:42:03.316] [info]: command: cmake-build-release-llvm-test-release/src/modle/modle sim -c modle/test/data/integration_tests/grch38.chrom.sizes -g modle/test/data/integration_tests/grch38_regions_of_interest.bed -b modle/test/data/integration_tests/grch38_h1_extrusion_barriers.bed.xz -o /tmp/modle-ci-SpASy09Lp3/out -r 20kb --verbose --target-contact-density 20 --ncells 2 --track-1d-lef-position --max-burnin-epochs 5000
[2023-05-11 11:42:03.316] [info]: simulation will use up to 16 out of 16 available CPU cores.
[2023-05-11 11:42:03.316] [info]: using --target-contact-density=20.00 as stopping criterion.
[2023-05-11 11:42:03.316] [info]: contact sampling strategy: tad-plus-loop-with-noise.
[2023-05-11 11:42:03.316] [info]: contact matrix resolution: 20000bp
[2023-05-11 11:42:03.316] [info]: importing chromosomes from "modle/test/data/integration_tests/grch38.chrom.sizes"...
[2023-05-11 11:42:03.317] [info]: imported 24 chromosomes in 697.651us.
[2023-05-11 11:42:03.317] [info]: importing genomic intervals from "modle/test/data/integration_tests/grch38_regions_of_interest.bed"...
[2023-05-11 11:42:03.317] [info]: imported 11 intervals in 319.171us.
[2023-05-11 11:42:03.317] [info]: importing extrusion barriers from "modle/test/data/integration_tests/grch38_h1_extrusion_barriers.bed.xz"...
[2023-05-11 11:42:03.361] [info]: imported 2087 barriers from "modle/test/data/integration_tests/grch38_h1_extrusion_barriers.bed.xz" in 43.849816ms.
[2023-05-11 11:42:03.362] [warning]: simulated size for the following 5 interval(s) is smaller than the simulation diagonal width (3000000 bp). Is this intended?
 - chr17:40200000-42800000: 2600000;
 - chr17:46800000-49300000: 2500000;
 - chr17:60200000-63100000: 2900000;
 - chr17:64600000-66200000: 1600000;
 - chr17:76800000-77200000: 400000;
[2023-05-11 11:42:03.362] [warning]: contact matrix for the following 10 interval(s) appears to be really small (less than 250000 pixels). Is this intended?
 - chr17:3400000-6500000: 23250 pixels
 - chr17:10800000-16100000: 39750 pixels
 - chr17:33500000-39800000: 47250 pixels
 - chr17:40200000-42800000: 19500 pixels
 - chr17:46800000-49300000: 18750 pixels
 - chr17:52100000-59500000: 55500 pixels
 - chr17:60200000-63100000: 21750 pixels
 - chr17:64600000-66200000: 12000 pixels
 - chr17:69100000-72900000: 28500 pixels
 - chr17:76800000-77200000: 3000 pixels
[2023-05-11 11:42:03.362] [info]: spawning IO thread...
[2023-05-11 11:42:03.362] [info]: spawning worker thread W0...
[2023-05-11 11:42:03.362] [info]: spawning worker thread W9...
[2023-05-11 11:42:03.362] [info]: spawning worker thread W2...
[2023-05-11 11:42:03.362] [debug]: [main]: submitting task #0 (chr17:3400000-6500000 cell #0)...
[2023-05-11 11:42:03.362] [info]: spawning worker thread W15...
[2023-05-11 11:42:03.362] [info]: spawning worker thread W4...
[2023-05-11 11:42:03.362] [info]: spawning worker thread W6...
[2023-05-11 11:42:03.362] [info]: spawning worker thread W7...
[2023-05-11 11:42:03.362] [info]: spawning worker thread W8...
[2023-05-11 11:42:03.362] [info]: spawning worker thread W11...
[2023-05-11 11:42:03.362] [info]: spawning worker thread W10...
[2023-05-11 11:42:03.362] [info]: spawning worker thread W13...
[2023-05-11 11:42:03.362] [info]: spawning worker thread W12...
[2023-05-11 11:42:03.362] [info]: spawning worker thread W14...
[2023-05-11 11:42:03.362] [info]: spawning worker thread W1...
[2023-05-11 11:42:03.362] [info]: spawning worker thread W3...
[2023-05-11 11:42:03.362] [debug]: [main]: submitting task #1 (chr17:3400000-6500000 cell #1)...
[2023-05-11 11:42:03.363] [debug]: [main]: submitting task #2 (chr17:10800000-16100000 cell #0)...
[2023-05-11 11:42:03.363] [debug]: [main]: submitting task #3 (chr17:10800000-16100000 cell #1)...
[2023-05-11 11:42:03.363] [debug]: [main]: submitting task #4 (chr17:33500000-39800000 cell #0)...
[2023-05-11 11:42:03.362] [info]: spawning worker thread W5...
[2023-05-11 11:42:03.363] [debug]: [main]: submitting task #5 (chr17:33500000-39800000 cell #1)...
[2023-05-11 11:42:03.363] [debug]: [main]: submitting task #6 (chr17:40200000-42800000 cell #0)...
[2023-05-11 11:42:03.363] [debug]: [main]: submitting task #7 (chr17:40200000-42800000 cell #1)...
[2023-05-11 11:42:03.363] [debug]: [main]: submitting task #8 (chr17:46800000-49300000 cell #0)...
[2023-05-11 11:42:03.363] [debug]: [main]: submitting task #9 (chr17:46800000-49300000 cell #1)...
[2023-05-11 11:42:03.363] [debug]: [main]: submitting task #10 (chr17:52100000-59500000 cell #0)...
[2023-05-11 11:42:03.363] [debug]: [main]: submitting task #11 (chr17:52100000-59500000 cell #1)...
[2023-05-11 11:42:03.363] [debug]: [main]: submitting task #12 (chr17:60200000-63100000 cell #0)...
[2023-05-11 11:42:03.363] [debug]: [main]: submitting task #13 (chr17:60200000-63100000 cell #1)...
[2023-05-11 11:42:03.363] [debug]: [main]: submitting task #14 (chr17:64600000-66200000 cell #0)...
[2023-05-11 11:42:03.363] [debug]: [main]: submitting task #15 (chr17:64600000-66200000 cell #1)...
[2023-05-11 11:42:03.363] [debug]: [main]: submitting task #16 (chr17:69100000-72900000 cell #0)...
[2023-05-11 11:42:03.363] [debug]: [main]: submitting task #17 (chr17:69100000-72900000 cell #1)...
[2023-05-11 11:42:03.363] [debug]: [main]: submitting task #18 (chr17:76800000-77200000 cell #0)...
[2023-05-11 11:42:03.363] [debug]: [main]: submitting task #19 (chr17:76800000-77200000 cell #1)...
[2023-05-11 11:42:03.363] [debug]: [W7]: begin processing task 11 (chr17:52100000-59500000 cell #1, 0x7dba39cdd272b1bf83640f92d867800050ed7654ab1cd71def39ddaa685862a4)...
[2023-05-11 11:42:03.363] [debug]: [W0]: begin processing task 1 (chr17:3400000-6500000 cell #1, 0xf989845f1d5e7366868d488424e5721f7a8c93f84c9e969febf7e8207c97060c)...
[2023-05-11 11:42:03.363] [info]: begin processing chr17:33500000-39800000: simulating ~11667 epochs across 2 cells using 126 LEFs and 137 barriers (~5833 epochs per cell)...
[2023-05-11 11:42:03.363] [debug]: [W12]: begin processing task 4 (chr17:33500000-39800000 cell #0, 0xacf52581c195c57e39e98ddbe38add37f6549f1ee8b840846975a121751c865b)...
[2023-05-11 11:42:03.363] [debug]: [W14]: begin processing task 3 (chr17:10800000-16100000 cell #1, 0x75fe0e29cac63fafb2f7295e5d2ada8143d515f87bdc6b45a329fcc1de2bad76)...
[2023-05-11 11:42:03.363] [info]: begin processing chr17:52100000-59500000: simulating ~11684 epochs across 2 cells using 148 LEFs and 106 barriers (~5842 epochs per cell)...
[2023-05-11 11:42:03.363] [debug]: [W2]: begin processing task 10 (chr17:52100000-59500000 cell #0, 0xc8352a6614980a95ad7d7142970c9e77216ea23c47e035be5eb926e3f770e34b)...
[2023-05-11 11:42:03.363] [debug]: [W10]: begin processing task 7 (chr17:40200000-42800000 cell #1, 0x77c8bde3ddf1faeeb8f5c4b09ed93c30392bfffe57974f1528f9d51d3eb7ad4d)...
[2023-05-11 11:42:03.363] [info]: begin processing chr17:10800000-16100000: simulating ~11691 epochs across 2 cells using 106 LEFs and 67 barriers (~5845 epochs per cell)...
[2023-05-11 11:42:03.363] [debug]: [W9]: begin processing task 2 (chr17:10800000-16100000 cell #0, 0x855f4ecc138b144575074359f6ee20516a35db461ce6a6986508bbafe169de0d)...
[2023-05-11 11:42:03.363] [debug]: [W1]: begin processing task 9 (chr17:46800000-49300000 cell #1, 0x35f79d9403df23f29508e06f1e7c8576d67f1e279c958778bf4da4a3b577f069)...
[2023-05-11 11:42:03.363] [info]: begin processing chr17:3400000-6500000: simulating ~11625 epochs across 2 cells using 62 LEFs and 80 barriers (~5812 epochs per cell)...
[2023-05-11 11:42:03.363] [debug]: [W3]: begin processing task 13 (chr17:60200000-63100000 cell #1, 0xa1f8b81598ef32e1bb5f939da04508662e1824121ae114051199c032087a6627)...
[2023-05-11 11:42:03.363] [info]: begin processing chr17:60200000-63100000: simulating ~11757 epochs across 2 cells using 58 LEFs and 35 barriers (~5878 epochs per cell)...
[2023-05-11 11:42:03.363] [debug]: [W6]: begin processing task 5 (chr17:33500000-39800000 cell #1, 0x635ba9dffc765bb40bf4cebcb302a6b3436beaf044a7fc9200f03ee5f9c538f0)...
[2023-05-11 11:42:03.363] [info]: begin processing chr17:64600000-66200000: simulating ~12000 epochs across 2 cells using 32 LEFs and 24 barriers (~6000 epochs per cell)...
[2023-05-11 11:42:03.363] [debug]: [W5]: begin processing task 15 (chr17:64600000-66200000 cell #1, 0x37aa169907199d414d7b15425a3a5cc6d5f7a22caebf44178323862575834445)...
[2023-05-11 11:42:03.363] [info]: begin processing chr17:40200000-42800000: simulating ~11818 epochs across 2 cells using 52 LEFs and 77 barriers (~5909 epochs per cell)...
[2023-05-11 11:42:03.363] [info]: begin processing chr17:46800000-49300000: simulating ~11719 epochs across 2 cells using 50 LEFs and 65 barriers (~5859 epochs per cell)...
[2023-05-11 11:42:03.363] [debug]: [W11]: begin processing task 8 (chr17:46800000-49300000 cell #0, 0xf9f36a5d81f9efbe2e7e9007e7e62b88d2483867a3636e959d992113c304e749)...
[2023-05-11 11:42:03.364] [debug]: allocating a 150x80 contact matrix...
[2023-05-11 11:42:03.363] [debug]: [main]: submitting task #20 (chr19:0-58617616 cell #0)...
[2023-05-11 11:42:03.363] [debug]: [W13]: begin processing task 12 (chr17:60200000-63100000 cell #0, 0x553b7b0c46f1c8bb42f6a71bffabe0bb71fc30ace94b900f4db4285c550f184d)...
[2023-05-11 11:42:03.364] [debug]: allocating a vector of size 80...
[2023-05-11 11:42:03.365] [debug]: allocating a 150x130 contact matrix...
[2023-05-11 11:42:03.365] [debug]: allocating a vector of size 130...
[2023-05-11 11:42:03.365] [debug]: allocating a 150x265 contact matrix...
[2023-05-11 11:42:03.363] [debug]: [W8]: begin processing task 14 (chr17:64600000-66200000 cell #0, 0x2c0d706739eb7bf9687d641e3c1716c96368835d83688ca6e31758f0fa227b26)...
[2023-05-11 11:42:03.365] [debug]: [main]: submitting task #21 (chr19:0-58617616 cell #1)...
[2023-05-11 11:42:03.365] [debug]: waiting for worker threads to return...
[2023-05-11 11:42:03.364] [debug]: [W4]: begin processing task 6 (chr17:40200000-42800000 cell #0, 0x6c4136a584f59c48aca20dcabbe48501963aad40b7fef2445473d046bb24f948)...
[2023-05-11 11:42:03.366] [debug]: allocating a vector of size 265...
[2023-05-11 11:42:03.366] [debug]: allocating a 150x145 contact matrix...
[2023-05-11 11:42:03.366] [debug]: allocating a vector of size 145...
[2023-05-11 11:42:03.367] [debug]: allocating a 150x155 contact matrix...
[2023-05-11 11:42:03.367] [debug]: allocating a vector of size 155...
[2023-05-11 11:42:03.368] [debug]: allocating a 150x125 contact matrix...
[2023-05-11 11:42:03.368] [debug]: allocating a vector of size 125...
[2023-05-11 11:42:03.368] [debug]: allocating a 150x315 contact matrix...
[2023-05-11 11:42:03.369] [debug]: allocating a vector of size 315...
[2023-05-11 11:42:03.373] [debug]: allocating a 150x370 contact matrix...
[2023-05-11 11:42:03.363] [debug]: [W15]: begin processing task 0 (chr17:3400000-6500000 cell #0, 0x6e705dd633f92ffa562fa5b8ad800ff4b323454db6b95b1434b52875d2961d9e)...
[2023-05-11 11:42:03.380] [debug]: allocating a vector of size 370...
[2023-05-11 11:42:03.468] [debug]: [W5]: finished processing task 15 (chr17:64600000-66200000 cell #1, 0x37aa169907199d414d7b15425a3a5cc6d5f7a22caebf44178323862575834445): collected 120000 interactions throughout 6851 epochs (148 burnin epochs)
[2023-05-11 11:42:03.468] [info]: begin processing chr17:69100000-72900000: simulating ~11633 epochs across 2 cells using 76 LEFs and 43 barriers (~5816 epochs per cell)...
[2023-05-11 11:42:03.468] [debug]: [W5]: begin processing task 16 (chr17:69100000-72900000 cell #0, 0xc35a95ef10d1fcb0bb6c84546f3b5c771ef2d412b3e39efe52a284b74dc37fbc)...
[2023-05-11 11:42:03.472] [debug]: [W8]: finished processing task 14 (chr17:64600000-66200000 cell #0, 0x2c0d706739eb7bf9687d641e3c1716c96368835d83688ca6e31758f0fa227b26): collected 120000 interactions throughout 7044 epochs (351 burnin epochs)
[2023-05-11 11:42:03.473] [debug]: [W8]: begin processing task 17 (chr17:69100000-72900000 cell #1, 0x3cfd0646005565a3616be293793ae8b47a183af2bd26de648e58238985c59fdd)...
[2023-05-11 11:42:03.477] [debug]: allocating a 150x190 contact matrix...
[2023-05-11 11:42:03.477] [debug]: allocating a vector of size 190...
[2023-05-11 11:42:03.524] [debug]: [W1]: finished processing task 9 (chr17:46800000-49300000 cell #1, 0x35f79d9403df23f29508e06f1e7c8576d67f1e279c958778bf4da4a3b577f069): collected 187500 interactions throughout 6531 epochs (359 burnin epochs)
[2023-05-11 11:42:03.524] [info]: begin processing chr17:76800000-77200000: simulating ~12000 epochs across 2 cells using 8 LEFs and 12 barriers (~6000 epochs per cell)...
[2023-05-11 11:42:03.524] [debug]: [W1]: begin processing task 18 (chr17:76800000-77200000 cell #0, 0x01bc4ddbec98829a30ac06832dd289e4378615b4cc7e1f7d4a52a8360d484329)...
[2023-05-11 11:42:03.525] [debug]: allocating a 150x20 contact matrix...
[2023-05-11 11:42:03.525] [debug]: allocating a vector of size 20...
[2023-05-11 11:42:03.529] [debug]: [W11]: finished processing task 8 (chr17:46800000-49300000 cell #0, 0xf9f36a5d81f9efbe2e7e9007e7e62b88d2483867a3636e959d992113c304e749): collected 187500 interactions throughout 6582 epochs (416 burnin epochs)
[2023-05-11 11:42:03.530] [debug]: [W11]: begin processing task 19 (chr17:76800000-77200000 cell #1, 0x90d2ff862cd23ec1073aee34c927b088d68a33877b76c3baca2d97868c1fc137)...
[2023-05-11 11:42:03.541] [debug]: [W13]: finished processing task 12 (chr17:60200000-63100000 cell #0, 0x553b7b0c46f1c8bb42f6a71bffabe0bb71fc30ace94b900f4db4285c550f184d): collected 217500 interactions throughout 6514 epochs (187 burnin epochs)
[2023-05-11 11:42:03.541] [info]: begin processing chr19:0-58617616: simulating ~11724 epochs across 2 cells using 1172 LEFs and 1441 barriers (~5862 epochs per cell)...
[2023-05-11 11:42:03.541] [debug]: [W13]: begin processing task 20 (chr19:0-58617616 cell #0, 0x9b93afbc98c209686dc41783d8e0039a0d048f6edd9719becab4bed89163951a)...
[2023-05-11 11:42:03.543] [debug]: [W3]: finished processing task 13 (chr17:60200000-63100000 cell #1, 0xa1f8b81598ef32e1bb5f939da04508662e1824121ae114051199c032087a6627): collected 217500 interactions throughout 6755 epochs (419 burnin epochs)
[2023-05-11 11:42:03.543] [debug]: [W3]: begin processing task 21 (chr19:0-58617616 cell #1, 0xbf6bff8c51e1ce920d4da20970373f04dc2177a111e26e0be5454e2985187a6f)...
[2023-05-11 11:42:03.550] [debug]: [W4]: finished processing task 6 (chr17:40200000-42800000 cell #0, 0x6c4136a584f59c48aca20dcabbe48501963aad40b7fef2445473d046bb24f948): collected 195000 interactions throughout 6510 epochs (230 burnin epochs)
[2023-05-11 11:42:03.552] [debug]: [W10]: finished processing task 7 (chr17:40200000-42800000 cell #1, 0x77c8bde3ddf1faeeb8f5c4b09ed93c30392bfffe57974f1528f9d51d3eb7ad4d): collected 195000 interactions throughout 6435 epochs (165 burnin epochs)
[2023-05-11 11:42:03.561] [debug]: [W1]: finished processing task 18 (chr17:76800000-77200000 cell #0, 0x01bc4ddbec98829a30ac06832dd289e4378615b4cc7e1f7d4a52a8360d484329): collected 30000 interactions throughout 9285 epochs (226 burnin epochs)
[2023-05-11 11:42:03.569] [debug]: [W11]: finished processing task 19 (chr17:76800000-77200000 cell #1, 0x90d2ff862cd23ec1073aee34c927b088d68a33877b76c3baca2d97868c1fc137): collected 30000 interactions throughout 9395 epochs (370 burnin epochs)
[2023-05-11 11:42:03.574] [debug]: [W0]: finished processing task 1 (chr17:3400000-6500000 cell #1, 0xf989845f1d5e7366868d488424e5721f7a8c93f84c9e969febf7e8207c97060c): collected 232500 interactions throughout 6361 epochs (263 burnin epochs)
[2023-05-11 11:42:03.581] [debug]: allocating a 150x2931 contact matrix...
[2023-05-11 11:42:03.583] [debug]: allocating a vector of size 2931...
[2023-05-11 11:42:03.607] [debug]: [W15]: finished processing task 0 (chr17:3400000-6500000 cell #0, 0x6e705dd633f92ffa562fa5b8ad800ff4b323454db6b95b1434b52875d2961d9e): collected 232500 interactions throughout 6429 epochs (334 burnin epochs)
[2023-05-11 11:42:03.608] [info]: [io] writing contacts for chr17:3400000-6500000 to file "/tmp/modle-ci-SpASy09Lp3/out.cool"...
[2023-05-11 11:42:03.630] [info]: [io]: written 465000 contacts for chr17:3400000-6500000 to file "/tmp/modle-ci-SpASy09Lp3/out.cool" in 22.840537ms (0.00M nnz out of 0.02M pixels).
[2023-05-11 11:42:03.630] [info]: [io]: writing 1D LEF occupancy profile for chr17:3400000-6500000 to file "/tmp/modle-ci-SpASy09Lp3/out_lef_1d_occupancy.bw"...
[2023-05-11 11:42:03.630] [info]: [io]: writing 1D LEF occupancy profile for chr17:3400000-6500000 to file "/tmp/modle-ci-SpASy09Lp3/out_lef_1d_occupancy.bw" took 29.038us
[2023-05-11 11:42:03.630] [debug]: deallocating a 150x155 contact matrix...
[2023-05-11 11:42:03.630] [debug]: deallocating a vector of size 155...
[2023-05-11 11:42:03.651] [debug]: [W4]: all tasks have been processed: returning!
[2023-05-11 11:42:03.654] [debug]: [W10]: all tasks have been processed: returning!
[2023-05-11 11:42:03.661] [debug]: [W1]: all tasks have been processed: returning!
[2023-05-11 11:42:03.669] [debug]: [W11]: all tasks have been processed: returning!
[2023-05-11 11:42:03.674] [debug]: [W0]: all tasks have been processed: returning!
[2023-05-11 11:42:03.683] [debug]: [W14]: finished processing task 3 (chr17:10800000-16100000 cell #1, 0x75fe0e29cac63fafb2f7295e5d2ada8143d515f87bdc6b45a329fcc1de2bad76): collected 397500 interactions throughout 6194 epochs (154 burnin epochs)
[2023-05-11 11:42:03.708] [debug]: [W15]: all tasks have been processed: returning!
[2023-05-11 11:42:03.717] [debug]: [W8]: finished processing task 17 (chr17:69100000-72900000 cell #1, 0x3cfd0646005565a3616be293793ae8b47a183af2bd26de648e58238985c59fdd): collected 285000 interactions throughout 6380 epochs (240 burnin epochs)
[2023-05-11 11:42:03.732] [debug]: [W5]: finished processing task 16 (chr17:69100000-72900000 cell #0, 0xc35a95ef10d1fcb0bb6c84546f3b5c771ef2d412b3e39efe52a284b74dc37fbc): collected 285000 interactions throughout 6679 epochs (529 burnin epochs)
[2023-05-11 11:42:03.732] [debug]: [W9]: finished processing task 2 (chr17:10800000-16100000 cell #0, 0x855f4ecc138b144575074359f6ee20516a35db461ce6a6986508bbafe169de0d): collected 397500 interactions throughout 6430 epochs (399 burnin epochs)
[2023-05-11 11:42:03.732] [info]: [io] writing contacts for chr17:10800000-16100000 to file "/tmp/modle-ci-SpASy09Lp3/out.cool"...
[2023-05-11 11:42:03.733] [debug]: [W12]: finished processing task 4 (chr17:33500000-39800000 cell #0, 0xacf52581c195c57e39e98ddbe38add37f6549f1ee8b840846975a121751c865b): collected 472500 interactions throughout 6337 epochs (399 burnin epochs)
[2023-05-11 11:42:03.744] [info]: [io]: written 795000 contacts for chr17:10800000-16100000 to file "/tmp/modle-ci-SpASy09Lp3/out.cool" in 11.809805ms (0.01M nnz out of 0.04M pixels).
[2023-05-11 11:42:03.744] [info]: [io]: writing 1D LEF occupancy profile for chr17:10800000-16100000 to file "/tmp/modle-ci-SpASy09Lp3/out_lef_1d_occupancy.bw"...
[2023-05-11 11:42:03.744] [info]: [io]: writing 1D LEF occupancy profile for chr17:10800000-16100000 to file "/tmp/modle-ci-SpASy09Lp3/out_lef_1d_occupancy.bw" took 203.41us
[2023-05-11 11:42:03.744] [debug]: deallocating a 150x265 contact matrix...
[2023-05-11 11:42:03.744] [debug]: deallocating a vector of size 265...
[2023-05-11 11:42:03.756] [debug]: [W6]: finished processing task 5 (chr17:33500000-39800000 cell #1, 0x635ba9dffc765bb40bf4cebcb302a6b3436beaf044a7fc9200f03ee5f9c538f0): collected 472500 interactions throughout 6145 epochs (205 burnin epochs)
[2023-05-11 11:42:03.756] [info]: [io] writing contacts for chr17:33500000-39800000 to file "/tmp/modle-ci-SpASy09Lp3/out.cool"...
[2023-05-11 11:42:03.757] [debug]: [W2]: finished processing task 10 (chr17:52100000-59500000 cell #0, 0xc8352a6614980a95ad7d7142970c9e77216ea23c47e035be5eb926e3f770e34b): collected 555000 interactions throughout 6320 epochs (353 burnin epochs)
[2023-05-11 11:42:03.768] [info]: [io]: written 945000 contacts for chr17:33500000-39800000 to file "/tmp/modle-ci-SpASy09Lp3/out.cool" in 11.162817ms (0.01M nnz out of 0.05M pixels).
[2023-05-11 11:42:03.768] [info]: [io]: writing 1D LEF occupancy profile for chr17:33500000-39800000 to file "/tmp/modle-ci-SpASy09Lp3/out_lef_1d_occupancy.bw"...
[2023-05-11 11:42:03.768] [info]: [io]: writing 1D LEF occupancy profile for chr17:33500000-39800000 to file "/tmp/modle-ci-SpASy09Lp3/out_lef_1d_occupancy.bw" took 160.875us
[2023-05-11 11:42:03.768] [debug]: deallocating a 150x315 contact matrix...
[2023-05-11 11:42:03.768] [debug]: deallocating a vector of size 315...
[2023-05-11 11:42:03.768] [info]: [io] writing contacts for chr17:40200000-42800000 to file "/tmp/modle-ci-SpASy09Lp3/out.cool"...
[2023-05-11 11:42:03.777] [info]: [io]: written 390000 contacts for chr17:40200000-42800000 to file "/tmp/modle-ci-SpASy09Lp3/out.cool" in 8.85221ms (0.00M nnz out of 0.02M pixels).
[2023-05-11 11:42:03.777] [info]: [io]: writing 1D LEF occupancy profile for chr17:40200000-42800000 to file "/tmp/modle-ci-SpASy09Lp3/out_lef_1d_occupancy.bw"...
[2023-05-11 11:42:03.777] [info]: [io]: writing 1D LEF occupancy profile for chr17:40200000-42800000 to file "/tmp/modle-ci-SpASy09Lp3/out_lef_1d_occupancy.bw" took 172.761us
[2023-05-11 11:42:03.777] [debug]: deallocating a 150x130 contact matrix...
[2023-05-11 11:42:03.777] [debug]: deallocating a vector of size 130...
[2023-05-11 11:42:03.777] [info]: [io] writing contacts for chr17:46800000-49300000 to file "/tmp/modle-ci-SpASy09Lp3/out.cool"...
[2023-05-11 11:42:03.783] [debug]: [W14]: all tasks have been processed: returning!
[2023-05-11 11:42:03.786] [info]: [io]: written 375000 contacts for chr17:46800000-49300000 to file "/tmp/modle-ci-SpASy09Lp3/out.cool" in 9.009251ms (0.00M nnz out of 0.02M pixels).
[2023-05-11 11:42:03.786] [info]: [io]: writing 1D LEF occupancy profile for chr17:46800000-49300000 to file "/tmp/modle-ci-SpASy09Lp3/out_lef_1d_occupancy.bw"...
[2023-05-11 11:42:03.786] [info]: [io]: writing 1D LEF occupancy profile for chr17:46800000-49300000 to file "/tmp/modle-ci-SpASy09Lp3/out_lef_1d_occupancy.bw" took 81.686us
[2023-05-11 11:42:03.786] [debug]: deallocating a 150x125 contact matrix...
[2023-05-11 11:42:03.786] [debug]: deallocating a vector of size 125...
[2023-05-11 11:42:03.797] [debug]: [W7]: finished processing task 11 (chr17:52100000-59500000 cell #1, 0x7dba39cdd272b1bf83640f92d867800050ed7654ab1cd71def39ddaa685862a4): collected 555000 interactions throughout 6293 epochs (323 burnin epochs)
[2023-05-11 11:42:03.797] [info]: [io] writing contacts for chr17:52100000-59500000 to file "/tmp/modle-ci-SpASy09Lp3/out.cool"...
[2023-05-11 11:42:03.809] [info]: [io]: written 1110000 contacts for chr17:52100000-59500000 to file "/tmp/modle-ci-SpASy09Lp3/out.cool" in 12.295641ms (0.01M nnz out of 0.06M pixels).
[2023-05-11 11:42:03.809] [info]: [io]: writing 1D LEF occupancy profile for chr17:52100000-59500000 to file "/tmp/modle-ci-SpASy09Lp3/out_lef_1d_occupancy.bw"...
[2023-05-11 11:42:03.810] [info]: [io]: writing 1D LEF occupancy profile for chr17:52100000-59500000 to file "/tmp/modle-ci-SpASy09Lp3/out_lef_1d_occupancy.bw" took 111.221us
[2023-05-11 11:42:03.810] [debug]: deallocating a 150x370 contact matrix...
[2023-05-11 11:42:03.810] [debug]: deallocating a vector of size 370...
[2023-05-11 11:42:03.810] [info]: [io] writing contacts for chr17:60200000-63100000 to file "/tmp/modle-ci-SpASy09Lp3/out.cool"...
[2023-05-11 11:42:03.817] [debug]: [W8]: all tasks have been processed: returning!
[2023-05-11 11:42:03.819] [info]: [io]: written 435000 contacts for chr17:60200000-63100000 to file "/tmp/modle-ci-SpASy09Lp3/out.cool" in 9.643798ms (0.01M nnz out of 0.02M pixels).
[2023-05-11 11:42:03.819] [info]: [io]: writing 1D LEF occupancy profile for chr17:60200000-63100000 to file "/tmp/modle-ci-SpASy09Lp3/out_lef_1d_occupancy.bw"...
[2023-05-11 11:42:03.819] [info]: [io]: writing 1D LEF occupancy profile for chr17:60200000-63100000 to file "/tmp/modle-ci-SpASy09Lp3/out_lef_1d_occupancy.bw" took 133.274us
[2023-05-11 11:42:03.819] [debug]: deallocating a 150x145 contact matrix...
[2023-05-11 11:42:03.819] [debug]: deallocating a vector of size 145...
[2023-05-11 11:42:03.819] [info]: [io] writing contacts for chr17:64600000-66200000 to file "/tmp/modle-ci-SpASy09Lp3/out.cool"...
[2023-05-11 11:42:03.827] [info]: [io]: written 240000 contacts for chr17:64600000-66200000 to file "/tmp/modle-ci-SpASy09Lp3/out.cool" in 7.99357ms (0.00M nnz out of 0.01M pixels).
[2023-05-11 11:42:03.827] [info]: [io]: writing 1D LEF occupancy profile for chr17:64600000-66200000 to file "/tmp/modle-ci-SpASy09Lp3/out_lef_1d_occupancy.bw"...
[2023-05-11 11:42:03.828] [info]: [io]: writing 1D LEF occupancy profile for chr17:64600000-66200000 to file "/tmp/modle-ci-SpASy09Lp3/out_lef_1d_occupancy.bw" took 82.026us
[2023-05-11 11:42:03.828] [debug]: deallocating a 150x80 contact matrix...
[2023-05-11 11:42:03.828] [debug]: deallocating a vector of size 80...
[2023-05-11 11:42:03.828] [info]: [io] writing contacts for chr17:69100000-72900000 to file "/tmp/modle-ci-SpASy09Lp3/out.cool"...
[2023-05-11 11:42:03.832] [debug]: [W5]: all tasks have been processed: returning!
[2023-05-11 11:42:03.832] [debug]: [W9]: all tasks have been processed: returning!
[2023-05-11 11:42:03.833] [debug]: [W12]: all tasks have been processed: returning!
[2023-05-11 11:42:03.838] [info]: [io]: written 570000 contacts for chr17:69100000-72900000 to file "/tmp/modle-ci-SpASy09Lp3/out.cool" in 10.105671ms (0.01M nnz out of 0.03M pixels).
[2023-05-11 11:42:03.838] [info]: [io]: writing 1D LEF occupancy profile for chr17:69100000-72900000 to file "/tmp/modle-ci-SpASy09Lp3/out_lef_1d_occupancy.bw"...
[2023-05-11 11:42:03.838] [info]: [io]: writing 1D LEF occupancy profile for chr17:69100000-72900000 to file "/tmp/modle-ci-SpASy09Lp3/out_lef_1d_occupancy.bw" took 67.837us
[2023-05-11 11:42:03.838] [debug]: deallocating a 150x190 contact matrix...
[2023-05-11 11:42:03.838] [debug]: deallocating a vector of size 190...
[2023-05-11 11:42:03.838] [info]: [io] writing contacts for chr17:76800000-77200000 to file "/tmp/modle-ci-SpASy09Lp3/out.cool"...
[2023-05-11 11:42:03.846] [info]: [io]: written 60000 contacts for chr17:76800000-77200000 to file "/tmp/modle-ci-SpASy09Lp3/out.cool" in 7.887678ms (0.00M nnz out of 0.00M pixels).
[2023-05-11 11:42:03.846] [info]: [io]: writing 1D LEF occupancy profile for chr17:76800000-77200000 to file "/tmp/modle-ci-SpASy09Lp3/out_lef_1d_occupancy.bw"...
[2023-05-11 11:42:03.846] [info]: [io]: writing 1D LEF occupancy profile for chr17:76800000-77200000 to file "/tmp/modle-ci-SpASy09Lp3/out_lef_1d_occupancy.bw" took 98.888us
[2023-05-11 11:42:03.846] [debug]: deallocating a 150x20 contact matrix...
[2023-05-11 11:42:03.846] [debug]: deallocating a vector of size 20...
[2023-05-11 11:42:03.857] [debug]: [W6]: all tasks have been processed: returning!
[2023-05-11 11:42:03.857] [debug]: [W2]: all tasks have been processed: returning!
[2023-05-11 11:42:03.897] [debug]: [W7]: all tasks have been processed: returning!
[2023-05-11 11:42:06.020] [debug]: [W3]: finished processing task 21 (chr19:0-58617616 cell #1, 0xbf6bff8c51e1ce920d4da20970373f04dc2177a111e26e0be5454e2985187a6f): collected 4396500 interactions throughout 6033 epochs (154 burnin epochs)
[2023-05-11 11:42:06.033] [debug]: [W13]: finished processing task 20 (chr19:0-58617616 cell #0, 0x9b93afbc98c209686dc41783d8e0039a0d048f6edd9719becab4bed89163951a): collected 4396500 interactions throughout 6153 epochs (274 burnin epochs)
[2023-05-11 11:42:06.033] [info]: [io] writing contacts for chr19:0-58617616 to file "/tmp/modle-ci-SpASy09Lp3/out.cool"...
[2023-05-11 11:42:06.066] [info]: [io]: written 8793000 contacts for chr19:0-58617616 to file "/tmp/modle-ci-SpASy09Lp3/out.cool" in 33.047737ms (0.10M nnz out of 0.44M pixels).
[2023-05-11 11:42:06.066] [info]: [io]: writing 1D LEF occupancy profile for chr19:0-58617616 to file "/tmp/modle-ci-SpASy09Lp3/out_lef_1d_occupancy.bw"...
[2023-05-11 11:42:06.066] [info]: [io]: writing 1D LEF occupancy profile for chr19:0-58617616 to file "/tmp/modle-ci-SpASy09Lp3/out_lef_1d_occupancy.bw" took 126.597us
[2023-05-11 11:42:06.066] [debug]: deallocating a 150x2931 contact matrix...
[2023-05-11 11:42:06.066] [debug]: deallocating a vector of size 2931...
[2023-05-11 11:42:06.121] [debug]: [W3]: all tasks have been processed: returning!
[2023-05-11 11:42:06.133] [debug]: [W13]: all tasks have been processed: returning!
[2023-05-11 11:42:06.133] [debug]: waiting for io threads to return...
[2023-05-11 11:42:06.133] [debug]: all background threads returned! Checking if any exception have been raised...
[2023-05-11 11:42:06.133] [debug]: context manager shutdown was successful.
[2023-05-11 11:42:06.133] [info]: simulation terminated without errors in 2.817521006s!

Bye.
Comparing /tmp/modle-ci-SpASy09Lp3/out.cool with modle/test/data/integration_tests/reference_001.cool...
Files are identical
Comparing modle/test/data/integration_tests/reference_001.bw with /tmp/modle-ci-SpASy09Lp3/out_lef_1d_occupancy.bw...
Files are identical

### PASS ###


[2023-05-11 11:44:11.676] [info]: Computing metric(s) for the following 24 intervals:
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
[2023-05-11 11:44:11.683] [warning]: Read 0 contacts for chr1:248956422. SKIPPING!
[2023-05-11 11:44:11.685] [warning]: Read 0 contacts for chr2:242193529. SKIPPING!
[2023-05-11 11:44:11.687] [warning]: Read 0 contacts for chr3:198295559. SKIPPING!
[2023-05-11 11:44:11.688] [warning]: Read 0 contacts for chr4:190214555. SKIPPING!
[2023-05-11 11:44:11.690] [warning]: Read 0 contacts for chr5:181538259. SKIPPING!
[2023-05-11 11:44:11.691] [warning]: Read 0 contacts for chr6:170805979. SKIPPING!
[2023-05-11 11:44:11.693] [warning]: Read 0 contacts for chr7:159345973. SKIPPING!
[2023-05-11 11:44:11.694] [warning]: Read 0 contacts for chr8:145138636. SKIPPING!
[2023-05-11 11:44:11.696] [warning]: Read 0 contacts for chr9:138394717. SKIPPING!
[2023-05-11 11:44:11.697] [warning]: Read 0 contacts for chr10:133797422. SKIPPING!
[2023-05-11 11:44:11.699] [warning]: Read 0 contacts for chr11:135086622. SKIPPING!
[2023-05-11 11:44:11.701] [warning]: Read 0 contacts for chr12:133275309. SKIPPING!
[2023-05-11 11:44:11.702] [warning]: Read 0 contacts for chr13:114364328. SKIPPING!
[2023-05-11 11:44:11.704] [warning]: Read 0 contacts for chr14:107043718. SKIPPING!
[2023-05-11 11:44:11.705] [warning]: Read 0 contacts for chr15:101991189. SKIPPING!
[2023-05-11 11:44:11.707] [warning]: Read 0 contacts for chr16:90338345. SKIPPING!
[2023-05-11 11:44:11.708] [warning]: Read 0 contacts for chr17:83257441. SKIPPING!
[2023-05-11 11:44:11.710] [warning]: Read 0 contacts for chr18:80373285. SKIPPING!
[2023-05-11 11:44:11.711] [warning]: Read 0 contacts for chr19:58617616. SKIPPING!
[2023-05-11 11:44:11.711] [info]: Reading contacts for chr20:64444167...
[2023-05-11 11:44:11.721] [info]: Read 56772 contacts for chr20:64444167 in 9.592302ms
[2023-05-11 11:44:11.722] [info]: Custom metric for vertical stripes from interval chr20:0-64444167 computed in 614.392us.
[2023-05-11 11:44:11.723] [info]: Custom metric for horizontal stripes from interval chr20:0-64444167 computed in 1.970983ms.
[2023-05-11 11:44:11.738] [info]: 2578 values have been written to files "out_custom_score_custom_metric_vertical.{tsv.gz,bw}" in 16.538823ms.
[2023-05-11 11:44:11.740] [info]: 2578 values have been written to files "out_custom_score_custom_metric_horizontal.{tsv.gz,bw}" in 17.294496ms.
[2023-05-11 11:44:11.741] [warning]: Read 0 contacts for chr21:46709983. SKIPPING!
[2023-05-11 11:44:11.741] [warning]: Read 0 contacts for chr22:50818468. SKIPPING!
[2023-05-11 11:44:11.741] [warning]: Read 0 contacts for chrX:156040895. SKIPPING!
[2023-05-11 11:44:11.741] [warning]: Read 0 contacts for chrY:57227415. SKIPPING!
[2023-05-11 11:44:11.741] [info]: DONE in 101.752652ms!
Comparing modle/test/data/integration_tests/4DNFI9GMP2J8_vs_4DNFIFJH2524_mt_eval_custom_score_custom_metric_horizontal.bw with /tmp/modle-tools-E2x3Ia8BMm/out_custom_score_custom_metric_horizontal.bw...
Files are identical
Comparing modle/test/data/integration_tests/4DNFI9GMP2J8_vs_4DNFIFJH2524_mt_eval_custom_score_custom_metric_vertical.bw with /tmp/modle-tools-E2x3Ia8BMm/out_custom_score_custom_metric_vertical.bw...
Files are identical
Comparing modle/test/data/integration_tests/4DNFI9GMP2J8_vs_4DNFIFJH2524_mt_eval_custom_score_custom_metric_horizontal.tsv.gz with /tmp/modle-tools-E2x3Ia8BMm/out_custom_score_custom_metric_horizontal.tsv.gz...
Files are identical
Comparing modle/test/data/integration_tests/4DNFI9GMP2J8_vs_4DNFIFJH2524_mt_eval_custom_score_custom_metric_vertical.tsv.gz with /tmp/modle-tools-E2x3Ia8BMm/out_custom_score_custom_metric_vertical.tsv.gz...
Files are identical

### PASS ###


[2023-05-11 11:44:34.756] [info]: transforming contacts from Cooler at URI "modle/test/data/integration_tests/4DNFI9GMP2J8_chr20_25kbp.cool"...
[2023-05-11 11:44:34.763] [warning]: read 0 contacts for chr1. SKIPPING!
[2023-05-11 11:44:34.765] [warning]: read 0 contacts for chr2. SKIPPING!
[2023-05-11 11:44:34.766] [warning]: read 0 contacts for chr3. SKIPPING!
[2023-05-11 11:44:34.768] [warning]: read 0 contacts for chr4. SKIPPING!
[2023-05-11 11:44:34.769] [warning]: read 0 contacts for chr5. SKIPPING!
[2023-05-11 11:44:34.771] [warning]: read 0 contacts for chr6. SKIPPING!
[2023-05-11 11:44:34.772] [warning]: read 0 contacts for chr7. SKIPPING!
[2023-05-11 11:44:34.774] [warning]: read 0 contacts for chr8. SKIPPING!
[2023-05-11 11:44:34.775] [warning]: read 0 contacts for chr9. SKIPPING!
[2023-05-11 11:44:34.777] [warning]: read 0 contacts for chr10. SKIPPING!
[2023-05-11 11:44:34.778] [warning]: read 0 contacts for chr11. SKIPPING!
[2023-05-11 11:44:34.780] [warning]: read 0 contacts for chr12. SKIPPING!
[2023-05-11 11:44:34.781] [warning]: read 0 contacts for chr13. SKIPPING!
[2023-05-11 11:44:34.783] [warning]: read 0 contacts for chr14. SKIPPING!
[2023-05-11 11:44:34.785] [warning]: read 0 contacts for chr15. SKIPPING!
[2023-05-11 11:44:34.786] [warning]: read 0 contacts for chr16. SKIPPING!
[2023-05-11 11:44:34.788] [warning]: read 0 contacts for chr17. SKIPPING!
[2023-05-11 11:44:34.789] [warning]: read 0 contacts for chr18. SKIPPING!
[2023-05-11 11:44:34.791] [warning]: read 0 contacts for chr19. SKIPPING!
[2023-05-11 11:44:34.791] [info]: processing contacts for chr20...
[2023-05-11 11:44:35.013] [info]: applying Gaussian blur with sigma=1 to contacts for chr20...
[2023-05-11 11:44:35.029] [info]: chr20 processing took 237.595086ms
[2023-05-11 11:44:35.159] [warning]: read 0 contacts for chr21. SKIPPING!
[2023-05-11 11:44:35.159] [warning]: read 0 contacts for chr22. SKIPPING!
[2023-05-11 11:44:35.159] [warning]: read 0 contacts for chrX. SKIPPING!
[2023-05-11 11:44:35.159] [warning]: read 0 contacts for chrY. SKIPPING!
[2023-05-11 11:44:35.159] [info]: DONE! Processed 24 chromosomes in 402.780193ms!
[2023-05-11 11:44:35.159] [info]: Transformed contacts have been saved to file "/tmp/modle-tools-vew8zaAIaC/out_blurred.cool"
[2023-05-11 11:44:35.283] [info]: transforming contacts from Cooler at URI "modle/test/data/integration_tests/4DNFI9GMP2J8_chr20_25kbp.cool"...
[2023-05-11 11:44:35.289] [warning]: read 0 contacts for chr1. SKIPPING!
[2023-05-11 11:44:35.290] [warning]: read 0 contacts for chr2. SKIPPING!
[2023-05-11 11:44:35.292] [warning]: read 0 contacts for chr3. SKIPPING!
[2023-05-11 11:44:35.294] [warning]: read 0 contacts for chr4. SKIPPING!
[2023-05-11 11:44:35.295] [warning]: read 0 contacts for chr5. SKIPPING!
[2023-05-11 11:44:35.297] [warning]: read 0 contacts for chr6. SKIPPING!
[2023-05-11 11:44:35.298] [warning]: read 0 contacts for chr7. SKIPPING!
[2023-05-11 11:44:35.299] [warning]: read 0 contacts for chr8. SKIPPING!
[2023-05-11 11:44:35.301] [warning]: read 0 contacts for chr9. SKIPPING!
[2023-05-11 11:44:35.303] [warning]: read 0 contacts for chr10. SKIPPING!
[2023-05-11 11:44:35.304] [warning]: read 0 contacts for chr11. SKIPPING!
[2023-05-11 11:44:35.306] [warning]: read 0 contacts for chr12. SKIPPING!
[2023-05-11 11:44:35.308] [warning]: read 0 contacts for chr13. SKIPPING!
[2023-05-11 11:44:35.310] [warning]: read 0 contacts for chr14. SKIPPING!
[2023-05-11 11:44:35.311] [warning]: read 0 contacts for chr15. SKIPPING!
[2023-05-11 11:44:35.312] [warning]: read 0 contacts for chr16. SKIPPING!
[2023-05-11 11:44:35.314] [warning]: read 0 contacts for chr17. SKIPPING!
[2023-05-11 11:44:35.315] [warning]: read 0 contacts for chr18. SKIPPING!
[2023-05-11 11:44:35.317] [warning]: read 0 contacts for chr19. SKIPPING!
[2023-05-11 11:44:35.317] [info]: processing contacts for chr20...
[2023-05-11 11:44:35.540] [info]: computing the difference of Gaussians for chr20 (sigma1=1; sigma2=1.6)...
[2023-05-11 11:44:35.582] [info]: chr20 processing took 265.136791ms
[2023-05-11 11:44:35.687] [warning]: read 0 contacts for chr21. SKIPPING!
[2023-05-11 11:44:35.687] [warning]: read 0 contacts for chr22. SKIPPING!
[2023-05-11 11:44:35.687] [warning]: read 0 contacts for chrX. SKIPPING!
[2023-05-11 11:44:35.687] [warning]: read 0 contacts for chrY. SKIPPING!
[2023-05-11 11:44:35.687] [info]: DONE! Processed 24 chromosomes in 404.878477ms!
[2023-05-11 11:44:35.687] [info]: Transformed contacts have been saved to file "/tmp/modle-tools-vew8zaAIaC/out_dog.cool"
Comparing /tmp/modle-tools-vew8zaAIaC/out_blurred.cool with modle/test/data/integration_tests/4DNFI9GMP2J8_chr20_25kbp_blurred.cool...
Files are identical
Comparing /tmp/modle-tools-vew8zaAIaC/out_dog.cool with modle/test/data/integration_tests/4DNFI9GMP2J8_chr20_25kbp_dog.cool...
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
- ___score___ is required and should be set to a number between 0 and 1 representing the average occupancy of the extrusion barrier site defined by the current line. Occupancies between 0.7-0.9 are usually a reasonable starting point.
  Refer to [MoDLE's paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02815-7) for more details regarding the extrusion barrier occupancy.
- ___strand___ is required and is used to define the extrusion barrier direction (should be one of `-`, `+` or `.`).
  As of MoDLE v1.0.1, extrusion barriers are modeled after CTCF barriers.
  Thus, the strand field should be populated with the direction of the CTCF binding site that is being defined.
  Barriers without strand information (i.e. with strand `.`) are ignored.

#### Running a simulation with default settings

```console
user@dev:/tmp/modle$ modle simulate \
    --chrom-sizes examples/data/hg38.chrom.sizes \
    --extrusion-barrier-file examples/data/hg38_extrusion_barriers.bed.xz \
    --output-prefix examples/out/hg38_default

[2023-05-11 11:47:15.718] [info]: running MoDLE v1.0.1
[2023-05-11 11:47:15.718] [info]: complete log will be written to file "examples/out/hg38_default.log"
[2023-05-11 11:47:15.718] [info]: writing simulation parameters to config file "examples/out/hg38_default_config.toml"
[2023-05-11 11:47:15.718] [info]: command: modle simulate --chrom-sizes examples/data/hg38.chrom.sizes --extrusion-barrier-file examples/data/hg38_extrusion_barriers.bed.xz --output-prefix examples/out/hg38_default
[2023-05-11 11:47:15.719] [info]: simulation will use up to 16 out of 16 available CPU cores.
[2023-05-11 11:47:15.719] [info]: using --target-contact-density=1.00 as stopping criterion.
[2023-05-11 11:47:15.719] [info]: contact sampling strategy: tad-plus-loop-with-noise.
[2023-05-11 11:47:15.719] [info]: contact matrix resolution: 5000bp
[2023-05-11 11:47:15.719] [info]: importing chromosomes from "examples/data/hg38.chrom.sizes"...
[2023-05-11 11:47:15.719] [info]: imported 24 chromosomes in 239.985us.
[2023-05-11 11:47:15.719] [info]: importing extrusion barriers from "examples/data/hg38_extrusion_barriers.bed.xz"...
[2023-05-11 11:47:15.761] [info]: imported 38815 barriers from "examples/data/hg38_extrusion_barriers.bed.xz" in 41.776052ms.
[2023-05-11 11:47:15.762] [info]: spawning worker thread W0...
[2023-05-11 11:47:15.762] [info]: spawning worker thread W1...
[2023-05-11 11:47:15.762] [info]: spawning worker thread W3...
...
[2023-05-11 11:47:15.764] [info]: begin processing chr1:0-248956422: simulating ~37485 epochs across 512 cells using 4979 LEFs and 3518 barriers (~73 epochs per cell)...
[2023-05-11 11:47:37.627] [info]: begin processing chr2:0-242193529: simulating ~37501 epochs across 512 cells using 4844 LEFs and 2974 barriers (~73 epochs per cell)...
[2023-05-11 11:47:43.592] [info]: [io] writing contacts for chr1:0-248956422 to file "examples/out/hg38_default.cool"...
[2023-05-11 11:47:46.182] [info]: [io]: written 29875200 contacts for chr1:0-248956422 to file "examples/out/hg38_default.cool" in 2.590415163s (3.47M nnz out of 29.88M pixels).
[2023-05-11 11:47:46.186] [info]: [io]: writing 1D LEF occupancy profile for chr1:0-248956422 to file "examples/out/hg38_default_lef_1d_occupancy.bw"...
[2023-05-11 11:47:46.196] [info]: [io]: writing 1D LEF occupancy profile for chr1:0-248956422 to file "examples/out/hg38_default_lef_1d_occupancy.bw" took 9.964882ms
[2023-05-11 11:47:57.200] [info]: begin processing chr3:0-198295559: simulating ~37474 epochs across 512 cells using 3966 LEFs and 2469 barriers (~73 epochs per cell)...
...
[2023-05-11 11:51:49.400] [info]: [io] writing contacts for chrY:0-57227415 to file "examples/out/hg38_default.cool"...
[2023-05-11 11:51:49.744] [info]: [io]: written 6867600 contacts for chrY:0-57227415 to file "examples/out/hg38_default.cool" in 344.473166ms (0.92M nnz out of 6.87M pixels).
[2023-05-11 11:51:49.744] [info]: [io]: writing 1D LEF occupancy profile for chrY:0-57227415 to file "examples/out/hg38_default_lef_1d_occupancy.bw"...
[2023-05-11 11:51:49.746] [info]: [io]: writing 1D LEF occupancy profile for chrY:0-57227415 to file "examples/out/hg38_default_lef_1d_occupancy.bw" took 1.565303ms
[2023-05-11 11:51:51.925] [info]: simulation terminated without errors in 4m36.206746164s!

Bye.
```

The above command will create folder `examples/out/` (if it doesn't already exist), and write the following files inside it:

```
examples/out
├── hg38_default_config.toml
├── hg38_default.cool
├── hg38_default_lef_1d_occupancy.bw
└── hg38_default.log
```

Output files produced by the above command are available [here](https://doi.org/10.5281/zenodo.7924253).

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

Output files produced by the above command are available [here](https://doi.org/10.5281/zenodo.7924253).

![Simulation w/ high contact density](examples/images/higlass_003.avif)

MoDLE's config files are in [TOML format](https://toml.io/en/).

Adding a line like `my-option=my_value` to the config file has the same effect as passing `--my-option=my_value` on the CLI.

For an up-to-date list of supported CLI options, please refer to MoDLE's help message:

```console
user@dev:/tmp/modle$ modle simulate --help

Simulate loop extrusion and write resulting molecular contacts in a .cool file.
Usage: modle simulate [OPTIONS]

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
    -g,--genomic-intervals      Path to BED3+ file with the genomic regions to be simulated.
                                Intervals listed in this file should be a subset of the chromosomes defined in the .chrom.sizes files.
                                Intervals referring to chromosomes not listed in the .chrom.sizes file are ignored.
    -f,--force                  Overwrite existing files (if any).
    -o,--output-prefix REQUIRED Output prefix.
                                Can be an absolute or relative path including the file name but without the extension.
                                Example: running modle sim -o /tmp/my_simulation ... yields the following files:
                                         - /tmp/my_simulation.cool
                                         - /tmp/my_simulation.log
                                         - /tmp/my_simulation_config.toml
                                         - /tmp/my_simulation_lef_1d_occupancy.bw
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
