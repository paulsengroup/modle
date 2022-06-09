<!--
Copyright (C) 2022 Roberto Rossini <roberros@uio.no>

SPDX-License-Identifier: MIT
-->

# README

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6628728.svg)](https://doi.org/10.5281/zenodo.6628728)

Test dataset(s) are archived on Zenodo ([10.5281/zenodo.6628728](https://zenodo.org/record/6628728)) and are downloaded by CMake during project configuration.

Datasets are extracted by CMake in this folder.

If you need to run MoDLE unit tests on a machine without internet access you should first manually download the tar archive from Zenodo and place it in the same folder as this README.

You can get the download link by running the following command from the repository root:

```bash
grep 'zenodo' cmake/FetchTestDataset.cmake
```
