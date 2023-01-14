# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

from conan import ConanFile
from conan.errors import ConanInvalidConfiguration
from conan.tools.build import check_min_cppstd
from conan.tools.scm import Version

required_conan_version = ">=1.53.0"


class MoDLE(ConanFile):
    name = "MoDLE"
    version = "1.0.0"
    homepage = "https://github.com/paulsengroup/modle"
    license = "MIT"
    author = "Roberto Rossini (roberros@uio.no)"
    settings = "os", "compiler", "build_type", "arch"
    requires = [
        "abseil/20220623.1@#266f2e66c0dcd428efb9b24c2a5c8f05",
        "bitflags/1.5.0@#f49a584785e75c40bf6970615c452b4e",
        "boost/1.81.0@#2a6097db2254b5ea31fcb54377982c4b",
        "bshoshany-thread-pool/3.3.0@#22e99aee6babc19e679754d95dad2de4",
        "bzip2/1.0.8@#464be69744fa6d48ed01928cfe470008",
        "catch2/3.2.1@#f4c25988de20c21a24a1a24dab87e9b5",
        "cli11/2.3.2@#8ccdf14fb1ad53532d498c16ae580b4b",
        "concurrentqueue/1.0.3@#c1cb7d960d8b64073643b45fa63f0bd1",
        "cpp-sort/1.14.0@#3453aaaf83c1dae4214ca3b5c4c3a5c8",
        "fast_float/3.8.1@#ef5291bee6d3ad59f986096df16f43e6",
        "fmt/9.1.0@#811e918ca4b4e0b9ddd6d5a2883efa82",
        "hdf5/1.12.2@#dc802e78ddfbfd6dac7a31cc004c8db0",  # Coolerpp
        "highfive/2.6.2@#829275c7502abf62045167adb79ee581",  # Coolerpp
        "libarchive/3.6.2@#43ce22888cd77d1306e332b41fce8d6c",
        "libbigwig/0.4.7@#c51cf59700963a74fe0089a81293c91c",
        "libcuckoo/0.3.1@#2fadf7043e85444c19c1ededf9e43653",
        "lz4/1.9.4@#8699136f62681eeba3dc033f85cc63e2",
        "lzo/2.10@#491362961bbb7678fdd623555a823b60",
        "range-v3/0.12.0@#abb9932b80d96eaf8fc992debe3870ed",
        "readerwriterqueue/1.0.6@#a95c8da3d68822dec4d4c13fff4b5c96",
        "spdlog/1.11.0@#1e0f4eb6338d05e4bd6fcc6bf4734172",
        "tomlplusplus/3.2.0@#86c72c809aa5489e6cf0ce4dc5db0002",
        "tsl-hopscotch-map/2.3.0@#497d3f41172cefe2df9ac17692c52734",  # Coolerpp
        "tsl-ordered-map/1.0.0@#dceec5bd0243a5c659939d0295dfee50",  # Coolerpp
        "xoshiro-cpp/1.1@#be8a2825995d67cf2fb26f01accb4f4d",
        "xxhash/0.8.1@#b60fcc5f9821c988a155935d87562e1d",
        "xz_utils/5.4.0@#394283edeac4658356894312aa0fdcbb",
        "zlib/1.2.13@#13c96f538b52e1600c40b88994de240f",
        "zstd/1.5.2@#485dd9c1f9245f0f362730e8b8031a17",
    ]

    generators = "cmake", "cmake_find_package", "cmake_find_package_multi"

    @property
    def _minimum_cpp_standard(self):
        return 17

    def validate(self):
        if self.settings.get_safe("compiler.cppstd"):
            check_min_cppstd(self, self._minimum_cpp_standard)

    def configure(self):
        if self.settings.compiler in ["clang", "gcc"]:
            self.settings.compiler.libcxx = "libstdc++11"

        # Set settings for dependencies
        self.options["boost"].system_no_deprecated = True
        self.options["boost"].asio_no_deprecated = True
        self.options["boost"].filesystem_no_deprecated = True
        self.options["boost"].zlib = True
        self.options["boost"].bzip2 = True
        self.options["boost"].lzma = True
        self.options["boost"].zstd = True
        self.options["boost"].without_chrono = True
        self.options["boost"].without_context = True
        self.options["boost"].without_contract = True
        self.options["boost"].without_coroutine = True
        self.options["boost"].without_date_time = True
        self.options["boost"].without_fiber = True
        self.options["boost"].without_graph = True
        self.options["boost"].without_graph_parallel = True
        self.options["boost"].without_json = True
        self.options["boost"].without_locale = True
        self.options["boost"].without_log = True
        self.options["boost"].without_mpi = True
        self.options["boost"].without_nowide = True
        self.options["boost"].without_program_options = True
        self.options["boost"].without_python = True
        self.options["boost"].without_serialization = False
        self.options["boost"].without_test = True
        self.options["boost"].without_thread = True
        self.options["boost"].without_timer = True
        self.options["boost"].without_type_erasure = True
        self.options["boost"].without_wave = True

        self.options["bizp2"].enable_executable = False

        # Coolerpp stuff
        self.options["highfive"].with_boost = False
        self.options["highfive"].with_eigen = False
        self.options["highfive"].with_opencv = False
        self.options["highfive"].with_xtensor = False

        self.options["libarchive"].with_zlib = True
        self.options["libarchive"].with_bzip2 = True
        self.options["libarchive"].with_lz4 = True
        self.options["libarchive"].with_lzo = True
        self.options["libarchive"].with_lzma = True
        self.options["libarchive"].with_zstd = True

        self.options["xxhash"].utility = False

    def imports(self):
        self.copy("license*", dst="licenses", folder=True, ignore_case=True)
