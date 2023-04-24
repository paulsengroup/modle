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
        "abseil/20230125.2@#24e0584aa2ae6d2fb9e9cd0027354ee2",
        "bitflags/1.5.0@#626da9d1913161321841f30caf9b994e",
        "boost/1.81.0@#05966cf97c97cdd334623324750f67a0",
        "bshoshany-thread-pool/3.3.0@#22e99aee6babc19e679754d95dad2de4",
        "bzip2/1.0.8@#411fc05e80d47a89045edc1ee6f23c1d",
        "catch2/3.3.2@#be6a2f0225146ba4fd8573ee9013e5ae",
        "cli11/2.3.2@#8ccdf14fb1ad53532d498c16ae580b4b",
        "concurrentqueue/1.0.3@#c1cb7d960d8b64073643b45fa63f0bd1",
        "cpp-sort/1.14.0@#3453aaaf83c1dae4214ca3b5c4c3a5c8",
        "fast_float/4.0.0@#90cac63d3ae321f6318b9abf8af9cbb1",
        "fmt/9.1.0@#e747928f85b03f48aaf227ff897d9634",
        "hdf5/1.14.0@#1192b687ae347bfcd40d3a7f620b5a57",  # Coolerpp
        "highfive/2.7.1@#a73bc6937c9add30c9d47a7a70a466eb",  # Coolerpp
        "libarchive/3.6.2@#97a1b90f7143a0b5010b3b7f8426c07a",
        "libbigwig/0.4.7@#c51cf59700963a74fe0089a81293c91c",
        "libcuckoo/0.3.1@#2fadf7043e85444c19c1ededf9e43653",
        "lz4/1.9.4@#bce1f314775b83c195dffc8e177ff368",
        "lzo/2.10@#850936b57f92160152bd9385b1770a69",
        "range-v3/0.12.0@#abb9932b80d96eaf8fc992debe3870ed",
        "readerwriterqueue/1.0.6@#aaa5ff6fac60c2aee591e9e51b063b83",
        "spdlog/1.11.0@#d0fdbaa523550b89156084bf42b41c90",
        "tomlplusplus/3.3.0@#ebb2a36577011fb1959b0de8c1509a6d",
        "tsl-hopscotch-map/2.3.0@#497d3f41172cefe2df9ac17692c52734",  # Coolerpp
        "tsl-ordered-map/1.1.0@#c8a6d6831f079d7fb012c46b5bcfa767",  # Coolerpp
        "xoshiro-cpp/1.1@#be8a2825995d67cf2fb26f01accb4f4d",
        "xxhash/0.8.1@#b60fcc5f9821c988a155935d87562e1d",
        "xz_utils/5.4.2@#b6ee8320403def553418874435445982",
        "zlib/1.2.13@#13c96f538b52e1600c40b88994de240f",
        "zstd/1.5.5@#93372fe14bb7883bd4de82914e0a1841",
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
