# Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

from conans import ConanFile, tools


class MoDLE(ConanFile):
    # Note: options are copied from CMake boolean options.
    # When turned off, CMake sometimes passes them as empty strings.
    options = {
        "enable_testing": ["ON", "OFF", True, False, ""]
    }

    default_options = {"enable_testing": "ON"}

    name = "MoDLE"
    version = "0.0.1"
    homepage = "https://github.com/robomics"
    license = "MIT"
    author = "Roberto Rossini (roberros@uio.no)"
    settings = "os", "compiler", "build_type", "arch"
    requires = ["abseil/20211102.0",
                "boost/1.78.0",
                "bzip2/1.0.8",
                "cli11/2.1.2",
                "concurrentqueue/1.0.3",
                "cpp-sort/1.12.1",
                "fmt/8.1.1",
                "hdf5/1.12.0",
                "libarchive/3.5.2",
                "lz4/1.9.3",
                "lzo/2.10",
                "readerwriterqueue/1.0.5",
                "spdlog/1.9.2",
                "tomlplusplus/2.5.0",
                "xxhash/0.8.0",
                "xz_utils/5.2.5",
                "zlib/1.2.11"]
    # "zstd/1.4.8"]

    generators = "cmake", "gcc", "txt", "cmake_find_package"

    def requirements(self):
        if bool(self.options.enable_testing):
            self.requires("catch2/2.13.8")

    def configure(self):
        if self.settings.compiler.cppstd:
            tools.check_min_cppstd(self, 17)

        # Set settings for dependencies
        self.options["boost"].system_no_deprecated = True
        self.options["boost"].asio_no_deprecated = True
        self.options["boost"].filesystem_no_deprecated = True
        self.options["boost"].zlib = True
        self.options["boost"].bzip2 = True
        self.options["boost"].lzma = True
        # self.options["boost"].zstd = True
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
        self.options["boost"].without_serialization = True
        self.options["boost"].without_test = True
        self.options["boost"].without_thread = True
        self.options["boost"].without_timer = True
        self.options["boost"].without_type_erasure = True
        self.options["boost"].without_wave = True

        self.options["spdlog"].header_only = True

        self.options["bizp2"].enable_executable = False

        self.options["libarchive"].with_zlib = True
        self.options["libarchive"].with_bzip2 = True
        self.options["libarchive"].with_lz4 = True
        self.options["libarchive"].with_lzo = True
        self.options["libarchive"].with_lzma = True
        # self.options["libarchive"].with_zstd = True

        self.options["xxhash"].utility = False

    def imports(self):
        self.copy("license*", dst="licenses", folder=True, ignore_case=True)
