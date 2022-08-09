# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

from conans import ConanFile, tools

required_conan_version = ">=1.45.0"


class MoDLE(ConanFile):
  name = "MoDLE"
  version = "1.0.0-rc.6"
  homepage = "https://github.com/paulsengroup/modle"
  license = "MIT"
  author = "Roberto Rossini (roberros@uio.no)"
  settings = "os", "compiler", "build_type", "arch"
  requires = ["abseil/20211102.0",
              "boost/1.79.0",
              "bshoshany-thread-pool/3.3.0",
              "bzip2/1.0.8",
              "catch2/3.1.0",
              "cli11/2.2.0",
              "concurrentqueue/1.0.3",
              "cpp-sort/1.13.0",
              "fast_float/3.5.1",
              "fmt/9.0.0",
              "hdf5/1.12.2",
              "libarchive/3.6.1",
              "libcuckoo/0.3.1",
              "lz4/1.9.3",
              "lzo/2.10",
              "range-v3/0.12.0",
              "readerwriterqueue/1.0.6",
              "spdlog/1.9.2",
              "tomlplusplus/3.1.0",
              "xoshiro-cpp/1.1",
              "xxhash/0.8.1",
              "xz_utils/5.2.5",
              "zlib/1.2.12",
              "zstd/1.5.2"]

  generators = "cmake", "cmake_find_package", "cmake_find_package_multi"

  def validate(self):
    if self.settings.compiler.get_safe("cppstd"):
      tools.check_min_cppstd(self, 17)

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

    self.options["libarchive"].with_zlib = True
    self.options["libarchive"].with_bzip2 = True
    self.options["libarchive"].with_lz4 = True
    self.options["libarchive"].with_lzo = True
    self.options["libarchive"].with_lzma = True
    self.options["libarchive"].with_zstd = True

    self.options["xxhash"].utility = False

  def imports(self):
    self.copy("license*", dst="licenses", folder=True, ignore_case=True)
