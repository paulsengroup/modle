# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

from conan import ConanFile
from conan.tools.build import check_min_cppstd
from conan.tools.scm import Version
from conan.errors import ConanInvalidConfiguration

required_conan_version = ">=1.51.3"


class MoDLE(ConanFile):
  name = "MoDLE"
  version = "1.0.0-rc.7"
  homepage = "https://github.com/paulsengroup/modle"
  license = "MIT"
  author = "Roberto Rossini (roberros@uio.no)"
  settings = "os", "compiler", "build_type", "arch"
  requires = ["abseil/20220623.0@#732381dc99db29b4cfd293684891da56",
              "bitflags/1.5.0@#f49a584785e75c40bf6970615c452b4e",
              "boost/1.80.0@#3264cfe783d4202f76c96b7d5399ff17",
              "bshoshany-thread-pool/3.3.0@#22e99aee6babc19e679754d95dad2de4",
              "bzip2/1.0.8@#7cbb1b682399e077c6b9b2d1bc52da81",
              "catch2/3.1.0@#13edd92657a2a23d108e6180d7f1026a",
              "cli11/2.2.0@#33cd38722fa134b15ae308dfb4e6c942",
              "concurrentqueue/1.0.3@#c1cb7d960d8b64073643b45fa63f0bd1",
              "cpp-sort/1.13.1@#32bc4841eedc6f06551d2317f889f920",
              "fast_float/3.5.1@#63ccdfa6e4dbc05de4bc598258b6a12f",
              "fmt/9.1.0@#78313935d914fe65715ad39687445de6",
              "hdf5/1.12.2@#b01e96ebe1e351ee1d65ae49a347c29c",
              "libarchive/3.6.1@#9e38393eaf7725e0396a8052dd7e29eb",
              "libbigwig/0.4.7@#c51cf59700963a74fe0089a81293c91c",
              "libcuckoo/0.3.1@#2fadf7043e85444c19c1ededf9e43653",
              "lz4/1.9.4@#8699136f62681eeba3dc033f85cc63e2",
              "lzo/2.10@#491362961bbb7678fdd623555a823b60",
              "range-v3/0.12.0@#ef25fe38953fd3847f5bf0dabe1df122",
              "readerwriterqueue/1.0.6@#a95c8da3d68822dec4d4c13fff4b5c96",
              "spdlog/1.9.2@#269bc16aef6a95887eb7c0f2e6aa4529",
              "tomlplusplus/3.2.0@#86c72c809aa5489e6cf0ce4dc5db0002",
              "xoshiro-cpp/1.1@#be8a2825995d67cf2fb26f01accb4f4d",
              "xxhash/0.8.1@#bbf2867eb2cfe186d9bc857da8e00752",
              "xz_utils/5.2.5@#dd8565c00d13a3a6fb2f01214f73a076",
              "zlib/1.2.13@#647c91ed13c0a6c9ea9add6ca968ea93",
              "zstd/1.5.2@#0183cee53cfff406eac0a68be4c4f418"]

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

    self.options["libarchive"].with_zlib = True
    self.options["libarchive"].with_bzip2 = True
    self.options["libarchive"].with_lz4 = True
    self.options["libarchive"].with_lzo = True
    self.options["libarchive"].with_lzma = True
    self.options["libarchive"].with_zstd = True

    self.options["xxhash"].utility = False

  def imports(self):
    self.copy("license*", dst="licenses", folder=True, ignore_case=True)
