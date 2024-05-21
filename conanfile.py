# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

from conan import ConanFile
from conan.tools.build import check_min_cppstd

required_conan_version = ">=1.53.0"


class MoDLEConan(ConanFile):
    name = "modle"
    description = "High-performance stochastic modeling of DNA loop extrusion interactions."
    license = "MIT"
    topics = ("modle", "bioinformatics")
    homepage = "https://github.com/paulsengroup/modle"
    url = "https://github.com/paulsengroup/modle"
    package_type = "application"
    settings = "os", "arch", "compiler", "build_type"

    options = {
        "shared": [True, False],
        "fPIC": [True, False],
    }

    default_options = {
        "shared": False,
        "fPIC": True,
    }

    generators = "CMakeDeps"

    @property
    def _min_cppstd(self):
        return 17

    def requirements(self):
        self.requires("abseil/20240116.2#932ead913d27bb5087303c2880a7b06d")
        self.requires("bitflags/1.5.0#626da9d1913161321841f30caf9b994e")
        self.requires("boost/1.85.0#7926babdab0a9779cc164d0af6c28d5e")
        self.requires("bshoshany-thread-pool/4.1.0#be1802a8768416a6c9b1393cf0ce5e9c")
        self.requires("bzip2/1.0.8#457c272f7da34cb9c67456dd217d36c4")
        self.requires("catch2/3.6.0#819bc5a82c2cb626916fc18ee1dbc45f")
        self.requires("cli11/2.4.2#1b431bda2fb2cd3efed633899abcd8cc")
        self.requires("concurrentqueue/1.0.4#1e48e1c712bcfd892087c9c622a51502")
        self.requires("cpp-sort/1.15.0#a48647a61f03d08b77d15fb5a8bfbe5e")
        self.requires("fast_float/6.1.1#e29acaa3d0543dee343abe3f6815346e")
        self.requires("fmt/10.2.1#9199a7a0611866dea5c8849a77467b25")
        self.requires("hictk/0.0.12#8e413cd45528da38b5a41ccffee41d6d")
        self.requires("libarchive/3.7.4#db39e5a5cddfd3a1884daaab4b7942fb")
        self.requires("libbigwig/0.4.7#3f34dc76212124688b984b781e6f853d")
        self.requires("libcuckoo/0.3.1#7e514d4c23a9aba3d8d80758824e9dc0")
        self.requires("lz4/1.9.4#1217a03c990b79aa34ed0faede18f534")
        self.requires("lzo/2.10#5725914235423c771cb1c6b607109b45")
        self.requires("range-v3/0.12.0#4c05d91d7b40e6b91b44b5345ac64408")
        self.requires("readerwriterqueue/1.0.6#aaa5ff6fac60c2aee591e9e51b063b83")
        self.requires("spdlog/1.14.1#972bbf70be1da4bc57ea589af0efde03", force=True)
        self.requires("tomlplusplus/3.4.0#92c93d4de6b0a6a2d0ae9b4430e09c9b")
        self.requires("xoshiro-cpp/1.1#20f566efb3e2bf6e1b813d4abfc5e62c")
        self.requires("xxhash/0.8.2#03fd1c9a839b3f9cdf5ea9742c312187")
        self.requires("xz_utils/5.4.5#51e5a6e6564f4ea3afd79def01f035ad", force=True)
        self.requires("zlib/1.3.1#f52e03ae3d251dec704634230cd806a2")
        self.requires("zstd/1.5.6#afefe79a309bc2a7b9f56c2093504c8b", force=True)

    def validate(self):
        if self.settings.get_safe("compiler.cppstd"):
            check_min_cppstd(self, self._min_cppstd)

    def configure(self):
        if self.settings.compiler in ["clang", "gcc"]:
            self.settings.compiler.libcxx = "libstdc++11"

        self.options["boost"].system_no_deprecated = True
        self.options["boost"].asio_no_deprecated = True
        self.options["boost"].filesystem_no_deprecated = True
        self.options["boost"].filesystem_version = 4
        self.options["boost"].zlib = True
        self.options["boost"].bzip2 = True
        self.options["boost"].lzma = True
        self.options["boost"].zstd = True
        self.options["boost"].without_atomic = False
        self.options["boost"].without_charconv = True
        self.options["boost"].without_chrono = True
        self.options["boost"].without_container = True
        self.options["boost"].without_context = True
        self.options["boost"].without_contract = True
        self.options["boost"].without_coroutine = True
        self.options["boost"].without_date_time = True
        self.options["boost"].without_exception = False
        self.options["boost"].without_fiber = True
        self.options["boost"].without_filesystem = False
        self.options["boost"].without_graph = True
        self.options["boost"].without_graph_parallel = True
        self.options["boost"].without_iostreams = False
        self.options["boost"].without_json = True
        self.options["boost"].without_locale = True
        self.options["boost"].without_log = True
        self.options["boost"].without_math = False
        self.options["boost"].without_mpi = True
        self.options["boost"].without_nowide = True
        self.options["boost"].without_program_options = True
        self.options["boost"].without_python = True
        self.options["boost"].without_random = False
        self.options["boost"].without_regex = False
        self.options["boost"].without_serialization = False
        self.options["boost"].without_stacktrace = True
        self.options["boost"].without_system = False
        self.options["boost"].without_test = True
        self.options["boost"].without_thread = True
        self.options["boost"].without_timer = True
        self.options["boost"].without_type_erasure = True
        self.options["boost"].without_url = True
        self.options["boost"].without_wave = True
        self.options["bzip2"].enable_executable = True
        self.options["fmt"].header_only = True
        self.options["hdf5"].enable_cxx = False
        self.options["hdf5"].hl = False
        self.options["hdf5"].threadsafe = False
        self.options["hdf5"].parallel = False
        # self.options["hictk"].with_arrow = False
        self.options["hictk"].with_eigen = False
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
        self.options["libbigwig"].with_curl = False
        self.options["spdlog"].header_only = True
        self.options["xxhash"].utility = False
        self.options["zstd"].build_programs = False
