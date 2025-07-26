# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

from conan import ConanFile
from conan.tools.build import check_min_cppstd

required_conan_version = ">=2.15.0"


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
        self.requires("bitflags/1.5.0#626da9d1913161321841f30caf9b994e")
        self.requires("boost/1.88.0#79153791f7764054ea4debb1c0f667c6")
        self.requires("bshoshany-thread-pool/5.0.0#d94da300363f0c35b8f41b2c5490c94d")
        self.requires("bzip2/1.0.8#00b4a4658791c1f06914e087f0e792f5")
        self.requires("catch2/3.8.1#141f4cd552b86c7278436c434473ae2f")
        self.requires("cli11/2.5.0#1b7c81ea2bff6279eb2150bbe06a200a")
        self.requires("concurrentqueue/1.0.4#1e48e1c712bcfd892087c9c622a51502")
        self.requires("cpp-sort/1.16.0#45e56dfc2a55a3e8c645576d4761dc73")
        self.requires("fast_float/8.0.2#846ad0ebab16bc265c511095c3b490e9")
        self.requires("fmt/11.2.0#579bb2cdf4a7607621beea4eb4651e0f")
        self.requires("hdf5/1.14.6#6f1acd01d23d00735fe97300f4d5980c", force=True)
        self.requires("highfive/2.10.0#c975a16d7fe3655c173f8a9aab16b416")
        self.requires("libarchive/3.8.1#5cf685686322e906cb42706ab7e099a8")
        self.requires("libbigwig/0.4.7#3d7ab50831d3bdc1bba52c01e405899d")
        self.requires("libcuckoo/0.3.1#7e514d4c23a9aba3d8d80758824e9dc0")
        self.requires("libdeflate/1.23#4994bea7cf7e93789da161fac8e26a53")
        self.requires("lz4/1.10.0#59fc63cac7f10fbe8e05c7e62c2f3504", force=True)
        self.requires("lzo/2.10#5725914235423c771cb1c6b607109b45")
        self.requires("nlohmann_json/3.12.0#2d634ab0ec8d9f56353e5ccef6d6612c")
        self.requires("parallel-hashmap/2.0.0#82acae64ffe2693fff5fb3f9df8e1746")
        self.requires("readerwriterqueue/1.0.6#aaa5ff6fac60c2aee591e9e51b063b83")
        self.requires("span-lite/0.11.0#519fd49fff711674cfed8cd17d4ed422")
        self.requires("spdlog/1.15.3#3ca0e9e6b83af4d0151e26541d140c86")
        self.requires("tomlplusplus/3.4.0#85dbfed71376fb8dc23cdcc0570e4727")
        self.requires("xoshiro-cpp/1.1#20f566efb3e2bf6e1b813d4abfc5e62c")
        self.requires("xxhash/0.8.3#681d36a0a6111fc56e5e45ea182c19cc")
        self.requires("xz_utils/5.4.5#b885d1d79c9d30cff3803f7f551dbe66")
        self.requires("zlib/1.3.1#b8bc2603263cf7eccbd6e17e66b0ed76")
        self.requires("zstd/1.5.7#fde461c0d847a22f16d3066774f61b11")

    def validate(self):
        if self.settings.get_safe("compiler.cppstd"):
            check_min_cppstd(self, self._min_cppstd)

    def configure(self):
        self.options["boost"].system_no_deprecated = True
        self.options["boost"].asio_no_deprecated = True
        self.options["boost"].filesystem_no_deprecated = True
        self.options["boost"].filesystem_use_std_fs = True
        self.options["boost"].filesystem_version = 4
        self.options["boost"].zlib = False
        self.options["boost"].bzip2 = False
        self.options["boost"].lzma = False
        self.options["boost"].zstd = False
        self.options["boost"].without_atomic = True
        self.options["boost"].without_charconv = True
        self.options["boost"].without_chrono = True
        self.options["boost"].without_cobalt = True
        # without_container is set to False to workaround https://github.com/conan-io/conan-center-index/issues/26890
        self.options["boost"].without_container = False
        self.options["boost"].without_context = True
        self.options["boost"].without_contract = True
        self.options["boost"].without_coroutine = True
        self.options["boost"].without_date_time = True
        self.options["boost"].without_exception = True
        self.options["boost"].without_fiber = True
        self.options["boost"].without_filesystem = True
        self.options["boost"].without_graph = True
        self.options["boost"].without_graph_parallel = True
        self.options["boost"].without_iostreams = True
        self.options["boost"].without_json = True
        self.options["boost"].without_locale = True
        self.options["boost"].without_log = True
        self.options["boost"].without_math = False
        self.options["boost"].without_mpi = True
        self.options["boost"].without_nowide = True
        self.options["boost"].without_process = True
        self.options["boost"].without_program_options = True
        self.options["boost"].without_python = True
        self.options["boost"].without_random = False
        self.options["boost"].without_regex = False
        self.options["boost"].without_serialization = True
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
