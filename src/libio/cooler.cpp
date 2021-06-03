#include "modle/cooler.hpp"

// IWYU pragma: no_include <H5DaccProp.h>
// IWYU pragma: no_include <H5DataSet.h>
// IWYU pragma: no_include <H5DataSpace.h>
// IWYU pragma: no_include <H5DcreatProp.h>
// IWYU pragma: no_include <H5Exception.h>
// IWYU pragma: no_include <H5File.h>
// IWYU pragma: no_include <H5Group.h>
// IWYU pragma: no_include <H5PredType.h>
// IWYU pragma: no_include <H5Public.h>
// IWYU pragma: no_include "H5Public.h"
// IWYU pragma: no_include "H5Ppublic.h"
// IWYU pragma: no_include <H5SPublic.h>
// IWYU pragma: no_include <H5Spublic.h>
// IWYU pragma: no_include <H5StrType.h>
// IWYU pragma: no_include <ext/alloc_traits.h>
// IWYU pragma: no_include <boost/core/checked_delete.hpp>
// IWYU pragma: no_include <boost/exception/detail/error_info_impl.hpp>
// IWYU pragma: no_include "hdf5_impl.hpp"

#include <H5Cpp.h>                                // IWYU pragma: keep
#include <absl/strings/match.h>                   // for EndsWithIgnoreCase, StartsWith
#include <absl/strings/str_cat.h>                 // for StrCat, StrAppend
#include <absl/strings/str_join.h>                // for StrJoin
#include <absl/strings/strip.h>                   // for StripPrefix, ConsumePrefix, StripS...
#include <absl/time/clock.h>                      // for Now
#include <absl/time/time.h>                       // for FormatTime, UTCTimeZone
#include <absl/types/span.h>                      // for Span, MakeConstSpan
#include <fmt/format.h>                           // for format, FMT_STRING, print
#include <fmt/ostream.h>                          // for formatbuf<>::int_type
#include <readerwriterqueue/readerwriterqueue.h>  // for BlockingReaderWriterQueue

#include <algorithm>    // for max, min, fill, find
#include <array>        // for array, array<>::value_type
#include <cassert>      // for assert
#include <chrono>       // for milliseconds
#include <cmath>        // for isnan, round
#include <cstddef>      // for size_t
#include <cstdint>      // for int64_t, uint32_t, uint8_t, uint_f...
#include <cstdio>       // for stderr
#include <filesystem>   // for operator<<, path
#include <iterator>     // for distance
#include <limits>       // for numeric_limits
#include <memory>       // for unique_ptr, allocator, make_unique
#include <sstream>      // for basic_stringbuf<>::int_type, basic...
#include <stdexcept>    // for runtime_error, logic_error
#include <string>       // for string, basic_string, operator!=
#include <string_view>  // for string_view, basic_string_view
#include <thread>       // IWYU pragma: keep
#include <utility>      // for pair, make_pair, move
#include <vector>       // for vector, vector<>::iterator

#include "modle/contacts.hpp"                    // for ContactMatrix
#include "modle/hdf5.hpp"                        // for construct_error_stack, has_dataset
#include "modle/suppress_compiler_warnings.hpp"  // for DISABLE_WARNING_POP, DISABLE_WARNI...
#include "modle/utils.hpp"                       // for throw_with_trace

namespace modle::cooler {

Cooler::Cooler(std::filesystem::path path_to_file, IO_MODE mode, size_t bin_size,
               size_t max_str_length, std::string_view assembly_name, Flavor flavor, bool validate,
               uint_fast8_t compression_lvl, size_t chunk_size, size_t cache_size)
    : STR_TYPE(generate_default_str_type(max_str_length)),
      _path_to_file(path_to_file),
      _mode(mode),
      _bin_size(bin_size),
      _assembly_name(assembly_name.data(), assembly_name.size()),
      _flavor(flavor),
      _fp(open_file(_path_to_file, _mode, _bin_size, max_str_length, _flavor, validate)),
      _groups(open_groups(*_fp, !this->is_read_only(), this->_bin_size)),
      _compression_lvl(compression_lvl),
      _chunk_size(chunk_size),
      _cache_size(cache_size),
      _cprop_str(this->is_read_only() ? nullptr
                                      : generate_default_cprop(_chunk_size, _compression_lvl,
                                                               Cooler::STR_TYPE, "\0")),
      _cprop_int32(this->is_read_only() ? nullptr
                                        : generate_default_cprop(_chunk_size, _compression_lvl,
                                                                 Cooler::INT32_TYPE, 0)),
      _cprop_int64(this->is_read_only() ? nullptr
                                        : generate_default_cprop(_chunk_size, _compression_lvl,
                                                                 Cooler::INT64_TYPE, 0L)),
      _aprop_str(generate_default_aprop(Cooler::STR_TYPE, _chunk_size, _cache_size)),
      _aprop_int32(generate_default_aprop(Cooler::INT32_TYPE, _chunk_size, _cache_size)),
      _aprop_int64(generate_default_aprop(Cooler::INT64_TYPE, _chunk_size, _cache_size)),
      _aprop_float64(
          generate_default_aprop(H5::PredType::NATIVE_DOUBLE, _chunk_size, _cache_size)) {
  assert(this->_flavor != UNK);
  if (this->_mode == READ_ONLY && this->_flavor == AUTO) {
    this->_flavor = Cooler::detect_file_flavor(*this->_fp);
  }
  if (this->_flavor == MCOOL) {
    assert(this->_bin_size != 0);
    absl::StrAppend(&this->_root_path, "resolutions/", this->_bin_size, "/");
  }
  if (this->is_read_only()) {
    if (this->_bin_size == 0) {  // i.e. file is cooler
      this->_bin_size =
          static_cast<size_t>(hdf5::read_attribute_int(*this->_fp, "bin-size", this->_root_path));
    }
    this->open_default_datasets();
    this->read_chrom_offset_idx();
    this->read_bin1_offset_idx();
  } else {
    this->init_default_datasets();
    this->write_metadata();
  }
}

Cooler::Cooler(const std::string &path_to_file, IO_MODE mode, size_t bin_size,
               size_t max_str_length, std::string_view assembly_name, Flavor flavor, bool validate,
               uint_fast8_t compression_lvl, size_t chunk_size, size_t cache_size)
    : Cooler::Cooler(std::filesystem::path(path_to_file), mode, bin_size, max_str_length,
                     assembly_name, flavor, validate, compression_lvl, chunk_size, cache_size) {}

Cooler::Cooler(std::string_view path_to_file, IO_MODE mode, size_t bin_size, size_t max_str_length,
               std::string_view assembly_name, Flavor flavor, bool validate,
               uint_fast8_t compression_lvl, size_t chunk_size, size_t cache_size)
    : Cooler::Cooler(std::filesystem::path(path_to_file), mode, bin_size, max_str_length,
                     assembly_name, flavor, validate, compression_lvl, chunk_size, cache_size) {}

Cooler::~Cooler() {
  try {
    if (this->_mode == WRITE_ONLY && this->_nchroms != 0 && this->_nbins != 0) {
      if (this->_fp) {
        auto &chrom_idx = this->_datasets[IDX_CHR];
        auto &bin1_idx = this->_datasets[IDX_BIN1];
        auto chrom_idx_offset = this->_dataset_file_offsets[IDX_CHR];
        auto bin1_idx_offset = this->_dataset_file_offsets[IDX_BIN1];

        const auto nbins = ++this->_nbins;
        (void)hdf5::write_number(nbins, chrom_idx, chrom_idx_offset);
        const auto nnz = ++this->_nnz;
        (void)hdf5::write_number(nnz, bin1_idx, bin1_idx_offset);
      } else {
        fmt::print(
            stderr,
            FMT_STRING("WARNING: Message for the developers: ~Cooler() for file '{}' was "
                       "called on a file that is supposed to be opened in write-mode, but the "
                       "file handle is actually already closed!"),
            this->_path_to_file);
      }
    }
  } catch (const H5::Exception &err) {
    fmt::print(stderr, FMT_STRING("The following error occurred while finalizing file '{}':\n{}"),
               this->_path_to_file, hdf5::construct_error_stack());
  }
}

bool Cooler::is_read_only() const { return this->_mode == Cooler::READ_ONLY; }

const std::filesystem::path &Cooler::get_path() const { return this->_path_to_file; }

size_t Cooler::get_nchroms() {
  assert(this->_fp);
  if (this->is_cool()) {
    return static_cast<size_t>(hdf5::read_attribute_int(*this->_fp, "nchroms"));
  }
  if (this->is_mcool()) {
    assert(this->_bin_size != 0);
    return static_cast<size_t>(hdf5::read_attribute_int(
        *this->_fp, "nchroms", absl::StrCat("/resolutions/", this->_bin_size)));
  }
  utils::throw_with_trace(std::logic_error("Unreachable code"));
}

void Cooler::get_chrom_names(std::vector<std::string> &buff) {
  assert(this->_fp);
  const auto nchroms = this->get_nchroms();
  buff.resize(nchroms);
  if (!this->_datasets.empty()) {
    const auto &d = this->_datasets[chrom_NAME];
    (void)hdf5::read_strings(d, buff, 0);
  } else {
    if (this->is_cool()) {
      auto d = this->_fp->openDataSet("/chroms/name", *this->_aprop_str);
      (void)hdf5::read_strings(d, buff, 0);
    } else if (this->is_mcool()) {
      assert(this->_bin_size != 0);
      auto d = this->_fp->openDataSet(
          absl::StrCat("/resolutions/", this->_bin_size, "/chroms/name"), *this->_aprop_str);
      (void)hdf5::read_strings(d, buff, 0);
    } else {
      assert(!this->is_scool());
      buff.clear();
    }
  }
}

std::vector<std::string> Cooler::get_chrom_names() {
  std::vector<std::string> buff;
  this->get_chrom_names(buff);
  return buff;
}

std::vector<std::pair<std::string, size_t>> Cooler::get_chroms() {
  std::vector<std::string> name_buff;
  std::vector<size_t> size_buff;
  this->get_chrom_names(name_buff);
  this->get_chrom_sizes(size_buff);
  assert(name_buff.size() == size_buff.size());  // NOLINT
  std::vector<std::pair<std::string, size_t>> buff(name_buff.size());
  for (auto i = 0UL; i < buff.size(); ++i) {
    buff[i].first = std::move(name_buff[i]);
    buff[i].second = size_buff[i];
  }
  return buff;
}

std::vector<int64_t> Cooler::get_chrom_sizes() {
  std::vector<int64_t> buff;
  this->get_chrom_sizes(buff);
  return buff;
}

bool Cooler::is_cool() const { return this->_flavor == COOL; }
bool Cooler::is_mcool() const { return this->_flavor == MCOOL; }
bool Cooler::is_scool() const { return this->_flavor == SCOOL; }

size_t Cooler::get_bin_size() const { return this->_bin_size; }

bool Cooler::has_contacts_for_chrom(std::string_view chrom_name, bool try_common_chrom_prefixes) {
  assert(this->_fp);  // NOLINT
  const auto chrom_idx = this->get_chrom_idx(chrom_name, try_common_chrom_prefixes);
  return this->has_contacts_for_chrom(chrom_idx);
}

bool Cooler::has_contacts_for_chrom(size_t chrom_idx) const {
  assert(this->_fp);                                   // NOLINT
  assert(this->is_read_only());                        // NOLINT
  assert(!this->_idx_bin1_offset.empty());             // NOLINT
  assert(!this->_idx_chrom_offset.empty());            // NOLINT
  assert(chrom_idx < this->_idx_chrom_offset.size());  // NOLINT

  const auto first_bin = static_cast<size_t>(this->_idx_chrom_offset[chrom_idx]);
  const auto last_bin = static_cast<size_t>(this->_idx_chrom_offset[chrom_idx + 1]);
  assert(last_bin >= first_bin);  // NOLINT

  return this->_idx_bin1_offset[first_bin] != this->_idx_bin1_offset[last_bin];
}

void Cooler::write_metadata() {
  if (this->is_read_only()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Caught attempt to write metadata to an HDF5 file that is opened in "
                               "read-only mode. File name: {}"),
                    this->_path_to_file.string()));
  }
  assert(this->_bin_size != 0);
  H5::DataSpace attr_space(H5S_SCALAR);
  int64_t int_buff{};
  std::string str_buff{};
  std::string name{};

  try {
    name = "format";
    str_buff = "HDF5::Cooler";
    hdf5::write_or_create_attribute(*this->_fp, name, str_buff);

    name = "format-version";
    int_buff = 3;
    hdf5::write_or_create_attribute(*this->_fp, name, int_buff);

    name = "bin-type";
    str_buff = "fixed";
    hdf5::write_or_create_attribute(*this->_fp, name, str_buff);

    name = "bin-size";
    int_buff = static_cast<int64_t>(this->_bin_size);
    hdf5::write_or_create_attribute(*this->_fp, name, int_buff);

    name = "storage-mode";
    str_buff = "symmetric-upper";
    hdf5::write_or_create_attribute(*this->_fp, name, str_buff);

    if (!this->_assembly_name.empty()) {
      name = "assembly-name";
      str_buff = this->_assembly_name;
      hdf5::write_or_create_attribute(*this->_fp, name, str_buff);
    }

    name = "generated-by";
    str_buff =
        fmt::format(FMT_STRING("ModLE-v{}.{}.{}"), 0, 0, 1);  // TODO make ModLE ver a tunable
    hdf5::write_or_create_attribute(*this->_fp, name, str_buff);

    name = "creation-date";
    str_buff = absl::FormatTime(absl::Now(), absl::UTCTimeZone());
    hdf5::write_or_create_attribute(*this->_fp, name, str_buff);
  } catch (const H5::Exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("The following error occurred while writing metadata to file {}: "
                               "error while writing attribute '{}':\n{}"),
                    this->_path_to_file, name, hdf5::construct_error_stack()));
  }
}

ContactMatrix<uint32_t> Cooler::cooler_to_cmatrix(std::string_view chrom_name, size_t nrows,
                                                  std::pair<size_t, size_t> chrom_boundaries,
                                                  bool try_common_chrom_prefixes,
                                                  bool prefer_using_balanced_counts) {
  assert(this->_fp);                         // NOLINT
  assert(!this->_datasets.empty());          // NOLINT
  assert(!this->_idx_bin1_offset.empty());   // NOLINT
  assert(!this->_idx_chrom_offset.empty());  // NOLINT

  const auto chrom_idx = this->get_chrom_idx(chrom_name, try_common_chrom_prefixes);
  const auto chrom_size = static_cast<size_t>(this->get_chrom_sizes()[chrom_idx]);
  if (chrom_boundaries.second > chrom_size) {
    chrom_boundaries.second = chrom_size;
  }
  const auto bin_range = std::make_pair(static_cast<size_t>(this->_idx_chrom_offset[chrom_idx]) +
                                            (chrom_boundaries.first / this->_bin_size),
                                        static_cast<size_t>(this->_idx_chrom_offset[chrom_idx]) +
                                            (chrom_boundaries.second / this->_bin_size));
  const auto bin1_offset_idx = this->get_bin1_offset_idx_for_chrom(chrom_idx, chrom_boundaries);
  double pxl_count_scaling_factor{1.0};

  if (prefer_using_balanced_counts &&
      hdf5::has_dataset(*this->_fp, "bins/weight", this->_root_path)) {
    const auto &d = this->_datasets[BIN_WEIGHT];
    uint8_t cis_only;
    try {
      hdf5::read_attribute(d, "cis_only", cis_only);
    } catch (const std::runtime_error &e) {
      if (absl::StrContains(e.what(), "cis_only")) {
        throw std::runtime_error(
            "File has a \"bins/weight\" dataset, but does not have an attribute named "
            "\"cis_only\". This most likely means that the file was generated by a very old "
            "version of Cooler. In order to proceed, you should rebalance the file using a "
            "recent "
            "version of Cooler, or in alternative remove the weight dataset (not recommended).");
      }
    }
    if (cis_only) {  // --cis-only balancing produces an array of sale factors
      std::vector<double> buff;
      hdf5::read_attribute(d, "scale", buff);
      pxl_count_scaling_factor = buff[chrom_idx];
    } else {  // standard or --trans-only balancing produces a single scale factor
      hdf5::read_attribute(d, "scale", pxl_count_scaling_factor);
    }
  }

  return cooler_to_cmatrix(bin_range, bin1_offset_idx, nrows, pxl_count_scaling_factor,
                           prefer_using_balanced_counts);
}

ContactMatrix<uint32_t> Cooler::cooler_to_cmatrix(std::string_view chrom_name,
                                                  size_t diagonal_width, size_t bin_size,
                                                  std::pair<size_t, size_t> chrom_boundaries,
                                                  bool try_common_chrom_prefixes,
                                                  bool prefer_using_balanced_counts) {
  assert(this->_bin_size != 0);
  if (bin_size != 0 && this->_bin_size != bin_size) {
    throw std::runtime_error(fmt::format(
        FMT_STRING(
            "Unable to read a Cooler file with bin size {} in a contact matrix with bin size {}"),
        this->_bin_size, bin_size));
  }
  const auto nrows = (diagonal_width / this->_bin_size) +
                     static_cast<size_t>(diagonal_width % this->_bin_size != 0);
  return cooler_to_cmatrix(chrom_name, nrows, chrom_boundaries, try_common_chrom_prefixes,
                           prefer_using_balanced_counts);
}

size_t Cooler::stream_contacts_for_chrom(
    moodycamel::BlockingReaderWriterQueue<Cooler::Pixel> &queue, std::string_view chrom_name,
    size_t nrows, std::pair<size_t, size_t> chrom_boundaries, bool try_common_chrom_prefixes,
    bool prefer_using_balanced_counts) {
  assert(this->_fp);                         // NOLINT
  assert(!this->_datasets.empty());          // NOLINT
  assert(!this->_idx_bin1_offset.empty());   // NOLINT
  assert(!this->_idx_chrom_offset.empty());  // NOLINT

  const auto chrom_idx = this->get_chrom_idx(chrom_name, try_common_chrom_prefixes);
  const auto chrom_size = static_cast<size_t>(this->get_chrom_sizes()[chrom_idx]);
  if (chrom_boundaries.second > chrom_size) {
    chrom_boundaries.second = chrom_size;
  }
  const auto bin_range = std::make_pair(static_cast<size_t>(this->_idx_chrom_offset[chrom_idx]) +
                                            (chrom_boundaries.first / this->_bin_size),
                                        static_cast<size_t>(this->_idx_chrom_offset[chrom_idx]) +
                                            (chrom_boundaries.second / this->_bin_size));
  const auto bin1_offset_idx = this->get_bin1_offset_idx_for_chrom(chrom_idx, chrom_boundaries);
  double pxl_count_scaling_factor{1.0};

  if (prefer_using_balanced_counts &&
      hdf5::has_dataset(*this->_fp, "bins/weight", this->_root_path)) {
    const auto &d = this->_datasets[BIN_WEIGHT];
    uint8_t cis_only;
    try {
      hdf5::read_attribute(d, "cis_only", cis_only);
    } catch (const std::runtime_error &e) {
      if (absl::StrContains(e.what(), "cis_only")) {
        throw std::runtime_error(
            "File has a \"bins/weight\" dataset, but does not have an attribute named "
            "\"cis_only\". This most likely means that the file was generated by a very old "
            "version of Cooler. In order to proceed, you should rebalance the file using a "
            "recent "
            "version of Cooler, or in alternative remove the weight dataset (not recommended).");
      }
    }
    if (cis_only) {  // --cis-only balancing produces an array of sale factors
      std::vector<double> buff;
      hdf5::read_attribute(d, "scale", buff);
      pxl_count_scaling_factor = buff[chrom_idx];
    } else {  // standard or --trans-only balancing produces a single scale factor
      hdf5::read_attribute(d, "scale", pxl_count_scaling_factor);
    }
  }

  return this->stream_contacts_for_chrom(queue, bin_range, bin1_offset_idx, nrows,
                                         pxl_count_scaling_factor, prefer_using_balanced_counts);
}

size_t Cooler::stream_contacts_for_chrom(
    moodycamel::BlockingReaderWriterQueue<Cooler::Pixel> &queue, std::string_view chrom_name,
    size_t diagonal_width, size_t bin_size, std::pair<size_t, size_t> chrom_boundaries,
    bool try_common_chrom_prefixes, bool prefer_using_balanced_counts) {
  const auto nrows =
      (diagonal_width / bin_size) + static_cast<size_t>(diagonal_width % bin_size == 0);
  return this->stream_contacts_for_chrom(queue, chrom_name, nrows, chrom_boundaries,
                                         try_common_chrom_prefixes, prefer_using_balanced_counts);
}

bool Cooler::validate_file_format(H5::H5File &f, Flavor expected_flavor, IO_MODE mode,
                                  size_t bin_size, bool throw_on_failure) {
  if (mode == WRITE_ONLY) {  // File is empty. Nothing to validate here
    return true;
  }

  try {
    const auto flavor = detect_file_flavor(f);
    if (expected_flavor != AUTO && expected_flavor != UNK && expected_flavor != flavor) {
      throw std::runtime_error(fmt::format(FMT_STRING("Expected format flavor {}, found {}"),
                                           flavor_to_string(expected_flavor),
                                           flavor_to_string(flavor)));
    }
    switch (flavor) {
      case COOL:
        return validate_cool_flavor(f, bin_size, "/", throw_on_failure);
      case MCOOL:
        return validate_multires_cool_flavor(f, bin_size, "/", throw_on_failure);
      case SCOOL:
        throw std::runtime_error("SCOOL flavor is not yet supported");
      default:
        assert(false);  // This code should be unreachable
        return false;
    }

  } catch (const std::runtime_error &e) {
    throw std::runtime_error(fmt::format(FMT_STRING("Cooler validation for file '{}' failed: {}"),
                                         f.getFileName(), e.what()));
  }
}

H5::StrType Cooler::generate_default_str_type(size_t max_str_length) {
  // Cooltools doesn't seem to properly handle variable length strings (H5T_VARIABLE)
  // For the time being we are forced to use fixed length, null-padded strings
  auto st = max_str_length > 0 ? H5::StrType(H5::PredType::C_S1, max_str_length)
                               : H5::StrType(H5::PredType::C_S1);
  st.setStrpad(H5T_STR_NULLPAD);
  st.setCset(H5T_CSET_ASCII);
  return st;
}

std::unique_ptr<H5::H5File> Cooler::open_file(const std::filesystem::path &path, IO_MODE mode,
                                              size_t bin_size, size_t max_str_length, Flavor flavor,
                                              bool validate) {
  (void)max_str_length;
#ifndef NDEBUG
  if (mode == WRITE_ONLY) {
    if (bin_size == 0) {
      utils::throw_with_trace(
          std::runtime_error("Cooler::open_file(): bin_size cannot be 0 when file is being opened "
                             "in WRITE_ONLY mode"));
    }
    if (max_str_length == 0) {
      utils::throw_with_trace(
          std::runtime_error("Cooler::open_file(): max_str_length cannot be 0 when file is being "
                             "opened in WRITE_ONLY mode"));
    }
  }
#endif
  try {
    H5::H5File f(path.c_str(), mode == READ_ONLY ? H5F_ACC_RDONLY : H5F_ACC_TRUNC);
    if (validate) {
      (void)validate_file_format(f, flavor, mode, bin_size, true);
    }
    return std::make_unique<H5::H5File>(f);
  } catch (const H5::Exception &e) {
    throw std::runtime_error(hdf5::construct_error_stack());
  }
}

std::vector<H5::Group> Cooler::open_groups(H5::H5File &f, bool create_if_not_exist,
                                           size_t bin_size) {
  std::vector<H5::Group> groups(4);  // NOLINT

  auto open_or_create_group = [&](Cooler::Groups g, const std::string &group_name) {
    groups[g] = create_if_not_exist && !f.nameExists(group_name) ? f.createGroup(group_name)
                                                                 : f.openGroup(group_name);
  };
  const std::string root_path = bin_size != 0 && hdf5::has_group(f, "/resolutions")
                                    ? absl::StrCat("/resolutions/", bin_size, "/")
                                    : "";
  try {
    open_or_create_group(chrom, absl::StrCat(root_path, "chroms"));
    open_or_create_group(BIN, absl::StrCat(root_path, "bins"));
    open_or_create_group(PXL, absl::StrCat(root_path, "pixels"));
    open_or_create_group(IDX, absl::StrCat(root_path, "indexes"));

    return groups;
  } catch (const H5::Exception &e) {
    throw std::runtime_error(hdf5::construct_error_stack());
  }
}

void Cooler::init_default_datasets() {
  constexpr hsize_t RANK{1};                 // i.e. number of dimensions
  constexpr hsize_t MAXDIMS{H5S_UNLIMITED};  // extensible dataset
  constexpr hsize_t BUFF_SIZE{1};            // Dummy buffer size

  this->_datasets.resize(DEFAULT_DATASETS_NR);
  this->_dataset_file_offsets.resize(DEFAULT_DATASETS_NR);
  std::fill(this->_dataset_file_offsets.begin(), this->_dataset_file_offsets.end(), 0);

  assert(this->_cprop_str);
  assert(this->_cprop_int32);
  assert(this->_cprop_int64);

  assert(this->_aprop_str);
  assert(this->_aprop_int32);
  assert(this->_aprop_int64);

  const auto &INT64 = this->INT64_TYPE;
  const auto &INT32 = this->INT32_TYPE;
  const auto &STR = this->STR_TYPE;

  const auto &cp64 = *this->_cprop_int64;
  const auto &cp32 = *this->_cprop_int32;
  const auto &cps = *this->_cprop_str;

  const auto &ap64 = *this->_aprop_int64;
  const auto &ap32 = *this->_aprop_int32;
  const auto &aps = *this->_aprop_str;

  auto &f = *this->_fp;
  auto &dset = this->_datasets;
  auto r = absl::StripSuffix(this->_root_path, "/");

  try {
    auto mem_space{H5::DataSpace(RANK, &BUFF_SIZE, &MAXDIMS)};
    // Do not change the order of these pushbacks
    dset[chrom_LEN] =
        f.createDataSet(absl::StrCat(r, "/chroms/length"), INT64, mem_space, cp64, ap64);
    dset[chrom_NAME] = f.createDataSet(absl::StrCat(r, "/chroms/name"), STR, mem_space, cps, aps);

    dset[BIN_CHROM] = f.createDataSet(absl::StrCat(r, "/bins/chrom"), INT64, mem_space, cp64, ap64);
    dset[BIN_START] = f.createDataSet(absl::StrCat(r, "/bins/start"), INT64, mem_space, cp64, ap64);
    dset[BIN_END] = f.createDataSet(absl::StrCat(r, "/bins/end"), INT64, mem_space, cp64, ap64);

    dset[PXL_B1] =
        f.createDataSet(absl::StrCat(r, "/pixels/bin1_id"), INT64, mem_space, cp64, ap64);
    dset[PXL_B2] =
        f.createDataSet(absl::StrCat(r, "/pixels/bin2_id"), INT64, mem_space, cp64, ap64);
    dset[PXL_COUNT] =
        f.createDataSet(absl::StrCat(r, "/pixels/count"), INT32, mem_space, cp32, ap32);

    dset[IDX_BIN1] =
        f.createDataSet(absl::StrCat(r, "/indexes/bin1_offset"), INT64, mem_space, cp64, ap64);
    dset[IDX_CHR] =
        f.createDataSet(absl::StrCat(r, "/indexes/chrom_offset"), INT64, mem_space, cp64, ap64);

  } catch (const H5::FileIException &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("An error occurred while initializing default Cooler dataset on file '{}': {}"),
        this->_path_to_file.string(), hdf5::construct_error_stack()));
  }
}

void Cooler::open_default_datasets() {
  assert(this->_fp);
  assert(this->_aprop_int32);
  assert(this->_aprop_int64);
  assert(this->_aprop_float64);
  assert(this->_aprop_str);

  auto &d = this->_datasets;
  auto &f = *this->_fp;
  const auto &ai32 = *this->_aprop_int32;
  const auto &ai64 = *this->_aprop_int64;
  const auto &af64 = *this->_aprop_float64;
  const auto &as = *this->_aprop_str;

  this->_datasets.resize(DEFAULT_DATASETS_NR);
  this->_dataset_file_offsets.resize(DEFAULT_DATASETS_NR);
  std::fill(this->_dataset_file_offsets.begin(), this->_dataset_file_offsets.end(), 0);
  try {
    // Do not change the order of these pushbacks
    d[chrom_LEN] = f.openDataSet(absl::StrCat(this->_root_path, "chroms/length"), ai64);
    d[chrom_NAME] = f.openDataSet(absl::StrCat(this->_root_path, "chroms/name"), as);

    d[BIN_CHROM] = f.openDataSet(absl::StrCat(this->_root_path, "bins/chrom"), ai64);
    d[BIN_START] = f.openDataSet(absl::StrCat(this->_root_path, "bins/start"), ai64);
    d[BIN_END] = f.openDataSet(absl::StrCat(this->_root_path, "bins/end"), ai64);
    if (hdf5::has_dataset(f, "bins/weight", this->_root_path)) {
      d[BIN_WEIGHT] = f.openDataSet(absl::StrCat(this->_root_path, "bins/weight"), af64);
    }

    d[PXL_B1] = f.openDataSet(absl::StrCat(this->_root_path, "pixels/bin1_id"), ai64);
    d[PXL_B2] = f.openDataSet(absl::StrCat(this->_root_path, "pixels/bin2_id"), ai64);
    d[PXL_COUNT] = f.openDataSet(absl::StrCat(this->_root_path, "pixels/count"), ai32);

    d[IDX_BIN1] = f.openDataSet(absl::StrCat(this->_root_path, "indexes/bin1_offset"), ai64);
    d[IDX_CHR] = f.openDataSet(absl::StrCat(this->_root_path, "indexes/chrom_offset"), ai64);

  } catch (const H5::FileIException &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING(
            "An error occurred while trying to open Cooler's default datasets of file '{}': {}"),
        this->_path_to_file, hdf5::construct_error_stack()));
  }
}

size_t Cooler::read_chrom_offset_idx() {
  const auto &d = this->_datasets[IDX_CHR];
  const auto buff_size = static_cast<hsize_t>(d.getSpace().getSimpleExtentNpoints());
  this->_idx_chrom_offset.resize(buff_size);

  const auto idx_size = hdf5::read_numbers(d, this->_idx_chrom_offset, 0);
  if (idx_size != this->_idx_chrom_offset.size()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("An error occurred while reading dataset 'indexes/chrom_offset' "
                               "from file '{}': expected to read {} numbers, but only read {}"),
                    this->_path_to_file, buff_size, idx_size));
  }
  return idx_size;
}

size_t Cooler::read_bin1_offset_idx() {
  const auto &d = this->_datasets[IDX_BIN1];
  const auto buff_size = static_cast<hsize_t>(d.getSpace().getSimpleExtentNpoints());
  this->_idx_bin1_offset.resize(buff_size);

  const auto idx_size = hdf5::read_numbers(d, this->_idx_bin1_offset, 0);
  if (idx_size != this->_idx_bin1_offset.size()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("An error occurred while reading dataset 'indexes/chrom_offset' "
                               "from file '{}': expected to read {} numbers, but only read {}"),
                    this->_path_to_file, buff_size, idx_size));
  }
  return idx_size;
}

absl::Span<const int64_t> Cooler::get_bin1_offset_idx_for_chrom(
    size_t chrom_idx, std::pair<size_t, size_t> chrom_subrange) {
  assert(!this->_idx_bin1_offset.empty());             // NOLINT
  assert(!this->_idx_chrom_offset.empty());            // NOLINT
  assert(chrom_idx < this->_idx_chrom_offset.size());  // NOLINT
  const auto chrom_start_bin = static_cast<size_t>(this->_idx_chrom_offset[chrom_idx]) +
                               (chrom_subrange.first / this->_bin_size);
  const auto chrom_end_bin =
      chrom_subrange.second == std::numeric_limits<decltype(chrom_subrange.second)>::max()
          ? static_cast<size_t>(this->_idx_chrom_offset[chrom_idx + 1])
          : static_cast<size_t>(this->_idx_chrom_offset[chrom_idx]) +
                (chrom_subrange.second / this->_bin_size);
  DISABLE_WARNING_PUSH
  DISABLE_WARNING_SIGN_COMPARE
  assert(chrom_end_bin <= this->_idx_chrom_offset[chrom_idx + 1]);
  assert(chrom_end_bin >= chrom_start_bin);  // NOLINT
  DISABLE_WARNING_POP

  return absl::MakeConstSpan(this->_idx_bin1_offset)
      .subspan(chrom_start_bin, chrom_end_bin - chrom_start_bin);
}

absl::Span<const int64_t> Cooler::get_bin1_offset_idx_for_chrom(
    std::string_view chrom_name, std::pair<size_t, size_t> chrom_subrange) {
  const auto chrom_idx = get_chrom_idx(chrom_name);
  return get_bin1_offset_idx_for_chrom(chrom_idx, chrom_subrange);
}

std::pair<int64_t, int64_t> Cooler::read_chrom_pixels_boundaries(std::string_view chrom_name) {
  const auto chrom_idx = get_chrom_idx(chrom_name);
  return read_chrom_pixels_boundaries(chrom_idx);
}

std::pair<int64_t, int64_t> Cooler::read_chrom_pixels_boundaries(size_t chrom_idx) {
  assert(chrom_idx < this->_idx_chrom_offset.size());  // NOLINT
  assert(!this->_idx_chrom_offset.empty());            // NOLINT
  const auto chrom_start_bin = static_cast<size_t>(this->_idx_chrom_offset[chrom_idx]);
  const auto chrom_end_bin = static_cast<size_t>(this->_idx_chrom_offset[chrom_idx + 1]);
  assert(chrom_end_bin >= chrom_start_bin);  // NOLINT

  std::pair<int64_t, int64_t> pixel_boundaries{};
  const auto &d = this->_datasets[IDX_BIN1];
  (void)hdf5::read_number(d, pixel_boundaries.first, chrom_start_bin);
  (void)hdf5::read_number(d, pixel_boundaries.second, chrom_end_bin);

  return pixel_boundaries;
}

ContactMatrix<uint32_t> Cooler::cooler_to_cmatrix(std::pair<hsize_t, hsize_t> bin_range,
                                                  absl::Span<const int64_t> bin1_offset_idx,
                                                  size_t nrows, double bias_scaling_factor,
                                                  bool prefer_using_balanced_counts) {
  if (this->_datasets.empty()) {
    this->open_default_datasets();
  }

  const auto &[first_bin, last_bin] = bin_range;
  assert(first_bin < last_bin);                            // NOLINT
  assert(last_bin <= first_bin + bin1_offset_idx.size());  // NOLINT
  ContactMatrix<uint32_t> cmatrix(nrows, last_bin - first_bin + 1);
  std::vector<int64_t> bin1_BUFF(nrows);
  std::vector<int64_t> bin2_BUFF(nrows);
  std::vector<int64_t> count_BUFF(nrows);
  std::vector<double> bin_weights;

  const auto &d = this->_datasets;
  if (prefer_using_balanced_counts &&
      hdf5::has_dataset(*this->_fp, "bins/weight", this->_root_path)) {
    bin_weights.resize(last_bin - first_bin + 1);
    (void)hdf5::read_numbers(d[BIN_WEIGHT], bin_weights, static_cast<hsize_t>(first_bin));
  }

  for (auto i = 1UL; i < bin1_offset_idx.size(); ++i) {
    const auto file_offset = static_cast<hsize_t>(bin1_offset_idx[i - 1]);
    const auto buff_size =
        std::min(static_cast<size_t>(bin1_offset_idx[i] - bin1_offset_idx[i - 1]), nrows);
    if (buff_size == 0) {
      continue;
    }

    bin1_BUFF.resize(buff_size);
    bin2_BUFF.resize(buff_size);
    count_BUFF.resize(buff_size);

    (void)hdf5::read_numbers(d[PXL_B1], bin1_BUFF, file_offset);
    (void)hdf5::read_numbers(d[PXL_B2], bin2_BUFF, file_offset);
    (void)hdf5::read_numbers(d[PXL_COUNT], count_BUFF, file_offset);

    assert(bin1_BUFF.size() == buff_size);   // NOLINT
    assert(bin2_BUFF.size() == buff_size);   // NOLINT
    assert(count_BUFF.size() == buff_size);  // NOLINT

    for (auto j = 0UL; j < buff_size; ++j) {
      assert(count_BUFF[j] != 0);  // NOLINT
      DISABLE_WARNING_PUSH
      DISABLE_WARNING_SIGN_CONVERSION
      DISABLE_WARNING_SIGN_COMPARE
      DISABLE_WARNING_CONVERSION
      const auto bin1 = bin1_BUFF[j] - first_bin;
      const auto bin2 = bin2_BUFF[j] - first_bin;
      if (bin2 >= i + nrows - 1 || bin2 >= bin1_offset_idx.size()) {
        break;
      }
      // fmt::print(stderr, "m[{}][{}]={} (nrows={}; ncols={})\n", bin2_BUFF[j] - first_bin,
      //           bin1_BUFF[j] - first_bin, count_BUFF[j], cmatrix.nrows(),
      //           cmatrix.ncols());
      if (bin_weights.empty()) {
        cmatrix.set(bin2, bin1, static_cast<uint32_t>(count_BUFF[j]));
      } else {
        // According to Cooler documentations, NaN means that a bin has been excluded by the
        // matrix balancing procedure. In this case we set the count to 0
        if (std::isnan(bin_weights[bin1]) || std::isnan(bin_weights[bin2])) {
          continue;  // Same as setting count to 0;
        }
        const auto bin1_bias = bin_weights[bin1];
        const auto bin2_bias = bin_weights[bin2];
        // See https://github.com/robomics/modle/issues/36 and
        // https://github.com/open2c/cooler/issues/35
        const auto count =
            static_cast<double>(count_BUFF[j]) / (bin1_bias * bin2_bias) / bias_scaling_factor;
        cmatrix.set(bin2, bin1, static_cast<uint32_t>(std::round(count)));
        DISABLE_WARNING_POP
      }
    }
  }

  return cmatrix;
}

ContactMatrix<uint32_t> Cooler::cooler_to_cmatrix(std::pair<hsize_t, hsize_t> bin_range,
                                                  const std::vector<int64_t> &bin1_offset_idx,
                                                  size_t nrows, double scaling_factor,
                                                  bool prefer_using_balanced_counts) {
  return this->cooler_to_cmatrix(bin_range, absl::MakeConstSpan(bin1_offset_idx), nrows,
                                 scaling_factor, prefer_using_balanced_counts);
}

size_t Cooler::stream_contacts_for_chrom(moodycamel::BlockingReaderWriterQueue<Pixel> &queue,
                                         std::pair<hsize_t, hsize_t> bin_range,
                                         absl::Span<const int64_t> bin1_offset_idx, size_t nrows,
                                         double bias_scaling_factor,
                                         bool prefer_using_balanced_counts) {
  if (this->_datasets.empty()) {
    this->open_default_datasets();
  }

  const auto &[first_bin, last_bin] = bin_range;
  assert(first_bin < last_bin);                            // NOLINT
  assert(last_bin <= first_bin + bin1_offset_idx.size());  // NOLINT
  std::vector<int64_t> bin1_BUFF(nrows);
  std::vector<int64_t> bin2_BUFF(nrows);
  std::vector<int64_t> count_BUFF(nrows);
  std::vector<double> bin_weights;
  size_t pixel_count = 0;

  const auto &d = this->_datasets;
  if (prefer_using_balanced_counts &&
      hdf5::has_dataset(*this->_fp, "bins/weight", this->_root_path)) {
    bin_weights.resize(last_bin - first_bin + 1);
    (void)hdf5::read_numbers(d[BIN_WEIGHT], bin_weights, static_cast<hsize_t>(first_bin));
  }

  for (auto i = 1UL; i < bin1_offset_idx.size(); ++i) {
    const auto file_offset = static_cast<hsize_t>(bin1_offset_idx[i - 1]);
    const auto buff_size =
        std::min(static_cast<size_t>(bin1_offset_idx[i] - bin1_offset_idx[i - 1]), nrows);
    if (buff_size == 0) {
      continue;
    }

    bin1_BUFF.resize(buff_size);
    bin2_BUFF.resize(buff_size);
    count_BUFF.resize(buff_size);

    (void)hdf5::read_numbers(d[PXL_B1], bin1_BUFF, file_offset);
    (void)hdf5::read_numbers(d[PXL_B2], bin2_BUFF, file_offset);
    (void)hdf5::read_numbers(d[PXL_COUNT], count_BUFF, file_offset);

    assert(bin1_BUFF.size() == buff_size);   // NOLINT
    assert(bin2_BUFF.size() == buff_size);   // NOLINT
    assert(count_BUFF.size() == buff_size);  // NOLINT

    for (auto j = 0UL; j < buff_size; ++j) {
      assert(count_BUFF[j] != 0);  // NOLINT
      DISABLE_WARNING_PUSH
      DISABLE_WARNING_SIGN_CONVERSION
      DISABLE_WARNING_SIGN_COMPARE
      DISABLE_WARNING_CONVERSION
      const auto bin1 = bin1_BUFF[j] - first_bin;
      const auto bin2 = bin2_BUFF[j] - first_bin;
      if (bin2 >= i + nrows - 1 || bin2 >= bin1_offset_idx.size()) {
        break;
      }
      // fmt::print(stderr, "m[{}][{}]={} (nrows={}; ncols={})\n", bin2_BUFF[j] - first_bin,
      //           bin1_BUFF[j] - first_bin, count_BUFF[j], cmatrix.nrows(),
      //           cmatrix.ncols());
      if (bin_weights.empty()) {
        while (!queue.try_emplace(Pixel{std::min(bin1, bin2), std::max(bin1, bin2),
                                        static_cast<size_t>(count_BUFF[j])})) {
          std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }
        ++pixel_count;
      } else {
        // According to Cooler documentations, NaN means that a bin has been excluded by the
        // matrix balancing procedure. In this case we set the count to 0
        if (std::isnan(bin_weights[bin1]) || std::isnan(bin_weights[bin2])) {
          continue;  // Same as setting count to 0;
        }
        const auto bin1_bias = bin_weights[bin1];
        const auto bin2_bias = bin_weights[bin2];
        // See https://github.com/robomics/modle/issues/36 and
        // https://github.com/open2c/cooler/issues/35
        const auto count =
            static_cast<double>(count_BUFF[j]) / (bin1_bias * bin2_bias) / bias_scaling_factor;
        while (!queue.try_emplace(Pixel{std::min(bin1, bin2), std::max(bin1, bin2),
                                        static_cast<size_t>(std::round(count))})) {
          std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }
        ++pixel_count;
        DISABLE_WARNING_POP
      }
    }
  }
  return pixel_count;
}

size_t Cooler::stream_contacts_for_chrom(moodycamel::BlockingReaderWriterQueue<Pixel> &queue,
                                         std::pair<hsize_t, hsize_t> bin_range,
                                         const std::vector<int64_t> &bin1_offset_idx, size_t nrows,
                                         double scaling_factor, bool prefer_using_balanced_counts) {
  return this->stream_contacts_for_chrom(queue, bin_range, absl::MakeConstSpan(bin1_offset_idx),
                                         nrows, scaling_factor, prefer_using_balanced_counts);
}

size_t Cooler::get_chrom_idx(std::string_view query_chrom_name, bool try_common_chrom_prefixes) {
  // Here's the issue: for a given genome assembly (say hg19), some tools name chromosome as
  // chrN, where N is the chromosome number, while others only use the number N. The purpose
  // of this lambda is to try to guess few reasonably common prefixes, look them up in the
  // Cooler file, then return the chromosome name variant that produced a hit together with
  // the chrom index

  try {
    const auto chrom_names = hdf5::read_strings(this->_datasets[chrom_NAME], 0);
    auto match = std::find(chrom_names.begin(), chrom_names.end(), query_chrom_name);
    if (match != chrom_names.end()) {
      return static_cast<size_t>(std::distance(chrom_names.begin(), match));
    }
    if (try_common_chrom_prefixes) {
      std::array<std::string_view, 3> queries;
      constexpr std::array<std::string_view, 3> prefixes = {"chrom", "chrom", "chrom"};
      for (auto i = 0U; i < prefixes.size(); ++i) {
        queries[i] = absl::StartsWith(query_chrom_name, prefixes[i])
                         ? absl::StripPrefix(query_chrom_name, prefixes[i])
                         : absl::StrCat(prefixes[i], query_chrom_name);
      }
      for (const auto &q : queries) {
        match = std::find(chrom_names.begin(), chrom_names.end(), q);
        if (match != chrom_names.end()) {
          return static_cast<size_t>(std::distance(chrom_names.begin(), match));
        }
      }
      throw std::runtime_error(
          fmt::format(FMT_STRING("Unable to find a chromosome named '{}'. The following chromosome "
                                 "name variants were searched: '{}'"),
                      query_chrom_name, absl::StrJoin(queries, "', '")));
    }
    throw std::runtime_error(
        fmt::format(FMT_STRING("Unable to find a chromosome named '{}'"), query_chrom_name));
  } catch (const std::exception &e) {
    if (absl::StartsWith(e.what(), "Unable to find a chromosome")) {
      throw;
    }
    throw std::runtime_error(fmt::format(FMT_STRING("The following error occurred while looking up "
                                                    "'{}' in dataset chroms/name of file '{}': {}"),
                                         query_chrom_name, this->_path_to_file, e.what()));
  }
}

std::string Cooler::flavor_to_string(Flavor f) {
  switch (f) {
    case UNK:
      return "unknown";
    case AUTO:
      return "auto";
    case COOL:
      return "COOL";
    case MCOOL:
      return "MCOOL";
    case SCOOL:
      return "SCOOL";
  }
  return "";
}

Cooler::Flavor Cooler::detect_file_flavor(H5::H5File &f) {
  if (const auto &fname = f.getFileName(); absl::EndsWithIgnoreCase(fname, ".cool")) {
    return COOL;
  } else if (absl::EndsWithIgnoreCase(fname, ".mcool")) {
    return MCOOL;
  } else if (absl::EndsWithIgnoreCase(fname, ".scool")) {
    return SCOOL;
  }

  std::string buff;
  if (hdf5::has_attribute(f, "format")) {
    hdf5::read_attribute(f, "format", buff);
  }
  if (const auto &fname = f.getFileName(); absl::EndsWithIgnoreCase(fname, "::cooler")) {
    return COOL;
  } else if (absl::EndsWithIgnoreCase(fname, "::mcool")) {
    return MCOOL;
  } else if (absl::EndsWithIgnoreCase(fname, "::scool")) {
    return SCOOL;
  }

  if (f.nameExists("/resolutions")) {
    return MCOOL;
  } else if (f.nameExists("/bins") && f.nameExists("/chroms")) {
    if (f.nameExists("/pixels") && f.nameExists("/indexes")) {
      return COOL;
    } else if (f.nameExists("/cells")) {
      return SCOOL;
    }
  }
  throw std::runtime_error("Unable to detect Cooler file flavor");
}

bool Cooler::validate_cool_flavor(H5::H5File &f, size_t bin_size, std::string_view root_path,
                                  bool throw_on_failure, bool check_version) {
  /* The following attributes are being checked:
   * format
   * format-version
   * bin-type
   * bin size
   * storage-mode
   */

  std::string str_buff{};
  int64_t int_buff{};

  constexpr int64_t min_format_ver = 2;
  constexpr int64_t max_format_ver = 3;

  try {
    // Checking attributes
    if (hdf5::has_attribute(f, "format", root_path)) {
      hdf5::read_attribute(f, "format", str_buff, root_path);
      if (!absl::EndsWithIgnoreCase(str_buff, "::cooler")) {
        if (!throw_on_failure) {
          return false;
        }
        throw std::runtime_error(fmt::format(
            FMT_STRING(
                "File is not in Cooler format: format attribute should be HDF5::Cooler, is '{}'"),
            str_buff));
      }
    } else {
      fmt::print(stderr, FMT_STRING("WARNING: missing attribute 'format' in file '{}'\n"),
                 f.getFileName());
    }
    str_buff.clear();

    if (check_version) {
      hdf5::read_attribute(f, "format-version", int_buff, root_path);
      if (int_buff < min_format_ver || int_buff > max_format_ver) {
        if (!throw_on_failure) {
          return false;
        }
        throw std::runtime_error(fmt::format(
            FMT_STRING("Expected format-version attribute to be between {} and {}, got {}"),
            min_format_ver, max_format_ver, int_buff));
      }
    }
    const auto format_version = static_cast<uint8_t>(int_buff);
    int_buff = 0;

    hdf5::read_attribute(f, "bin-type", str_buff, root_path);
    if (str_buff != "fixed") {
      if (!throw_on_failure) {
        return false;
      }
      throw std::runtime_error(
          fmt::format(FMT_STRING("Expected bin-type attribute to be 'fixed', got '{}'"), str_buff));
    }
    str_buff.clear();

    if (bin_size != 0) {
      hdf5::read_attribute(f, "bin-size", int_buff, root_path);
      if (static_cast<int64_t>(bin_size) != int_buff) {
        if (!throw_on_failure) {
          return false;
        }
        throw std::runtime_error(fmt::format(
            FMT_STRING("Expected bin-size attribute to be {}, got {}"), bin_size, int_buff));
      }
    }
    int_buff = 0;

    if (format_version > 2) {
      hdf5::read_attribute(f, "storage-mode", str_buff, root_path);
      if (str_buff != "symmetric-upper") {
        if (!throw_on_failure) {
          return false;
        }
        throw std::runtime_error(fmt::format(
            FMT_STRING("Expected storage-mode attribute to be 'symmetric-upper', got '{}'"),
            str_buff));
      }
    }
    str_buff.clear();

    // Checking groups
    auto check_group = [&](std::string_view name) {
      if (!hdf5::has_group(f, name, root_path)) {
        if (!throw_on_failure) {
          return false;
        }
        throw std::runtime_error(fmt::format(FMT_STRING("Missing mandatory group '{}/{}'"),
                                             absl::StripPrefix(root_path, "/"), name));
      }
      return true;
    };

    bool ok = true;
    ok = ok && check_group("chroms");
    ok = ok && check_group("bins");
    ok = ok && check_group("pixels");
    ok = ok && check_group("indexes");

    if (!ok) {
      return false;
    }

    // Checking datasets
    auto check_dset = [&](std::string_view name) {
      if (!hdf5::has_dataset(f, "chroms/length", root_path)) {
        if (!throw_on_failure) {
          return false;
        }
        throw std::runtime_error(fmt::format(FMT_STRING("Missing mandatory dataset '{}/{}'"),
                                             absl::StripPrefix(root_path, "/"), name));
      }
      return true;
    };

    ok = ok && check_dset("chroms/length");
    ok = ok && check_dset("chroms/name");
    ok = ok && check_dset("bins/chrom");
    ok = ok && check_dset("bins/start");
    ok = ok && check_dset("bins/end");
    ok = ok && check_dset("pixels/bin1_id");
    ok = ok && check_dset("pixels/bin2_id");
    ok = ok && check_dset("pixels/count");
    ok = ok && check_dset("indexes/bin1_offset");
    ok = ok && check_dset("indexes/chrom_offset");

    if (!ok) {
      return false;
    }

  } catch (const std::runtime_error &e) {
    if (std::string_view msg{e.what()};
        absl::ConsumePrefix(&msg, "Unable to find an attribute named '")) {
      if (!throw_on_failure) {
        return false;
      }
      auto attr_name = msg.substr(0, msg.find_first_of("'"));
      throw std::runtime_error(
          fmt::format(FMT_STRING("Missing mandatory attribute '{}'{}"), attr_name,
                      root_path == "/" ? "" : absl::StrCat(" from path '", root_path, "'")));
    }
    throw;
  }
  return true;
}

bool Cooler::validate_multires_cool_flavor(H5::H5File &f, size_t bin_size,
                                           std::string_view root_path, bool throw_on_failure) {
  constexpr int64_t min_format_ver = 2;
  constexpr int64_t max_format_ver = 3;

  if (bin_size == 0) {
    throw std::runtime_error(
        "A bin size other than 0 is required when calling "
        "Cooler::validate_multires_cool_flavor()");
  }

  if (hdf5::has_attribute(f, "format", root_path)) {
    if (const auto format = hdf5::read_attribute_str(f, "format", root_path);
        !absl::EndsWithIgnoreCase(format, "::mcool")) {
      if (!throw_on_failure) {
        return false;
      }
      throw std::runtime_error(
          fmt::format(FMT_STRING("File is not in Multires-Cooler (MCool) format: format attribute "
                                 "should be HDF5::MCOOL, is '{}'"),
                      format));
    }
  } else {
    fmt::print(stderr, FMT_STRING("WARNING: missing attribute 'format' in file '{}'\n"),
               f.getFileName());
  }

  if (const auto format_ver = hdf5::read_attribute_int(f, "format-version", root_path);
      format_ver < min_format_ver || format_ver > max_format_ver) {
    if (!throw_on_failure) {
      return false;
    }
    throw std::runtime_error(
        fmt::format(FMT_STRING("Expected format-version attribute to be between {} and {}, got {}"),
                    min_format_ver, max_format_ver, format_ver));
  }

  if (!hdf5::has_group(f, "resolutions", root_path)) {
    if (!throw_on_failure) {
      return false;
    }
    throw std::runtime_error(fmt::format(FMT_STRING("Missing mandatory group '{}/resolutions'"),
                                         absl::StripPrefix(root_path, "/")));
  }

  if (!hdf5::has_group(f, absl::StrCat("resolutions/", bin_size), root_path)) {
    if (!throw_on_failure) {
      return false;
    }
    throw std::runtime_error(
        fmt::format(FMT_STRING("Missing group for resolution {} bp"), bin_size));
  }

  return validate_cool_flavor(
      f, bin_size, absl::StrCat(absl::StripPrefix(root_path, "/"), "/resolutions/", bin_size),
      throw_on_failure, false);
}

}  // namespace modle::cooler
