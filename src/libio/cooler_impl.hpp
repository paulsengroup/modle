#pragma once

#include <H5Cpp.h>
#include <absl/strings/match.h>
#include <absl/types/span.h>
#include <fmt/format.h>

#include <algorithm>
#include <cstdint>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include "modle/contacts.hpp"
#include "modle/hdf5.hpp"

namespace modle::cooler {

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

H5::StrType Cooler::generate_default_str_type(std::size_t max_str_length) {
  // Cooltools doesn't seem to properly handle variable length strings (H5T_VARIABLE)
  // For the time being we are forced to use fixed length, null-padded strings
  auto st = H5::StrType(H5::PredType::C_S1, max_str_length);
  st.setStrpad(H5T_STR_NULLPAD);
  st.setCset(H5T_CSET_ASCII);
  return st;
}

template <typename T1, typename T2>
std::unique_ptr<H5::DSetCreatPropList> Cooler::generate_default_cprop(hsize_t chunk_size,
                                                                      uint8_t compression_lvl,
                                                                      T1 type, T2 fill_value) {
  static_assert(
      std::is_same_v<T2, int64_t> || std::is_same_v<T2, int32_t> ||
          std::is_constructible_v<H5std_string, T2>,
      "fill_value should have one of the following types: int64_t, int32_t or std::string.");
  static_assert(
      (std::is_constructible_v<H5std_string, T2> && std::is_same_v<T1, H5::StrType>) ||
          std::is_same_v<T1, H5::PredType>,
      "Incompatible data type for variables type and fill_value: if T2 is string "
      "constructible, then T1 must be H5::StrType, else T2 is integral and type is H5::PredType");

  H5::DSetCreatPropList prop{};
  prop.setChunk(1, &chunk_size);
  prop.setDeflate(compression_lvl);
  prop.setFillValue(type, &fill_value);

  return std::make_unique<H5::DSetCreatPropList>(prop);
}

template <typename T>
std::unique_ptr<H5::DSetAccPropList> Cooler::generate_default_aprop(T type, hsize_t chunk_size,
                                                                    hsize_t cache_size) {
  static_assert(std::is_same_v<T, H5::StrType> || std::is_same_v<T, H5::PredType>,
                "type should be of type H5::StrType or H5::PredType");
  H5::DSetAccPropList prop{};
  // https://support.hdfgroup.org/HDF5/doc/Advanced/Chunking/index.html
  constexpr std::size_t default_multiplier{100};
  constexpr double rdcc_w0{0.99};
  if constexpr (std::is_same_v<T, H5::StrType>) {
    prop.setChunkCache(default_multiplier * (cache_size / chunk_size), cache_size, rdcc_w0);
  } else {
    if (type == H5::PredType::NATIVE_INT64) {  // int64_t
      prop.setChunkCache(default_multiplier * (cache_size / chunk_size), cache_size, rdcc_w0);
    } else if (type == H5::PredType::NATIVE_INT) {  // int32_t
      prop.setChunkCache(default_multiplier * (cache_size / (chunk_size + chunk_size)), cache_size,
                         rdcc_w0);
    } else {
      throw std::runtime_error(
          "Cooler::generate_default_aprop(), type should have type H5::StrType or be one of "
          "H5::PredType:NATIVE_INT, H5::PredType::NATIVE_INT64");
    }
  }
  return std::make_unique<H5::DSetAccPropList>(prop);
}

std::unique_ptr<H5::H5File> Cooler::open_file(const std::filesystem::path &path, IO_MODE mode,
                                              std::size_t bin_size, Flavor flavor, bool validate) {
  if (mode == WRITE_ONLY && bin_size == 0) {
    throw std::runtime_error(
        "Cooler::open_file(): bin_size cannot be 0 when file is being opened in WRITE_ONLY mode");
  }
  try {
    H5::H5File f(path.c_str(), mode == READ_ONLY ? H5F_ACC_RDONLY : H5F_ACC_TRUNC);
    if (validate) {
      (void)validate_file_format(f, flavor, mode, bin_size, true);
    }
    return std::make_unique<H5::H5File>(f);
  } catch (const H5::Exception &e) {
    throw std::runtime_error(e.getDetailMsg());
  }
}

Cooler::Flavor Cooler::detect_file_flavor(H5::H5File &f) {
  Cooler::Flavor flavor{UNK};

  if (const auto &fname = f.getFileName(); absl::EndsWithIgnoreCase(fname, ".cool")) {
    flavor = COOL;
  } else if (absl::EndsWithIgnoreCase(fname, ".mcool")) {
    flavor = MCOOL;
  } else if (absl::EndsWithIgnoreCase(fname, ".scool")) {
    flavor = SCOOL;
  }
  if (flavor != UNK) {
    return flavor;
  }

  std::string buff;
  if (hdf5::has_attribute(f, "format")) {
    hdf5::read_attribute(f, "format", buff);
  }
  if (const auto &fname = f.getFileName(); absl::EndsWithIgnoreCase(fname, "::cooler")) {
    flavor = COOL;
  } else if (absl::EndsWithIgnoreCase(fname, "::mcool")) {
    flavor = MCOOL;
  } else if (absl::EndsWithIgnoreCase(fname, "::scool")) {
    flavor = SCOOL;
  }
  if (flavor != UNK) {
    return flavor;
  }

  if (f.nameExists("/resolutions")) {
    flavor = MCOOL;
  } else if (f.nameExists("/bins") && f.nameExists("/chroms")) {
    if (f.nameExists("/pixels") && f.nameExists("/indexes")) {
      flavor = COOL;
    } else if (f.nameExists("/cells")) {
      flavor = SCOOL;
    }
  }

  if (flavor == UNK) {
    throw std::runtime_error("Unable to detect Cooler file flavor");
  }

  return flavor;
}

bool Cooler::validate_file_format(H5::H5File &f, Flavor expected_flavor, IO_MODE mode,
                                  std::size_t bin_size, bool throw_on_failure) {
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

bool Cooler::validate_cool_flavor(H5::H5File &f, std::size_t bin_size, std::string_view root_path,
                                  bool throw_on_failure) {
  /* The following attributes are being checked:
   * format
   * format-version
   * bin-type
   * bin size
   * storage-mode
   */

  std::string str_buff;
  int64_t int_buff;

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

    hdf5::read_attribute(f, "format-version", int_buff, root_path);
    if (int_buff < min_format_ver || int_buff > max_format_ver) {
      if (!throw_on_failure) {
        return false;
      }
      throw std::runtime_error(fmt::format(
          FMT_STRING("Expected format-version attribute to be between {} and {}, got {}"),
          min_format_ver, max_format_ver, int_buff));
    }
    const auto format_version = static_cast<uint8_t>(int_buff);

    hdf5::read_attribute(f, "bin-type", str_buff, root_path);
    if (str_buff != "fixed") {
      if (!throw_on_failure) {
        return false;
      }
      throw std::runtime_error(
          fmt::format(FMT_STRING("Expected bin-type attribute to be 'fixed', got '{}'"), str_buff));
    }

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

    // Checking groups
    auto check_group = [&](std::string_view name) {
      if (!hdf5::group_exists(f, name, root_path)) {
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
      if (!hdf5::dataset_exists(f, "chroms/length", root_path)) {
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

std::vector<H5::Group> Cooler::open_groups(H5::H5File &f, bool create_if_not_exist) {
  std::vector<H5::Group> groups(4);  // NOLINT
  auto open_or_create_group = [&](Cooler::Groups g, const std::string &group_name) {
    groups[g] = create_if_not_exist && !f.nameExists(group_name) ? f.createGroup(group_name)
                                                                 : f.openGroup(group_name);
  };
  try {
    open_or_create_group(CHR, "chroms");
    open_or_create_group(BIN, "bins");
    open_or_create_group(PXL, "pixels");
    open_or_create_group(IDX, "indexes");

    return groups;
  } catch (const H5::Exception &e) {
    throw std::runtime_error(e.getDetailMsg());
  }
}

bool Cooler::validate_multires_cool_flavor(H5::H5File &f, std::size_t bin_size,
                                           std::string_view root_path, bool throw_on_failure) {
  constexpr int64_t min_format_ver = 2;
  constexpr int64_t max_format_ver = 3;

  if (bin_size == 0) {
    throw std::runtime_error(
        "A bin size other than 0 is required when calling Cooler::validate_multires_cool_flavor()");
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

  if (!hdf5::group_exists(f, "resolutions", root_path)) {
    if (!throw_on_failure) {
      return false;
    }
    throw std::runtime_error(fmt::format(FMT_STRING("Missing mandatory group '{}/resolutions'"),
                                         absl::StripPrefix(root_path, "/")));
  }

  if (!hdf5::group_exists(f, absl::StrCat("resolutions/", bin_size), root_path)) {
    if (!throw_on_failure) {
      return false;
    }
    throw std::runtime_error(
        fmt::format(FMT_STRING("Missing group for resolution {} bp"), bin_size));
  }

  return validate_cool_flavor(
      f, bin_size, absl::StrCat(absl::StripPrefix(root_path, "/"), "/resolutions/", bin_size),
      throw_on_failure);
}

Cooler::Cooler(std::string_view path_to_file, IO_MODE mode, std::size_t bin_size,
               std::size_t max_str_length, std::string_view assembly_name, Flavor flavor,
               bool validate, uint8_t compression_lvl, std::size_t chunk_size,
               std::size_t cache_size)
    : STR_TYPE(generate_default_str_type(max_str_length)),
      _path_to_file(std::string{path_to_file.data(), path_to_file.size()}),
      _mode(mode),
      _bin_size(bin_size),
      _assembly_name(assembly_name.data(), assembly_name.size()),
      _flavor(flavor),
      _fp(open_file(_path_to_file, _mode, _bin_size, _flavor, validate)),
      _groups(open_groups(*_fp, !this->is_read_only())),
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
      _aprop_int64(generate_default_aprop(Cooler::INT64_TYPE, _chunk_size, _cache_size)) {
  assert(this->_flavor != UNK);
  if (this->_mode == READ_ONLY && this->_flavor == AUTO) {
    this->_flavor = Cooler::detect_file_flavor(*this->_fp);
  }
}

bool Cooler::is_read_only() const { return this->_mode == Cooler::READ_ONLY; }

void Cooler::write_metadata() {
  if (this->is_read_only()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Caught attempt to write metadata to an HDF5 file that is opened in "
                               "read-only mode. File name: {}"),
                    this->_path_to_file.string()));
  }
  assert(this->_bin_size != 0);
  H5::StrType METADATA_STR_TYPE(H5::PredType::C_S1, H5T_VARIABLE);
  METADATA_STR_TYPE.setCset(H5T_CSET_UTF8);

  H5::DataSpace attr_space(H5S_SCALAR);
  int32_t int_buff{};
  std::string str_buff{};

  auto att = this->_fp->createAttribute("format", METADATA_STR_TYPE, attr_space);
  str_buff = "HDF5::Cooler";
  att.write(METADATA_STR_TYPE, str_buff);

  att = this->_fp->createAttribute("format-version", H5::PredType::NATIVE_INT64, attr_space);
  int_buff = 3;
  att.write(H5::PredType::NATIVE_INT, &int_buff);

  att = this->_fp->createAttribute("bin-type", METADATA_STR_TYPE, attr_space);
  str_buff = "fixed";
  att.write(METADATA_STR_TYPE, str_buff);

  att = this->_fp->createAttribute("bin-size", H5::PredType::NATIVE_INT64, attr_space);
  att.write(H5::PredType::NATIVE_INT, &this->_bin_size);

  att = this->_fp->createAttribute("storage-mode", METADATA_STR_TYPE, attr_space);
  str_buff = "symmetric-upper";
  att.write(METADATA_STR_TYPE, str_buff);

  if (!this->_assembly_name.empty()) {
    str_buff = this->_assembly_name;
    att = this->_fp->createAttribute("assembly-name", METADATA_STR_TYPE, attr_space);
    att.write(METADATA_STR_TYPE, &str_buff);
  }

  att = this->_fp->createAttribute("generated-by", METADATA_STR_TYPE, attr_space);
  str_buff = fmt::format(FMT_STRING("ModLE-v{}.{}.{}"), 0, 0, 1);  // TODO make ModLE ver a tunable
  att.write(METADATA_STR_TYPE, str_buff);

  att = this->_fp->createAttribute("creation-date", METADATA_STR_TYPE, attr_space);
  str_buff = absl::FormatTime(absl::Now(), absl::UTCTimeZone());
  att.write(METADATA_STR_TYPE, str_buff);
}

void Cooler::init_default_datasets() {
  constexpr hsize_t RANK{1};                 // i.e. number of dimensions
  constexpr hsize_t MAXDIMS{H5S_UNLIMITED};  // extensible dataset
  constexpr hsize_t BUFF_SIZE{1};            // Dummy buffer size

  this->_datasets.clear();
  this->_datasets.reserve(10);  // NOLINT

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

  auto &dset = this->_datasets;

  try {
    auto mem_space{H5::DataSpace(RANK, &BUFF_SIZE, &MAXDIMS)};
    // Do not change the order of these pushbacks
    dset.push_back(this->_fp->createDataSet("chroms/length", INT64, mem_space, cp64, ap64));
    dset.push_back(this->_fp->createDataSet("chroms/name", STR, mem_space, cps, aps));

    dset.push_back(this->_fp->createDataSet("bins/chrom", INT64, mem_space, cp64, ap64));
    dset.push_back(this->_fp->createDataSet("bins/start", INT64, mem_space, cp64, ap64));
    dset.push_back(this->_fp->createDataSet("bins/end", INT64, mem_space, cp64, ap64));

    dset.push_back(this->_fp->createDataSet("pixels/bin1_id", INT64, mem_space, cp64, ap64));
    dset.push_back(this->_fp->createDataSet("pixels/bin2_id", INT64, mem_space, cp64, ap64));
    dset.push_back(this->_fp->createDataSet("pixels/count", INT32, mem_space, cp32, ap32));

    dset.push_back(this->_fp->createDataSet("indexes/bin1_offset", INT64, mem_space, cp64, ap64));
    dset.push_back(this->_fp->createDataSet("indexes/chrom_offset", INT64, mem_space, cp64, ap64));

  } catch (const H5::FileIException &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("An error occurred while initializing default Cooler dataset on file '{}': {}"),
        this->_path_to_file, hdf5::construct_error_stack()));
  }
}

template <typename I>
void Cooler::write_cmatrix_to_file(const ContactMatrix<I> &cmatrix, std::string_view chr_name,
                                   uint64_t chr_start, uint64_t chr_end, uint64_t chr_length) {
  // Casting away the constness from the ptr is needed in order to make things work with slices.
  // This shouldn't cause any problem, as in the end the contact matrix is accesses through a const
  // slice NOLINTNEXTLINE
  const std::vector<ContactMatrix<I> *> cmatrices{const_cast<ContactMatrix<I> *>(&cmatrix)};
  std::string chr_name_{chr_name.data(), chr_name.size()};
  write_cmatrix_to_file(absl::MakeConstSpan(cmatrices), absl::MakeConstSpan(&chr_name_, 1),
                        absl::MakeConstSpan(&chr_start, 1), absl::MakeConstSpan(&chr_end, 1),
                        absl::MakeConstSpan(&chr_length, 1));
}

template <typename I>
void Cooler::write_cmatrix_to_file(const std::vector<ContactMatrix<I>> &cmatrices,
                                   const std::vector<std::string> &chr_names,
                                   const std::vector<uint64_t> &chr_starts,
                                   const std::vector<uint64_t> &chr_ends,
                                   const std::vector<uint64_t> &chr_sizes) {
  std::vector<ContactMatrix<I> *> v(cmatrices.size());
  std::transform(cmatrices.begin(), cmatrices.end(), v.begin(), [](const auto &m) { return &m; });
  write_cmatrix_to_file(v, chr_names, chr_starts, chr_ends, chr_sizes);
}

template <typename I>
void Cooler::write_cmatrix_to_file(const std::vector<ContactMatrix<I> *> &cmatrices,
                                   const std::vector<std::string> &chr_names,
                                   const std::vector<uint64_t> &chr_starts,
                                   const std::vector<uint64_t> &chr_ends,
                                   const std::vector<uint64_t> &chr_sizes) {
  Cooler::write_cmatrix_to_file(absl::MakeConstSpan(cmatrices), absl::MakeConstSpan(chr_names),
                                absl::MakeConstSpan(chr_starts), absl::MakeConstSpan(chr_ends),
                                absl::MakeConstSpan(chr_sizes));
}

template <typename I>
void Cooler::write_cmatrix_to_file(absl::Span<ContactMatrix<I> *const> cmatrices,
                                   absl::Span<const std::string> chr_names,
                                   absl::Span<const uint64_t> chr_starts,
                                   absl::Span<const uint64_t> chr_ends,
                                   absl::Span<const uint64_t> chr_sizes) {
  static_assert(std::is_integral_v<I>, "I should be an integral type.");

  assert(this->_bin_size != 0);

  const auto n_chromosomes = cmatrices.size();
  if (n_chromosomes != chr_names.size() || n_chromosomes != chr_starts.size() ||
      n_chromosomes != chr_ends.size() || n_chromosomes != chr_sizes.size()) {
    utils::throw_with_trace(std::runtime_error(fmt::format(
        FMT_STRING("Message for the DEVS: The vectors passed to function "
                   "cooler::write_modle_cmatrices_to_cooler() should all have "
                   "the same size!\n  cmatrices.size()={}\n  chr_names.size()={}\n  "
                   "chr_starts={}\n  chr_ends={}\n  chr_sizes={}\n"),
        cmatrices.size(), chr_names.size(), chr_starts.size(), chr_ends.size(), chr_sizes.size())));
  }

  auto t0 = absl::Now();

  this->init_default_datasets();

  std::size_t chr_name_foffset = 0;
  std::size_t chr_length_foffset = 0;

  constexpr std::size_t BUFF_SIZE = 1024 * 1024 / sizeof(int64_t);  // 1 MB
  std::vector<int64_t> bin_pos_buff;
  std::vector<int32_t> bin_chrom_buff;
  std::vector<int64_t> pixel_b1_idx_buff;
  std::vector<int64_t> pixel_b2_idx_buff;
  std::vector<int64_t> pixel_count_buff;
  std::vector<int64_t> idx_bin1_offset_buff;
  std::vector<int64_t> idx_chrom_offset_buff;

  bin_pos_buff.reserve(BUFF_SIZE);
  bin_chrom_buff.reserve(BUFF_SIZE);
  pixel_b1_idx_buff.reserve(BUFF_SIZE);
  pixel_b2_idx_buff.reserve(BUFF_SIZE);
  pixel_count_buff.reserve(BUFF_SIZE);
  idx_bin1_offset_buff.reserve(BUFF_SIZE);
  idx_chrom_offset_buff.reserve(cmatrices.size());

  hsize_t pixel_b1_idx_h5df_offset = 0;
  hsize_t pixel_b2_idx_h5df_offset = 0;
  hsize_t pixel_count_h5df_offset = 0;
  hsize_t idx_bin1_offset_h5df_offset = 0;

  auto &d = this->_datasets;

  auto write_pixels_to_file = [&]() {
    pixel_b1_idx_h5df_offset =
        hdf5::write_numbers(pixel_b1_idx_buff, d[PXL_B1], pixel_b1_idx_h5df_offset);
    pixel_b2_idx_h5df_offset =
        hdf5::write_numbers(pixel_b2_idx_buff, d[PXL_B2], pixel_b2_idx_h5df_offset);
    pixel_count_h5df_offset =
        hdf5::write_numbers(pixel_count_buff, d[PXL_COUNT], pixel_count_h5df_offset);

    pixel_b1_idx_buff.clear();
    pixel_b2_idx_buff.clear();
    pixel_count_buff.clear();
  };

  int64_t nnz = 0;
  int64_t offset = 0;
  int64_t nbins = 0;

  this->write_metadata();

  for (auto chr_idx = 0UL; chr_idx < n_chromosomes; ++chr_idx) {
    const auto &cmatrix = cmatrices[chr_idx];
    const auto &chr_name = chr_names[chr_idx];
    const auto &chr_total_len = chr_sizes[chr_idx];
    const auto &chr_start = chr_starts[chr_idx];
    const auto &chr_end = chr_ends[chr_idx];
    const auto chr_simulated_len = chr_end - chr_start;
    chr_name_foffset = hdf5::write_str(chr_name, d[CHR_NAME], STR_TYPE, chr_name_foffset);
    chr_length_foffset = hdf5::write_number(chr_total_len, d[CHR_LEN], chr_length_foffset);
    DISABLE_WARNING_PUSH
    DISABLE_WARNING_CONVERSION
    DISABLE_WARNING_SIGN_CONVERSION
    fmt::print(stderr, FMT_STRING("Writing contacts for '{}' ({:.2f} Mbp)..."), chr_name,
               chr_simulated_len / 1.0e6);  // NOLINT
    const auto t1 = absl::Now();
    idx_chrom_offset_buff.push_back(nbins);
    nbins = this->write_bins(chr_idx, chr_simulated_len, this->_bin_size,  // NOLINT
                             bin_chrom_buff, bin_pos_buff, nbins);
    DISABLE_WARNING_POP

    for (auto i = 0UL; i < chr_start / this->_bin_size; ++i) {
      idx_bin1_offset_buff.push_back(nnz);
      if (idx_bin1_offset_buff.size() == BUFF_SIZE) {
        idx_bin1_offset_h5df_offset =
            hdf5::write_numbers(idx_bin1_offset_buff, d[IDX_BIN1], idx_bin1_offset_h5df_offset);
        idx_bin1_offset_buff.clear();
      }
    }

    offset += static_cast<int64_t>(chr_start / this->_bin_size);
    for (auto i = 0UL; i < cmatrix->n_cols(); ++i) {
      idx_bin1_offset_buff.push_back(nnz);
      for (auto j = i; j < i + cmatrix->n_rows() && j < cmatrix->n_cols(); ++j) {
        if (const auto m = cmatrix->get(i, j); m != 0) {
          pixel_b1_idx_buff.push_back(offset + static_cast<int64_t>(i));
          pixel_b2_idx_buff.push_back(offset + static_cast<int64_t>(j));
#ifndef NDEBUG
          if (offset + static_cast<int64_t>(i) > offset + static_cast<int64_t>(j)) {
            utils::throw_with_trace(std::runtime_error(
                fmt::format(FMT_STRING("b1 > b2: b1={}; b2={}; offset={}; m={}\n"),
                            pixel_b1_idx_buff.back(), pixel_b2_idx_buff.back(), offset, m)));
          }
#endif
          pixel_count_buff.push_back(m);
          ++nnz;
          if (pixel_b1_idx_buff.size() == BUFF_SIZE) {
            write_pixels_to_file();
          }
        }
      }
      if (idx_bin1_offset_buff.size() == BUFF_SIZE) {
        idx_bin1_offset_h5df_offset =
            hdf5::write_numbers(idx_bin1_offset_buff, d[IDX_BIN1], idx_bin1_offset_h5df_offset);
        idx_bin1_offset_buff.clear();
      }
    }
    for (auto i = 0UL; i < (chr_total_len - chr_end) / this->_bin_size; ++i) {
      idx_bin1_offset_buff.push_back(nnz);
      if (idx_bin1_offset_buff.size() == BUFF_SIZE) {
        idx_bin1_offset_h5df_offset =
            hdf5::write_numbers(idx_bin1_offset_buff, d[IDX_BIN1], idx_bin1_offset_h5df_offset);
        idx_bin1_offset_buff.clear();
      }
    }

    offset = nbins;
    fmt::print(stderr, FMT_STRING(" DONE in {}!\n"), absl::FormatDuration(absl::Now() - t1));
  }

  if (!pixel_b1_idx_buff.empty()) {
    write_pixels_to_file();
  }

  idx_chrom_offset_buff.push_back(nbins);
  idx_bin1_offset_buff.push_back(nnz);

  (void)hdf5::write_numbers(idx_chrom_offset_buff, d[IDX_CHR], 0);
  (void)hdf5::write_numbers(idx_bin1_offset_buff, d[IDX_BIN1], idx_bin1_offset_h5df_offset);

  H5::DataSpace att_space(H5S_SCALAR);
  auto att = this->_fp->createAttribute("nchroms", INT64_TYPE, att_space);
  att.write(INT64_TYPE, &n_chromosomes);
  att = this->_fp->createAttribute("nbins", INT64_TYPE, att_space);
  att.write(INT64_TYPE, &nbins);
  att = this->_fp->createAttribute("nnz", INT64_TYPE, att_space);
  att.write(INT64_TYPE, &nnz);

  fmt::print(stderr, "DONE! Saved {} pixels in {}.\n", nnz, absl::FormatDuration(absl::Now() - t0));
}

hsize_t Cooler::write_bins(int32_t chrom, int64_t length, int64_t bin_size,
                           std::vector<int32_t> &buff32, std::vector<int64_t> &buff64,
                           hsize_t file_offset, hsize_t buff_size) {
  assert(bin_size != 0);
  buff32.resize(buff_size);
  buff64.resize(buff_size);

  std::fill(buff32.begin(), buff32.end(), chrom);

  int64_t start = 0;
  int64_t end = std::min(length, bin_size);
  hsize_t bins_processed = 0;
  const int64_t nbins = (length / bin_size) + 1;
  const int64_t nchunks = (nbins / static_cast<int64_t>(buff_size)) + 1;
  auto &d = this->_datasets;

  for (auto i = 0; i < nchunks; ++i) {
    const auto chunk_size = std::min(buff_size, static_cast<hsize_t>(nbins) - bins_processed);
    if (chunk_size != buff_size) {
      buff32.resize(chunk_size);
      buff64.resize(chunk_size);
    }
    (void)hdf5::write_numbers(buff32, d[BIN_CHROM], file_offset);

    std::generate(buff64.begin(), buff64.end(), [&]() { return (start += bin_size) - bin_size; });
    (void)hdf5::write_numbers(buff64, d[BIN_START], file_offset);

    std::generate(buff64.begin(), buff64.end(), [&]() { return (end += bin_size) - bin_size; });
    if (chunk_size != buff_size) {
      buff64.back() = length;
    }
    file_offset = hdf5::write_numbers(buff64, d[BIN_END], file_offset);
    bins_processed += chunk_size;
  }
  assert(buff64.back() == length);

  return file_offset;
}

void Cooler::open_cooler_datasets() {
  assert(this->_fp);
  assert(this->_aprop_int32);
  assert(this->_aprop_int64);
  assert(this->_aprop_str);

  auto &d = this->_datasets;
  const auto &a32 = *this->_aprop_int32;
  const auto &a64 = *this->_aprop_int64;
  const auto &as = *this->_aprop_str;

  this->_datasets.resize(10);

  try {
    // Do not change the order of these pushbacks
    d[CHR_LEN] = this->_fp->openDataSet("chroms/length", a64);
    d[CHR_NAME] = this->_fp->openDataSet("chroms/name", as);

    d[BIN_CHROM] = this->_fp->openDataSet("bins/chrom", a64);
    d[BIN_START] = this->_fp->openDataSet("bins/start", a64);
    d[BIN_END] = this->_fp->openDataSet("bins/end", a64);

    d[PXL_B1] = this->_fp->openDataSet("pixels/bin1_id", a64);
    d[PXL_B2] = this->_fp->openDataSet("pixels/bin2_id", a64);
    d[PXL_COUNT] = this->_fp->openDataSet("pixels/count", a32);

    d[IDX_BIN1] = this->_fp->openDataSet("indexes/bin1_offset", a64);
    d[IDX_CHR] = this->_fp->openDataSet("indexes/chrom_offset", a64);

  } catch (const H5::FileIException &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING(
            "An error occurred while trying to open Cooler's default dataset of file '{}': {}"),
        this->_path_to_file, hdf5::construct_error_stack()));
  }
}

ContactMatrix<uint32_t> Cooler::cooler_to_cmatrix(int64_t bin_offset,
                                                  const std::vector<int64_t> &bin1_offset_idx,
                                                  std::size_t nrows) {
  return this->cooler_to_cmatrix(bin_offset, absl::MakeConstSpan(bin1_offset_idx), nrows);
}

ContactMatrix<uint32_t> Cooler::cooler_to_cmatrix(std::string_view chr_name,
                                                  std::size_t diagonal_width, std::size_t bin_size,
                                                  bool try_common_chr_prefixes) {
  assert(this->_bin_size != 0);
  if (bin_size != 0 && this->_bin_size != bin_size) {
    throw std::runtime_error(fmt::format(
        FMT_STRING(
            "Unable to read a Cooler file with bin size {} in a contact matrix of bin size {}"),
        this->_bin_size, bin_size));
  }
  const auto nrows = (diagonal_width / this->_bin_size) +
                     static_cast<std::size_t>(diagonal_width % this->_bin_size != 0);
  return cooler_to_cmatrix(chr_name, nrows, try_common_chr_prefixes);
}

ContactMatrix<uint32_t> Cooler::cooler_to_cmatrix(std::string_view chr_name, std::size_t nrows,
                                                  bool try_common_chr_prefixes) {
  // TODO: Add check for Cooler bin size

  {
    std::scoped_lock lock(this->_mutex);
    if (this->_datasets.empty()) {
      this->open_cooler_datasets();
    }
    if (this->_idx_chrom_offset.empty()) {
      this->read_chr_offset_idx();
    }
    if (this->_idx_bin1_offset.empty()) {
      this->read_bin1_offset_idx();
    }
  }

  auto search_chr_offset_idx = [&]() {
    if (!try_common_chr_prefixes) {
      return std::make_pair(chr_name, this->get_chr_idx(chr_name));
    }
    // Here's the issue: for a given genome assembly (say hg19), some tools name chromosome as
    // chrN, where N is the chromosome number, while others only use the number N. The purpose
    // of this lambda is to try to guess few reasonably common prefixes, look them up in the
    // Cooler file, then return the chromosome name variant that produced a hit together with
    // the chr index
    std::vector<std::string_view> chr_names{chr_name};
    for (const auto &prefix : {"chr", "CHR", "Chr"}) {
      if (absl::StartsWith(chr_name, prefix)) {
        chr_names.emplace_back(absl::StripPrefix(chr_name, prefix));
      } else {
        chr_names.emplace_back(absl::StrCat(prefix, chr_name));
      }
    }
    // TODO: Consider reading the full array of chr names and then look up all the combinations.
    // This would allows us to avoid using exceptions, but will increase memory usage when
    // processing highly fragmented assemblies
    chr_name = "";  // Just to make sure we are not using this variable instead of actual_chr_name
    for (const auto name : chr_names) {
      try {
        return std::make_pair(name, this->get_chr_idx(name));
      } catch (const std::runtime_error &e) {
        if (!absl::StartsWith(e.what(), "Unable to find a chromosome named")) {
          throw;
        }
      }
    }
    throw std::runtime_error(
        fmt::format(FMT_STRING("Unable to find a chromosome named '{}'. The following chromosome "
                               "name variants were searched: '{}'"),
                    chr_name, absl::StrJoin(chr_names, "', '")));
  };

  const auto [actual_chr_name, chr_idx] = search_chr_offset_idx();
  const auto bin1_offset_idx = this->get_bin1_offset_idx_for_chr(chr_idx);
  const auto bin_offset = this->_idx_chrom_offset[chr_idx];

  return cooler_to_cmatrix(bin_offset, bin1_offset_idx, nrows);
}

ContactMatrix<uint32_t> Cooler::cooler_to_cmatrix(int64_t bin_offset,
                                                  absl::Span<const int64_t> bin1_offset_idx,
                                                  std::size_t nrows) {
  ContactMatrix<uint32_t> cmatrix(nrows, bin1_offset_idx.size());
  std::vector<int64_t> bin1_BUFF(nrows);
  std::vector<int64_t> bin2_BUFF(nrows);
  std::vector<int64_t> count_BUFF(nrows);
  const auto &d = this->_datasets;

  for (auto i = 1UL; i < bin1_offset_idx.size(); ++i) {
    const auto file_offset = static_cast<hsize_t>(bin1_offset_idx[i - 1]);
    const auto buff_size =  // Figure out why buff size is always small
        std::min(static_cast<std::size_t>(bin1_offset_idx[i] - bin1_offset_idx[i - 1]), nrows);
    if (buff_size == 0) {
      continue;
    }
    // fmt::print(stderr, "i={}; buff_size={}\n", i, buff_size);

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
      if (const auto bin2_id = bin2_BUFF[j] - bin_offset;
          bin2_id >= i + nrows - 1 || bin2_id >= bin1_offset_idx.size()) {
        break;
      }
      // fmt::print(stderr, "m[{}][{}]={} (nrows={}; ncols={})\n", bin2_BUFF[j] - bin_offset,
      //           bin1_BUFF[j] - bin_offset, count_BUFF[j], cmatrix.n_rows(),
      //           cmatrix.n_cols());
      cmatrix.set(bin2_BUFF[j] - bin_offset, bin1_BUFF[j] - bin_offset, count_BUFF[j]);

      DISABLE_WARNING_POP
    }
  }
  return cmatrix;
}

std::size_t Cooler::get_chr_idx(std::string_view chr_name) {
  try {
    const auto chr_names = hdf5::read_strings(this->_datasets[CHR_NAME], 0);
    const auto match = std::find(chr_names.begin(), chr_names.end(), chr_name);
    if (match != chr_names.end()) {
      return static_cast<std::size_t>(std::distance(chr_names.begin(), match));
    }
    throw std::runtime_error(
        fmt::format(FMT_STRING("Unable to find a chromosome named '{}'"), chr_name));
  } catch (const std::exception &e) {
    if (absl::StartsWith(e.what(), "Unable to find a chromosome")) {
      throw;
    }
    throw std::runtime_error(fmt::format(FMT_STRING("The following error occurred while looking up "
                                                    "'{}' in dataset chroms/name of file '{}': {}"),
                                         chr_name, this->_path_to_file, e.what()));
  }
}

std::size_t Cooler::read_chr_offset_idx() {
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

std::size_t Cooler::read_bin1_offset_idx() {
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

absl::Span<const int64_t> Cooler::get_bin1_offset_idx_for_chr(std::string_view chr_name) {
  const auto chr_idx = get_chr_idx(chr_name);
  return get_bin1_offset_idx_for_chr(chr_idx);
}

absl::Span<const int64_t> Cooler::get_bin1_offset_idx_for_chr(std::size_t chr_idx) {
  assert(!this->_idx_bin1_offset.empty());           // NOLINT
  assert(!this->_idx_chrom_offset.empty());          // NOLINT
  assert(chr_idx < this->_idx_chrom_offset.size());  // NOLINT
  const auto chr_start_bin = static_cast<std::size_t>(this->_idx_chrom_offset[chr_idx]);
  const auto chr_end_bin = static_cast<std::size_t>(this->_idx_chrom_offset[chr_idx + 1]);
  assert(chr_end_bin >= chr_start_bin);  // NOLINT

  return absl::MakeConstSpan(this->_idx_bin1_offset)
      .subspan(chr_start_bin, chr_end_bin - chr_start_bin);
}

std::pair<int64_t, int64_t> Cooler::read_chrom_pixels_boundaries(std::string_view chr_name) {
  const auto chr_idx = get_chr_idx(chr_name);
  return read_chrom_pixels_boundaries(chr_idx);
}

std::pair<int64_t, int64_t> Cooler::read_chrom_pixels_boundaries(std::size_t chr_idx) {
  assert(chr_idx < this->_idx_chrom_offset.size());  // NOLINT
  assert(!this->_idx_chrom_offset.empty());          // NOLINT
  const auto chr_start_bin = static_cast<std::size_t>(this->_idx_chrom_offset[chr_idx]);
  const auto chr_end_bin = static_cast<std::size_t>(this->_idx_chrom_offset[chr_idx + 1]);
  assert(chr_end_bin >= chr_start_bin);  // NOLINT

  std::pair<int64_t, int64_t> pixel_boundaries{};
  const auto &d = this->_datasets[IDX_BIN1];
  (void)hdf5::read_number(d, pixel_boundaries.first, chr_start_bin);
  (void)hdf5::read_number(d, pixel_boundaries.second, chr_end_bin);

  return pixel_boundaries;
}

bool Cooler::is_cool() const { return this->_flavor == COOL; }
bool Cooler::is_mcool() const { return this->_flavor == MCOOL; }
bool Cooler::is_scool() const { return this->_flavor == SCOOL; }

std::size_t Cooler::get_n_chroms() {
  assert(this->_fp);
  if (this->is_cool()) {
    return static_cast<std::size_t>(hdf5::read_attribute_int(*this->_fp, "nchroms"));
  }
  if (this->is_mcool()) {
    assert(this->_bin_size != 0);
    return static_cast<std::size_t>(hdf5::read_attribute_int(
        *this->_fp, "nchroms", absl::StrCat("/resolutions", this->_bin_size)));
  }
  throw std::runtime_error(
      "Message for the devs: Cooler::get_n_chroms() executed code that should be unreachable!");
}

void Cooler::get_chr_names(std::vector<std::string> &buff) {
  assert(this->_fp);
  const auto nchroms = this->get_n_chroms();
  buff.resize(nchroms);
  if (!this->_datasets.empty()) {
    const auto &d = this->_datasets[CHR_NAME];
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
std::vector<std::string> Cooler::get_chr_names() {
  std::vector<std::string> buff;
  this->get_chr_names(buff);
  return buff;
}

template <typename I>
void Cooler::get_chr_sizes(std::vector<I> &buff) {
  static_assert(std::is_integral_v<I>, "buff should be a vector of integers.");
  assert(this->_fp);
  const auto nchroms = this->get_n_chroms();
  buff.resize(nchroms);
  if (!this->_datasets.empty()) {
    const auto &d = this->_datasets[CHR_LEN];
    (void)hdf5::read_numbers(d, buff, 0);
  } else {
    if (this->is_cool()) {
      auto d = this->_fp->openDataSet("/chroms/length", *this->_aprop_str);
      (void)hdf5::read_numbers(d, buff, 0);
    } else if (this->is_mcool()) {
      assert(this->_bin_size != 0);
      auto d = this->_fp->openDataSet(
          absl::StrCat("/resolutions/", this->_bin_size, "/chroms/length"), *this->_aprop_str);
      (void)hdf5::read_numbers(d, buff, 0);
    } else {
      assert(!this->is_scool());
      buff.clear();
    }
  }
}

std::vector<int64_t> Cooler::get_chr_sizes() {
  std::vector<int64_t> buff;
  this->get_chr_sizes(buff);
  return buff;
}

}  // namespace modle::cooler
