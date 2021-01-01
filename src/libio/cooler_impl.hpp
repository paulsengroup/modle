#pragma once

#include <H5Cpp.h>
#include <absl/strings/match.h>
#include <fmt/format.h>

#include <algorithm>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <type_traits>

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
}

H5::StrType Cooler::generate_default_str_type() {
  auto st = H5::StrType(H5::PredType::C_S1);
  st.setStrpad(H5T_STR_NULLTERM);
  st.setCset(H5T_CSET_ASCII);
  return st;
}

template <typename T1, typename T2>
std::unique_ptr<H5::DSetCreatPropList> Cooler::generate_default_cprop(hsize_t chunk_size,
                                                                      uint8_t compression_lvl,
                                                                      T1 type, T2 fill_value) {
  static_assert(
      std::is_same_v<T2, int64_t> || std::is_same_v<T2, int32_t> || std::is_same_v<T2, char>,
      "fill_value should have one of the following types: int64_t, int32_t char.");
  static_assert((std::is_same_v<T2, char> && std::is_same_v<T1, H5::StrType>) ||
                    std::is_same_v<T1, H5::PredType>,
                "Incompatible data type for variables type and fill_value: if T2 is char, then T1 "
                "must be H5::StrType, else T2 is integral and type is H5::PredType");

  H5::DSetCreatPropList prop{};
  prop.setChunk(1, &chunk_size);
  prop.setDeflate(compression_lvl);
  prop.setFillValue(type, &fill_value);  // TODO: make sure this works for STR_TYPE, otherwise
                                         // use const H5std_string FILL_VAL_STR{'\0'};

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
  try {
    H5::H5File f(path.c_str(), mode == READ_ONLY ? H5F_ACC_RDONLY : H5F_ACC_CREAT);
    if (validate) {
      validate_file_format(f, flavor, mode, bin_size, true);
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
  hdf5::read_attribute(f, "format", buff);
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

    hdf5::read_attribute(f, "format-version", int_buff, root_path);
    if (int_buff < min_format_ver && int_buff > max_format_ver) {
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
    auto check_int_dset = [&](std::string_view name) {
      if (!hdf5::dataset_exists(f, "chroms/length", H5::PredType::NATIVE_INT, root_path)) {
        if (!throw_on_failure) {
          return false;
        }
        throw std::runtime_error(fmt::format(FMT_STRING("Missing mandatory dataset '{}/{}'"),
                                             absl::StripPrefix(root_path, "/"), name));
      }
    };

    auto check_str_dset = [&](std::string_view name) {
      if (!hdf5::dataset_exists(f, "chroms/length", H5::PredType::C_S1, root_path)) {
        if (!throw_on_failure) {
          return false;
        }
        throw std::runtime_error(fmt::format(FMT_STRING("Missing mandatory dataset '{}/{}'"),
                                             absl::StripPrefix(root_path, "/"), name));
      }
    };

    ok = ok && check_int_dset("chroms/length");
    ok = ok && check_str_dset("chroms/name");
    ok = ok && check_int_dset("bins/chrom");
    ok = ok && check_int_dset("bins/start");
    ok = ok && check_int_dset("bins/end");
    ok = ok && check_int_dset("pixels/bin1_id");
    ok = ok && check_int_dset("pixels/bin2_id");
    ok = ok && check_int_dset("pixels/count");
    ok = ok && check_int_dset("indexes/bin1_offset");
    ok = ok && check_int_dset("indexes/chrom_offset");

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

  if (const auto format = hdf5::read_attribute_str(f, "format", root_path);
      absl::EndsWithIgnoreCase(format, "::mcool")) {
    if (!throw_on_failure) {
      return false;
    }
    throw std::runtime_error(
        fmt::format(FMT_STRING("File is not in Multires-Cooler (MCool) format: format attribute "
                               "should be HDF5::MCOOL, is '{}'"),
                    format));
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

Cooler::Cooler(std::string_view path_to_file, IO_MODE mode, std::size_t bin_size, Flavor flavor,
               bool validate, uint8_t compression_lvl, std::size_t chunk_size,
               std::size_t cache_size)
    : _path_to_file(std::string{path_to_file.data(), path_to_file.size()}),
      _mode(mode),
      _bin_size(bin_size),
      _flavor(flavor),
      _fp(open_file(_path_to_file, _mode, _bin_size, _flavor, validate)),
      _groups(open_groups(*_fp, !this->is_read_only())),
      _compression_lvl(compression_lvl),
      _chunk_size(chunk_size),
      _cache_size(cache_size),
      _cprop_str(this->is_read_only() ? nullptr
                                      : generate_default_cprop(_chunk_size, _compression_lvl,
                                                               Cooler::STR_TYPE, '\0')),
      _cprop_int32(this->is_read_only() ? nullptr
                                        : generate_default_cprop(_chunk_size, _compression_lvl,
                                                                 Cooler::INT32_TYPE, 0)),
      _cprop_int64(this->is_read_only() ? nullptr
                                        : generate_default_cprop(_chunk_size, _compression_lvl,
                                                                 Cooler::INT64_TYPE, 0L)),
      _aprop_str(this->is_read_only()
                     ? generate_default_aprop(Cooler::STR_TYPE, _chunk_size, _cache_size)
                     : nullptr),
      _aprop_int32(this->is_read_only()
                       ? generate_default_aprop(Cooler::INT32_TYPE, _chunk_size, _cache_size)
                       : nullptr),
      _aprop_int64(this->is_read_only()
                       ? generate_default_aprop(Cooler::INT64_TYPE, _chunk_size, _cache_size)
                       : nullptr) {}

bool Cooler::is_read_only() const { return this->_mode == Cooler::READ_ONLY; }

void Cooler::write_metadata() {
  if (this->is_read_only()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Caught attempt to write metadata to an HDF5 file that is opened in "
                               "read-only mode. File name: {}"),
                    this->_path_to_file.string()));
  }
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

  if (!assembly_name.empty()) {
    str_buff = std::string{assembly_name};
    att = this->_fp->createAttribute("assembly-name", METADATA_STR_TYPE, attr_space);
    att.write(METADATA_STR_TYPE, &str_buff);
  }

  att = this->_fp->createAttribute("generated-by", METADATA_STR_TYPE, attr_space);
  str_buff = fmt::format("ModLE-v{}", "0.0.1");  // TODO make ModLE ver a tunable
  att.write(METADATA_STR_TYPE, str_buff);

  att = this->_fp->createAttribute("creation-date", METADATA_STR_TYPE, attr_space);
  str_buff = absl::FormatTime(absl::Now(), absl::UTCTimeZone());
  att.write(METADATA_STR_TYPE, str_buff);
}

void Cooler::write_cmatrix_to_file(const ContactMatrix<I> &cmatrix, bool wipe_before_writing) {}

}  // namespace modle::cooler