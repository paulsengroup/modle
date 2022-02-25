// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <absl/strings/str_cat.h>  // for StrCat, StrAppend
#include <absl/types/variant.h>    // for visit
#include <fmt/ostream.h>
#include <spdlog/spdlog.h>  // for error, warn

#include <boost/filesystem/path.hpp>  // for path
#include <cassert>                    // for assert
#include <exception>                  // for exception
#include <memory>                     // for make_unique
#include <tuple>                      // for ignore

#include "modle/common/common.hpp"  // for i32, i64, usize
#include "modle/hdf5/hdf5.hpp"      // for read_attribute, read_numbers, wri...

namespace modle::cooler {

template <class N>
Cooler<N>::Cooler(boost::filesystem::path path_to_file, IO_MODE mode, usize bin_size,
                  usize max_str_length, std::string_view assembly_name, FLAVOR flavor,
                  bool validate, u8f compression_lvl, usize chunk_size, usize cache_size)
    : STR_TYPE(generate_default_str_type(max_str_length)),
      _path_to_file(std::move(path_to_file)),
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
      _cprop_int32(this->is_read_only()
                       ? nullptr
                       : generate_default_cprop(_chunk_size, _compression_lvl,
                                                H5::PredType::NATIVE_INT32, i32(0))),
      _cprop_int64(this->is_read_only()
                       ? nullptr
                       : generate_default_cprop(_chunk_size, _compression_lvl,
                                                H5::PredType::NATIVE_INT64, i64(0))),
      _cprop_float64(this->is_read_only()
                         ? nullptr
                         : generate_default_cprop(_chunk_size, _compression_lvl,
                                                  H5::PredType::NATIVE_DOUBLE, double(0))),
      _aprop_str(generate_default_aprop(Cooler::STR_TYPE, _chunk_size, _cache_size)),
      _aprop_int32(generate_default_aprop(H5::PredType::NATIVE_INT32, _chunk_size, _cache_size)),
      _aprop_int64(generate_default_aprop(H5::PredType::NATIVE_INT64, _chunk_size, _cache_size)),
      _aprop_float64(
          generate_default_aprop(H5::PredType::NATIVE_DOUBLE, _chunk_size, _cache_size)) {
  assert(this->_flavor != FLAVOR::UNK);
  if (this->is_read_only() && this->_flavor == FLAVOR::AUTO) {
    this->_flavor = Cooler::detect_file_flavor(*this->_fp);
  }
  if (this->is_mcool()) {
    assert(this->_bin_size != 0);
    absl::StrAppend(&this->_root_path, "resolutions/", this->_bin_size, "/");
  }
  if (this->is_read_only()) {
    if (this->_bin_size == 0) {  // i.e. file is cooler
      assert(this->is_cool());
      DISABLE_WARNING_PUSH
      DISABLE_WARNING_USELESS_CAST
      this->_bin_size =
          static_cast<usize>(hdf5::read_attribute_int(*this->_fp, "bin-size", this->_root_path));
      this->_nnz = static_cast<i64>(hdf5::read_attribute_int(*this->_fp, "nnz", this->_root_path));
      if (hdf5::has_attribute(*this->_fp, "sum", this->_root_path)) {
        using SumT = decltype(this->_sum);
        if constexpr (IS_FP) {
          this->_sum =
              static_cast<SumT>(hdf5::read_attribute<double>(*this->_fp, "sum", this->_root_path));
        } else {
          this->_sum =
              static_cast<SumT>(hdf5::read_attribute_int(*this->_fp, "sum", this->_root_path));
        }
      }
      DISABLE_WARNING_POP
    }
    this->open_default_datasets();
    this->read_chrom_offset_idx();
    this->read_bin1_offset_idx();
  } else {
    this->_buff = std::make_unique<InternalBuffers>();
    this->init_default_datasets();
    this->write_metadata();
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
template <class N>
Cooler<N>::~Cooler() {
  try {
    if (!this->is_read_only() && this->_nchroms != 0 && this->_nbins != 0) {
      if (this->_fp) {
        auto &chrom_idx = this->_datasets[IDX_CHR];
        auto &bin1_idx = this->_datasets[IDX_BIN1];
        auto chrom_idx_offset = this->_dataset_file_offsets[IDX_CHR];
        auto bin1_idx_offset = this->_dataset_file_offsets[IDX_BIN1];

        if (chrom_idx_offset != 0) {
          decltype(this->_nbins) buff;
          std::ignore = hdf5::read_number(chrom_idx, buff, chrom_idx_offset - 1);
          if (buff != this->_nbins) {
            assert(buff < this->_nbins);
            std::ignore = hdf5::write_number(this->_nbins, chrom_idx, chrom_idx_offset);
          }
        }
        if (bin1_idx_offset != 0) {
          decltype(this->_nnz) buff;
          std::ignore = hdf5::read_number(bin1_idx, buff, bin1_idx_offset - 1);
          if (buff != this->_nnz) {
            assert(buff < this->_nnz);
            std::ignore = hdf5::write_number(this->_nnz, bin1_idx, bin1_idx_offset);
          }
        }
      } else {
        spdlog::error(
            FMT_STRING(
                "Message for the developers: ~Cooler() for file {} was called on a closed "
                "file handle. This should never happen, as whenever Cooler objects are created "
                "in write-mode the file handle is supposed to be closed upon object "
                "destruction!"),
            this->_path_to_file);
      }
    }
  } catch (const H5::Exception &e) {
    spdlog::error(FMT_STRING("The following error occurred while finalizing file {}: {}"),
                  this->_path_to_file, hdf5::construct_error_stack(e));
    spdlog::error(FMT_STRING("The content of file {} may be corrupted or incomplete."),
                  this->_path_to_file);
  } catch (const std::exception &e) {
    spdlog::error(
        FMT_STRING("Caught the following unhandled exception while finalizing file {}: {}"),
        this->_path_to_file, e.what());
    spdlog::error(FMT_STRING("The content of file {} may be corrupted or incomplete."),
                  this->_path_to_file);
  }
}

template <class N>
constexpr bool Cooler<N>::is_read_only() const noexcept {
  return this->_mode == IO_MODE::READ_ONLY;
}

template <class N>
const boost::filesystem::path &Cooler<N>::get_path() const {
  return this->_path_to_file;
}

template <class N>
usize Cooler<N>::get_nchroms() {
  assert(this->_fp);
  if (this->is_cool()) {
    return static_cast<usize>(hdf5::read_attribute_int(*this->_fp, "nchroms"));
  }
  if (this->is_mcool()) {
    assert(this->_bin_size != 0);
    return static_cast<usize>(hdf5::read_attribute_int(
        *this->_fp, "nchroms", absl::StrCat("/resolutions/", this->_bin_size)));
  }
  MODLE_UNREACHABLE_CODE;
}

template <class N>
void Cooler<N>::get_chrom_names(std::vector<std::string> &buff) {
  assert(this->_fp);
  const auto nchroms = this->get_nchroms();
  buff.resize(nchroms);
  if (!this->_datasets.empty()) {
    const auto &d = this->_datasets[CHROM_NAME];
    std::ignore = hdf5::read_strings(d, buff, 0);
  } else {
    if (this->is_cool()) {
      auto d = hdf5::open_dataset(*this->_fp, "/chroms/name", *this->_aprop_str);
      std::ignore = hdf5::read_strings(d, buff, 0);
    } else if (this->is_mcool()) {
      assert(this->_bin_size != 0);
      auto d = hdf5::open_dataset(*this->_fp,
                                  absl::StrCat("/resolutions/", this->_bin_size, "/chroms/name"),
                                  *this->_aprop_str);
      std::ignore = hdf5::read_strings(d, buff, 0);
    } else {
      assert(!this->is_scool());
      buff.clear();
    }
  }
}

template <class N>
std::vector<std::string> Cooler<N>::get_chrom_names() {
  std::vector<std::string> buff;
  this->get_chrom_names(buff);
  return buff;
}

template <class N>
template <class I>
void Cooler<N>::get_chrom_sizes(std::vector<I> &buff) {
  static_assert(std::is_integral_v<I>, "buff should be a vector of integers.");
  assert(this->_fp);
  const auto nchroms = this->get_nchroms();
  buff.resize(nchroms);
  if (!this->_datasets.empty()) {
    const auto &d = this->_datasets[CHROM_LEN];
    std::ignore = hdf5::read_numbers(d, buff, 0);
  } else {
    if (this->is_cool()) {
      auto d = hdf5::open_dataset(*this->_fp, "/chroms/length", *this->_aprop_str);
      std::ignore = hdf5::read_numbers(d, buff, 0);
    } else if (this->is_mcool()) {
      assert(this->_bin_size != 0);
      auto d = hdf5::open_dataset(*this->_fp,
                                  absl::StrCat("/resolutions/", this->_bin_size, "/chroms/length"),
                                  *this->_aprop_str);
      std::ignore = hdf5::read_numbers(d, buff, 0);
    } else {
      assert(!this->is_scool());
      buff.clear();
    }
  }
}

template <class N>
std::vector<i64> Cooler<N>::get_chrom_sizes() {
  std::vector<i64> buff;
  this->get_chrom_sizes(buff);
  return buff;
}

template <class N>
std::vector<std::pair<std::string, usize>> Cooler<N>::get_chroms() {
  std::vector<std::string> name_buff;
  std::vector<usize> size_buff;
  this->get_chrom_names(name_buff);
  this->get_chrom_sizes(size_buff);
  assert(name_buff.size() == size_buff.size());
  std::vector<std::pair<std::string, usize>> buff(name_buff.size());
  for (usize i = 0; i < buff.size(); ++i) {
    buff[i].first = std::move(name_buff[i]);
    buff[i].second = size_buff[i];
  }
  return buff;
}

template <class N>
constexpr bool Cooler<N>::is_cool() const noexcept {
  return this->_flavor == FLAVOR::COOL;
}
template <class N>
constexpr bool Cooler<N>::is_mcool() const noexcept {
  return this->_flavor == FLAVOR::MCOOL;
}
template <class N>
constexpr bool Cooler<N>::is_scool() const noexcept {
  return this->_flavor == FLAVOR::SCOOL;
}

template <class N>
constexpr usize Cooler<N>::get_bin_size() const noexcept {
  return this->_bin_size;
}

template <class N>
bool Cooler<N>::has_contacts_for_chrom(std::string_view chrom_name,
                                       bool try_common_chrom_prefixes) {
  assert(this->_fp);
  const auto chrom_idx = this->get_chrom_idx(chrom_name, try_common_chrom_prefixes);
  return this->has_contacts_for_chrom(chrom_idx);
}

template <class N>
bool Cooler<N>::has_contacts_for_chrom(usize chrom_idx) const {
  assert(this->_fp);
  assert(this->is_read_only());
  assert(!this->_idx_bin1_offset.empty());
  assert(!this->_idx_chrom_offset.empty());
  assert(chrom_idx < this->_idx_chrom_offset.size());

  const auto first_bin = static_cast<usize>(this->_idx_chrom_offset[chrom_idx]);
  const auto last_bin = static_cast<usize>(this->_idx_chrom_offset[chrom_idx + 1]);
  assert(last_bin >= first_bin);

  return this->_idx_bin1_offset[first_bin] != this->_idx_bin1_offset[last_bin];
}

template <class N>
bool Cooler<N>::validate_file_format(H5::H5File &f, FLAVOR expected_flavor, IO_MODE mode,
                                     usize bin_size, bool throw_on_failure) {
  if (mode == IO_MODE::WRITE_ONLY) {  // File is empty. Nothing to validate here
    return true;
  }

  try {
    const auto flavor = detect_file_flavor(f);
    if (expected_flavor != FLAVOR::AUTO && expected_flavor != FLAVOR::UNK &&
        expected_flavor != flavor) {
      throw std::runtime_error(fmt::format(FMT_STRING("Expected format flavor {}, found {}"),
                                           flavor_to_string(expected_flavor),
                                           flavor_to_string(flavor)));
    }
    switch (flavor) {
      case FLAVOR::COOL:
        return validate_cool_flavor(f, bin_size, "/", throw_on_failure);
      case FLAVOR::MCOOL:
        return validate_multires_cool_flavor(f, bin_size, "/", throw_on_failure);
      case FLAVOR::SCOOL:
        throw std::runtime_error("SCOOL flavor is not yet supported");
      default:
        MODLE_UNREACHABLE_CODE;
    }

  } catch (const std::runtime_error &e) {
    throw std::runtime_error(fmt::format(FMT_STRING("Cooler validation for file \"{}\" failed: {}"),
                                         hdf5::get_file_name(f), e.what()));
  }
}

template <class N>
H5::StrType Cooler<N>::generate_default_str_type(usize max_str_length) {
  // Cooltools doesn't seem to properly handle variable length strings (H5T_VARIABLE)
  // For the time being we are forced to use fixed length, null-padded strings
  auto st = max_str_length > 0 ? H5::StrType(H5::PredType::C_S1, max_str_length)
                               : H5::StrType(H5::PredType::C_S1);
  st.setStrpad(H5T_STR_NULLPAD);
  st.setCset(H5T_CSET_ASCII);
  return st;
}

template <class N>
template <class T1, class T2>
std::unique_ptr<H5::DSetCreatPropList> Cooler<N>::generate_default_cprop(
    hsize_t chunk_size, u8 compression_lvl, [[maybe_unused]] T1 type,
    [[maybe_unused]] T2 fill_value) {
  static_assert(
      std::is_same_v<T2, i64> || std::is_same_v<T2, i32> || std::is_same_v<T2, double> ||
          std::is_constructible_v<H5std_string, T2>,
      "fill_value should have one of the following types: i64, i32, double or std::string.");
  static_assert(
      (std::is_constructible_v<H5std_string, T2> && std::is_same_v<T1, H5::StrType>) ||
          std::is_same_v<T1, H5::PredType>,
      "Incompatible data type for variables type and fill_value: if T2 is string "
      "constructible, then T1 must be H5::StrType, else T2 is integral and type is H5::PredType");

  return std::make_unique<H5::DSetCreatPropList>(
      hdf5::generate_creat_prop_list(chunk_size, compression_lvl, type, fill_value));
}

template <class N>
template <class T>
std::unique_ptr<H5::DSetAccPropList> Cooler<N>::generate_default_aprop([[maybe_unused]] T type,
                                                                       hsize_t chunk_size,
                                                                       hsize_t cache_size) {
  static_assert(std::is_same_v<T, H5::StrType> || std::is_same_v<T, H5::PredType>,
                "type should be of type H5::StrType or H5::PredType");
  [[maybe_unused]] const auto rdcc_w0_streaming{0.99};
  [[maybe_unused]] const auto rdcc_w0_random{0.01};

  chunk_size = [&]() {
    if constexpr (std::is_same_v<T, H5::StrType>) {
      return chunk_size;
    }
    if (type == H5::PredType::NATIVE_INT64 || type == H5::PredType::NATIVE_UINT64 ||
        type == H5::PredType::NATIVE_DOUBLE) {
      return chunk_size;
    }
    if (type == H5::PredType::NATIVE_INT32 || type == H5::PredType::NATIVE_UINT32) {
      return chunk_size / 2;
    }
    if (type == H5::PredType::NATIVE_INT16 || type == H5::PredType::NATIVE_UINT16) {
      return chunk_size / 4;
    }
    if (type == H5::PredType::NATIVE_INT8 || type == H5::PredType::NATIVE_UINT8) {
      return chunk_size / 8;
    }
    throw std::runtime_error(
        "Cooler<N>::generate_default_aprop() was called with an unsupported H5 PredType/StrType");
  }();

  return std::make_unique<H5::DSetAccPropList>(
      hdf5::generate_acc_prop_list(type, chunk_size, cache_size, rdcc_w0_streaming));
}

template <class N>
std::unique_ptr<H5::H5File> Cooler<N>::open_file(const boost::filesystem::path &path, IO_MODE mode,
                                                 usize bin_size,
                                                 [[maybe_unused]] usize max_str_length,
                                                 FLAVOR flavor, bool validate) {
  assert(!path.has_parent_path() || boost::filesystem::is_directory(path.parent_path()));
  if constexpr (utils::ndebug_not_defined()) {
    if (mode == IO_MODE::WRITE_ONLY) {
      if (bin_size == 0) {
        throw std::runtime_error(
            "Cooler::open_file(): bin_size cannot be 0 when file is being opened "
            "in WRITE_ONLY mode");
      }
      if (max_str_length == 0) {
        throw std::runtime_error(
            "Cooler::open_file(): max_str_length cannot be 0 when file is being "
            "opened in WRITE_ONLY mode");
      }
    }
  }
  try {
    auto f = mode == IO_MODE::READ_ONLY ? hdf5::open_file_for_reading(path)
                                        : hdf5::open_file_for_writing(path);
    if (validate) {
      std::ignore = validate_file_format(f, flavor, mode, bin_size, true);
    }
    return std::make_unique<H5::H5File>(f);
  } catch (const H5::Exception &e) {
    const auto error_msg = hdf5::construct_error_stack();
    if (absl::StrContains(error_msg, "Unable to open file")) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("Unable to open file {} for {}"), path,
          mode == IO_MODE::READ_ONLY
              ? "reading"
              : "writing. Please ensure that the file is not in use by another application"));
    }
    throw std::runtime_error(hdf5::construct_error_stack(e));
  }
}

template <class N>
std::vector<H5::Group> Cooler<N>::open_groups(H5::H5File &f, bool create_if_not_exist,
                                              usize bin_size) {
  std::vector<H5::Group> groups(4);

  const std::string root_path = bin_size != 0 && hdf5::has_group(f, "/resolutions")
                                    ? absl::StrCat("/resolutions/", bin_size, "/")
                                    : "";
  try {
    if (create_if_not_exist) {
      groups[chrom] = hdf5::open_or_create_group(f, absl::StrCat(root_path, "chroms"));
      groups[BIN] = hdf5::open_or_create_group(f, absl::StrCat(root_path, "bins"));
      groups[PXL] = hdf5::open_or_create_group(f, absl::StrCat(root_path, "pixels"));
      groups[IDX] = hdf5::open_or_create_group(f, absl::StrCat(root_path, "indexes"));
    } else {
      groups[chrom] = hdf5::open_group(f, absl::StrCat(root_path, "chroms"));
      groups[BIN] = hdf5::open_group(f, absl::StrCat(root_path, "bins"));
      groups[PXL] = hdf5::open_group(f, absl::StrCat(root_path, "pixels"));
      groups[IDX] = hdf5::open_group(f, absl::StrCat(root_path, "indexes"));
    }

    return groups;
  } catch ([[maybe_unused]] const H5::Exception &e) {
    throw std::runtime_error(hdf5::construct_error_stack());
  }
}

template <class N>
usize Cooler<N>::read_chrom_offset_idx() {
  const auto &d = this->_datasets[IDX_CHR];
  const auto buff_size = static_cast<hsize_t>(d.getSpace().getSimpleExtentNpoints());
  this->_idx_chrom_offset.resize(buff_size);

  const auto idx_size = hdf5::read_numbers(d, this->_idx_chrom_offset, 0);
  if (idx_size != this->_idx_chrom_offset.size()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("An error occurred while reading dataset 'indexes/chrom_offset' "
                               "from file \"{}\": expected to read {} numbers, but only read {}"),
                    this->_path_to_file, buff_size, idx_size));
  }
  return idx_size;
}

template <class N>
usize Cooler<N>::read_bin1_offset_idx() {
  const auto &d = this->_datasets[IDX_BIN1];
  const auto buff_size = static_cast<hsize_t>(d.getSpace().getSimpleExtentNpoints());
  this->_idx_bin1_offset.resize(buff_size);

  const auto idx_size = hdf5::read_numbers(d, this->_idx_bin1_offset, 0);
  if (idx_size != this->_idx_bin1_offset.size()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("An error occurred while reading dataset 'indexes/chrom_offset' "
                               "from file \"{}\": expected to read {} numbers, but only read {}"),
                    this->_path_to_file, buff_size, idx_size));
  }
  return idx_size;
}

template <class N>
absl::Span<const i64> Cooler<N>::get_bin1_offset_idx_for_chrom(
    usize chrom_idx, std::pair<usize, usize> chrom_subrange) {
  assert(!this->_idx_bin1_offset.empty());
  assert(!this->_idx_chrom_offset.empty());
  assert(chrom_idx < this->_idx_chrom_offset.size());
  const auto chrom_start_bin = static_cast<usize>(this->_idx_chrom_offset[chrom_idx]) +
                               (chrom_subrange.first / this->_bin_size);
  const auto chrom_end_bin =
      chrom_subrange.second == std::numeric_limits<decltype(chrom_subrange.second)>::max()
          ? static_cast<usize>(this->_idx_chrom_offset[chrom_idx + 1])
          : static_cast<usize>(this->_idx_chrom_offset[chrom_idx]) +
                ((chrom_subrange.second + this->_bin_size - 1) / this->_bin_size);
  DISABLE_WARNING_PUSH
  DISABLE_WARNING_SIGN_COMPARE
  assert(chrom_end_bin <= this->_idx_chrom_offset[chrom_idx + 1]);
  assert(chrom_end_bin >= chrom_start_bin);
  DISABLE_WARNING_POP

  return absl::MakeConstSpan(this->_idx_bin1_offset)
      .subspan(chrom_start_bin, chrom_end_bin - chrom_start_bin + 1);
}

template <class N>
absl::Span<const i64> Cooler<N>::get_bin1_offset_idx_for_chrom(
    std::string_view chrom_name, std::pair<usize, usize> chrom_subrange) {
  const auto chrom_idx = get_chrom_idx(chrom_name);
  return get_bin1_offset_idx_for_chrom(chrom_idx, chrom_subrange);
}

template <class N>
std::pair<i64, i64> Cooler<N>::read_chrom_pixels_boundaries(std::string_view chrom_name) {
  const auto chrom_idx = get_chrom_idx(chrom_name);
  return read_chrom_pixels_boundaries(chrom_idx);
}

template <class N>
std::pair<i64, i64> Cooler<N>::read_chrom_pixels_boundaries(usize chrom_idx) {
  assert(chrom_idx < this->_idx_chrom_offset.size());
  assert(!this->_idx_chrom_offset.empty());
  const auto chrom_start_bin = static_cast<usize>(this->_idx_chrom_offset[chrom_idx]);
  const auto chrom_end_bin = static_cast<usize>(this->_idx_chrom_offset[chrom_idx + 1]);
  assert(chrom_end_bin >= chrom_start_bin);

  std::pair<i64, i64> pixel_boundaries{};
  const auto &d = this->_datasets[IDX_BIN1];
  std::ignore = hdf5::read_number(d, pixel_boundaries.first, chrom_start_bin);
  std::ignore = hdf5::read_number(d, pixel_boundaries.second, chrom_end_bin);

  return pixel_boundaries;
}

template <class N>
void Cooler<N>::init_default_datasets() {
  this->_datasets.resize(DEFAULT_DATASETS_NR);
  this->_dataset_file_offsets.resize(DEFAULT_DATASETS_NR);
  std::fill(this->_dataset_file_offsets.begin(), this->_dataset_file_offsets.end(), 0);

  assert(this->_cprop_str);
  assert(this->_cprop_int32);
  assert(this->_cprop_int64);

  assert(this->_aprop_str);
  assert(this->_aprop_int32);
  assert(this->_aprop_int64);

  const auto &cpi64 = *this->_cprop_int64;
  [[maybe_unused]] const auto &cpi32 = *this->_cprop_int32;
  [[maybe_unused]] const auto &cpf64 = *this->_cprop_float64;
  const auto &cps = *this->_cprop_str;

  const auto &api64 = *this->_aprop_int64;
  [[maybe_unused]] const auto &api32 = *this->_aprop_int32;
  [[maybe_unused]] const auto &apf64 = *this->_aprop_float64;
  const auto &aps = *this->_aprop_str;

  auto &f = *this->_fp;
  auto &dset = this->_datasets;
  auto r = absl::StripSuffix(this->_root_path, "/");

  try {
    // Do not change the order of these pushbacks
    dset[CHROM_LEN] = hdf5::create_dataset(f, absl::StrCat(r, "/chroms/length"),
                                           H5::PredType::NATIVE_INT64, cpi64, api64);
    dset[CHROM_NAME] = hdf5::create_dataset(f, absl::StrCat(r, "/chroms/name"), STR_TYPE, cps, aps);

    dset[BIN_CHROM] = hdf5::create_dataset(f, absl::StrCat(r, "/bins/chrom"),
                                           H5::PredType::NATIVE_INT64, cpi64, api64);
    dset[BIN_START] = hdf5::create_dataset(f, absl::StrCat(r, "/bins/start"),
                                           H5::PredType::NATIVE_INT64, cpi64, api64);
    dset[BIN_END] = hdf5::create_dataset(f, absl::StrCat(r, "/bins/end"),
                                         H5::PredType::NATIVE_INT64, cpi64, api64);

    dset[PXL_B1] = hdf5::create_dataset(f, absl::StrCat(r, "/pixels/bin1_id"),
                                        H5::PredType::NATIVE_INT64, cpi64, api64);
    dset[PXL_B2] = hdf5::create_dataset(f, absl::StrCat(r, "/pixels/bin2_id"),
                                        H5::PredType::NATIVE_INT64, cpi64, api64);
    if constexpr (std::is_floating_point_v<N>) {
      dset[PXL_COUNT] = hdf5::create_dataset(f, absl::StrCat(r, "/pixels/count"),
                                             H5::PredType::NATIVE_DOUBLE, cpf64, apf64);
    } else {
      dset[PXL_COUNT] = hdf5::create_dataset(f, absl::StrCat(r, "/pixels/count"),
                                             H5::PredType::NATIVE_INT32, cpi32, api32);
    }

    dset[IDX_BIN1] = hdf5::create_dataset(f, absl::StrCat(r, "/indexes/bin1_offset"),
                                          H5::PredType::NATIVE_INT64, cpi64, api64);
    dset[IDX_CHR] = hdf5::create_dataset(f, absl::StrCat(r, "/indexes/chrom_offset"),
                                         H5::PredType::NATIVE_INT64, cpi64, api64);

  } catch ([[maybe_unused]] const H5::FileIException &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING(
            "An error occurred while initializing default Cooler dataset on file \"{}\": {}"),
        this->_path_to_file.string(), hdf5::construct_error_stack()));
  }
}

template <class N>
void Cooler<N>::open_default_datasets() {
  assert(this->_fp);
  assert(this->_aprop_int32);
  assert(this->_aprop_int64);
  assert(this->_aprop_float64);
  assert(this->_aprop_str);

  auto &d = this->_datasets;
  auto &f = *this->_fp;
  [[maybe_unused]] const auto &ai32 = *this->_aprop_int32;
  const auto &ai64 = *this->_aprop_int64;
  [[maybe_unused]] const auto &af64 = *this->_aprop_float64;
  const auto &achrom_name = *this->_aprop_str;
  const auto &achrom_size = *this->_aprop_str;

  this->_datasets.resize(DEFAULT_DATASETS_NR);
  this->_dataset_file_offsets.resize(DEFAULT_DATASETS_NR);
  std::fill(this->_dataset_file_offsets.begin(), this->_dataset_file_offsets.end(), 0);
  try {
    d[CHROM_LEN] =
        hdf5::open_dataset(f, absl::StrCat(this->_root_path, "chroms/length"), achrom_size);
    d[CHROM_NAME] =
        hdf5::open_dataset(f, absl::StrCat(this->_root_path, "chroms/name"), achrom_name);

    d[BIN_CHROM] = hdf5::open_dataset(f, absl::StrCat(this->_root_path, "bins/chrom"), ai64);
    d[BIN_START] = hdf5::open_dataset(f, absl::StrCat(this->_root_path, "bins/start"), ai64);
    d[BIN_END] = hdf5::open_dataset(f, absl::StrCat(this->_root_path, "bins/end"), ai64);
    if (hdf5::has_dataset(f, "bins/weight", this->_root_path)) {
      d[BIN_WEIGHT] = hdf5::open_dataset(f, absl::StrCat(this->_root_path, "bins/weight"), af64);
    }

    d[PXL_B1] = hdf5::open_dataset(f, absl::StrCat(this->_root_path, "pixels/bin1_id"), ai64);
    d[PXL_B2] = hdf5::open_dataset(f, absl::StrCat(this->_root_path, "pixels/bin2_id"), ai64);
    if constexpr (std::is_floating_point_v<N>) {
      d[PXL_COUNT] = hdf5::open_dataset(f, absl::StrCat(this->_root_path, "pixels/count"), af64);
    } else {
      d[PXL_COUNT] = hdf5::open_dataset(f, absl::StrCat(this->_root_path, "pixels/count"), ai32);
    }

    d[IDX_BIN1] =
        hdf5::open_dataset(f, absl::StrCat(this->_root_path, "indexes/bin1_offset"), ai64);
    d[IDX_CHR] =
        hdf5::open_dataset(f, absl::StrCat(this->_root_path, "indexes/chrom_offset"), ai64);

  } catch ([[maybe_unused]] const H5::FileIException &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("An error occurred while trying to open Cooler's default datasets "
                               "of file \"{}\": {}"),
                    this->_path_to_file, hdf5::construct_error_stack()));
  }
}

template <class N>
usize Cooler<N>::get_chrom_idx(std::string_view query_chrom_name, bool try_common_chrom_prefixes) {
  // Here's the issue: for a given genome assembly (say hg19), some tools name chromosome as
  // chrN, where N is the chromosome number, while others only use the number N. The purpose
  // of this function is to try to guess few reasonably common prefixes, look them up in the
  // Cooler file, then return the chromosome name variant that produced a hit together with
  // the chrom index

  try {
    const auto chrom_names = hdf5::read_strings(this->_datasets[CHROM_NAME], 0);
    auto match = std::find(chrom_names.begin(), chrom_names.end(), query_chrom_name);
    if (match != chrom_names.end()) {
      return static_cast<usize>(std::distance(chrom_names.begin(), match));
    }
    if (try_common_chrom_prefixes) {
      std::array<std::string, 3> queries;
      constexpr std::array<std::string_view, 3> prefixes = {"chr", "CHR", "Chr"};
      std::transform(prefixes.begin(), prefixes.end(), queries.begin(), [&](const auto prefix) {
        if (absl::StartsWith(query_chrom_name, prefix)) {
          return std::string{absl::StripPrefix(query_chrom_name, prefix)};
        }
        return absl::StrCat(prefix, query_chrom_name);
      });

      for (const auto &q : queries) {
        match = std::find(chrom_names.begin(), chrom_names.end(), q);
        if (match != chrom_names.end()) {
          return static_cast<usize>(std::distance(chrom_names.begin(), match));
        }
      }
      throw std::runtime_error(fmt::format(
          FMT_STRING("Unable to find a chromosome named \"{}\". The following chromosome "
                     "name variants were searched: \"{}\""),
          query_chrom_name, fmt::join(queries, "', '")));
    }
    throw std::runtime_error(
        fmt::format(FMT_STRING("Unable to find a chromosome named \"{}\""), query_chrom_name));
  } catch (const std::exception &e) {
    if (absl::StartsWith(e.what(), "Unable to find a chromosome")) {
      throw;
    }
    throw std::runtime_error(
        fmt::format(FMT_STRING("The following error occurred while looking up "
                               "\"{}\" in dataset chroms/name of file \"{}\": {}"),
                    query_chrom_name, this->_path_to_file, e.what()));
  }
}

template <class N>
std::string Cooler<N>::flavor_to_string(FLAVOR f) {
  switch (f) {
    case FLAVOR::UNK:
      return "unknown";
    case FLAVOR::AUTO:
      return "auto";
    case FLAVOR::COOL:
      return "COOL";
    case FLAVOR::MCOOL:
      return "MCOOL";
    case FLAVOR::SCOOL:
      return "SCOOL";
  }
  return "";
}

template <class N>
typename Cooler<N>::FLAVOR Cooler<N>::detect_file_flavor(H5::H5File &f) {
  if (const auto &fname = hdf5::get_file_name(f); absl::EndsWithIgnoreCase(fname, ".cool")) {
    return FLAVOR::COOL;  // NOLINTNEXTLINE(readability-else-after-return)
  } else if (absl::EndsWithIgnoreCase(fname, ".mcool")) {
    return FLAVOR::MCOOL;
  } else if (absl::EndsWithIgnoreCase(fname, ".scool")) {
    return FLAVOR::SCOOL;
  }

  std::string buff;
  if (hdf5::has_attribute(f, "format")) {
    hdf5::read_attribute(f, "format", buff);
  }
  if (absl::EndsWithIgnoreCase(buff, "::cooler")) {
    return FLAVOR::COOL;
  }
  if (absl::EndsWithIgnoreCase(buff, "::mcool")) {
    return FLAVOR::MCOOL;
  }
  if (absl::EndsWithIgnoreCase(buff, "::scool")) {
    return FLAVOR::SCOOL;
  }

  if (hdf5::has_group(f, "/resolutions")) {
    return FLAVOR::MCOOL;
  }
  if (hdf5::has_group(f, "/bins") && hdf5::has_group(f, "/chroms")) {
    if (hdf5::has_group(f, "/pixels") && hdf5::has_group(f, "/indexes")) {
      return FLAVOR::COOL;
    }
    if (hdf5::has_group(f, "/cells")) {
      return FLAVOR::SCOOL;
    }
  }
  throw std::runtime_error("Unable to detect Cooler file flavor");
}

// I am not sure there's much to be done to reduce the complexity of this function
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
template <class N>
bool Cooler<N>::validate_cool_flavor(H5::H5File &f, usize bin_size, std::string_view root_path,
                                     bool throw_on_failure, bool check_version) {
  /* The following attributes are being checked:
   * format
   * format-version
   * bin-type
   * bin size
   * storage-mode
   */

  std::string str_buff{};
  i64 int_buff{};

  constexpr i64 min_format_ver = 2;
  constexpr i64 max_format_ver = 3;

  try {
    // Checking attributes
    if (hdf5::has_attribute(f, "format", root_path)) {
      hdf5::read_attribute(f, "format", str_buff, root_path);
      if (!absl::EndsWithIgnoreCase(str_buff, "::cooler")) {
        if (!throw_on_failure) {
          return false;
        }
        throw std::runtime_error(
            fmt::format(FMT_STRING("File is not in Cooler format: format attribute should be "
                                   "HDF5::Cooler, is \"{}\""),
                        str_buff));
      }
    } else {
      spdlog::warn(FMT_STRING("WARNING: missing attribute 'format' in file \"{}\"\n"),
                   hdf5::get_file_name(f));
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
    const auto format_version = static_cast<u8>(int_buff);
    int_buff = 0;

    hdf5::read_attribute(f, "bin-type", str_buff, root_path);
    if (str_buff != "fixed") {
      if (!throw_on_failure) {
        return false;
      }
      throw std::runtime_error(fmt::format(
          FMT_STRING("Expected bin-type attribute to be 'fixed', got \"{}\""), str_buff));
    }
    str_buff.clear();

    if (bin_size != 0) {
      hdf5::read_attribute(f, "bin-size", int_buff, root_path);
      if (static_cast<i64>(bin_size) != int_buff) {
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
            FMT_STRING("Expected storage-mode attribute to be 'symmetric-upper', got \"{}\""),
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
      auto attr_name = msg.substr(0, msg.find_first_of('\''));
      throw std::runtime_error(
          fmt::format(FMT_STRING("Missing mandatory attribute \"{}\"{}"), attr_name,
                      root_path == "/" ? "" : absl::StrCat(" from path '", root_path, "'")));
    }
    throw;
  }
  return true;
}

template <class N>
bool Cooler<N>::validate_multires_cool_flavor(H5::H5File &f, usize bin_size,
                                              std::string_view root_path, bool throw_on_failure) {
  constexpr i64 min_format_ver = 2;
  constexpr i64 max_format_ver = 3;

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
                                 "should be HDF5::MCOOL, is \"{}\""),
                      format));
    }
  } else {
    spdlog::warn(FMT_STRING("WARNING: missing attribute 'format' in file \"{}\"\n"),
                 hdf5::get_file_name(f));
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

template <class N>
Cooler<N>::InternalBuffers::InternalBuffers(usize buff_size) {
  bin_pos_buff.reserve(buff_size);
  bin_chrom_buff.reserve(buff_size);
  pixel_b1_idx_buff.reserve(buff_size);
  pixel_b2_idx_buff.reserve(buff_size);
  pixel_count_buff.reserve(buff_size);
  idx_bin1_offset_buff.reserve(buff_size);
  idx_chrom_offset_buff.reserve(buff_size);
}

template <class N>
usize Cooler<N>::InternalBuffers::capacity() const {
  return this->bin_pos_buff.capacity();
}

}  // namespace modle::cooler
