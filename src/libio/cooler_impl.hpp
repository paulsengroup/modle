#pragma once

#include <H5Cpp.h>
#include <absl/strings/match.h>
#include <absl/types/span.h>
#include <fmt/format.h>

#include <algorithm>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include "modle/contacts.hpp"
#include "modle/hdf5.hpp"

namespace modle::cooler {

Cooler::Cooler(std::filesystem::path path_to_file, IO_MODE mode, std::size_t bin_size,
               std::size_t max_str_length, std::string_view assembly_name, Flavor flavor,
               bool validate, uint8_t compression_lvl, std::size_t chunk_size,
               std::size_t cache_size)
    : STR_TYPE(generate_default_str_type(max_str_length)),
      _path_to_file(path_to_file),
      _mode(mode),
      _bin_size(bin_size),
      _assembly_name(assembly_name.data(), assembly_name.size()),
      _flavor(flavor),
      _fp(open_file(_path_to_file, _mode, _bin_size, _flavor, validate)),
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
    this->read_chr_offset_idx();
    this->read_bin1_offset_idx();
  } else {
    this->init_default_datasets();
    this->write_metadata();
  }
}

Cooler::Cooler(const std::string &path_to_file, IO_MODE mode, std::size_t bin_size,
               std::size_t max_str_length, std::string_view assembly_name, Flavor flavor,
               bool validate, uint8_t compression_lvl, std::size_t chunk_size,
               std::size_t cache_size)
    : Cooler::Cooler(std::filesystem::path(path_to_file), mode, bin_size, max_str_length,
                     assembly_name, flavor, validate, compression_lvl, chunk_size, cache_size) {}

Cooler::Cooler(std::string_view path_to_file, IO_MODE mode, std::size_t bin_size,
               std::size_t max_str_length, std::string_view assembly_name, Flavor flavor,
               bool validate, uint8_t compression_lvl, std::size_t chunk_size,
               std::size_t cache_size)
    : Cooler::Cooler(std::filesystem::path(path_to_file), mode, bin_size, max_str_length,
                     assembly_name, flavor, validate, compression_lvl, chunk_size, cache_size) {}

Cooler::~Cooler() {
  if (this->_mode == WRITE_ONLY && this->_nchroms != 0 && this->_nbins != 0) {
    if (this->_fp) {
      auto &chrom_idx = this->_datasets[IDX_CHR];
      auto &bin1_idx = this->_datasets[IDX_BIN1];
      auto chrom_idx_offset = this->_dataset_file_offsets[IDX_CHR];
      auto bin1_idx_offset = this->_dataset_file_offsets[IDX_BIN1];

      (void)hdf5::write_number(++this->_nbins, chrom_idx, chrom_idx_offset);
      const auto nnz = ++this->_nnz;
      (void)hdf5::write_number(nnz, bin1_idx, bin1_idx_offset);
    } else {
      fmt::print(stderr,
                 FMT_STRING("WARNING: Message for the developers: ~Cooler() for file '{}' was "
                            "called on a file that is supposed to be opened in write-mode, but the "
                            "file handle is actually already closed!"),
                 this->_path_to_file);
    }
  }
}

bool Cooler::is_read_only() const { return this->_mode == Cooler::READ_ONLY; }

std::size_t Cooler::get_nchroms() {
  assert(this->_fp);
  if (this->is_cool()) {
    return static_cast<std::size_t>(hdf5::read_attribute_int(*this->_fp, "nchroms"));
  }
  if (this->is_mcool()) {
    assert(this->_bin_size != 0);
    return static_cast<std::size_t>(hdf5::read_attribute_int(
        *this->_fp, "nchroms", absl::StrCat("/resolutions/", this->_bin_size)));
  }
  utils::throw_with_trace("Unreachable code");
}

void Cooler::get_chr_names(std::vector<std::string> &buff) {
  assert(this->_fp);
  const auto nchroms = this->get_nchroms();
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
  const auto nchroms = this->get_nchroms();
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

bool Cooler::is_cool() const { return this->_flavor == COOL; }
bool Cooler::is_mcool() const { return this->_flavor == MCOOL; }
bool Cooler::is_scool() const { return this->_flavor == SCOOL; }

std::size_t Cooler::get_bin_size() const { return this->_bin_size; }

bool Cooler::has_contacts_for_chr(std::string_view chr_name, bool try_common_chr_prefixes) {
  assert(this->_fp);  // NOLINT
  const auto chr_idx = this->get_chr_idx(chr_name, try_common_chr_prefixes);
  return this->has_contacts_for_chr(chr_idx);
}

bool Cooler::has_contacts_for_chr(std::size_t chr_idx) const {
  assert(this->_fp);                                 // NOLINT
  assert(this->is_read_only());                      // NOLINT
  assert(!this->_idx_bin1_offset.empty());           // NOLINT
  assert(!this->_idx_chrom_offset.empty());          // NOLINT
  assert(chr_idx < this->_idx_chrom_offset.size());  // NOLINT

  const auto first_bin = static_cast<std::size_t>(this->_idx_chrom_offset[chr_idx]);
  const auto last_bin = static_cast<std::size_t>(this->_idx_chrom_offset[chr_idx + 1]);
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
  H5::StrType METADATA_STR_TYPE(H5::PredType::C_S1, H5T_VARIABLE);
  METADATA_STR_TYPE.setCset(H5T_CSET_UTF8);

  H5::DataSpace attr_space(H5S_SCALAR);
  int64_t int_buff{};
  std::string str_buff{};

  auto att = this->_fp->createAttribute("format", METADATA_STR_TYPE, attr_space);
  str_buff = "HDF5::Cooler";
  att.write(METADATA_STR_TYPE, str_buff);

  att = this->_fp->createAttribute("format-version", hdf5::getH5_type<decltype(int_buff)>(),
                                   attr_space);
  int_buff = 3;
  att.write(H5::PredType::NATIVE_INT64, &int_buff);

  att = this->_fp->createAttribute("bin-type", METADATA_STR_TYPE, attr_space);
  str_buff = "fixed";
  att.write(METADATA_STR_TYPE, str_buff);

  att = this->_fp->createAttribute("bin-size", hdf5::getH5_type<decltype(int_buff)>(), attr_space);
  int_buff = static_cast<int64_t>(this->_bin_size);
  att.write(H5::PredType::NATIVE_INT64, &int_buff);

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

template <typename I1, typename I2>
void Cooler::write_or_append_cmatrix_to_file(const ContactMatrix<I1> &cmatrix,
                                             std::string_view chr_name, I2 chr_start, I2 chr_end,
                                             I2 chr_length, bool quiet) {
  // Casting away the constness from the ptr is needed in order to make things work with slices.
  // This shouldn't cause any problem, as in the end the contact matrix is accesses through a
  // const slice NOLINTNEXTLINE
  const std::vector<ContactMatrix<I1> *> cmatrices{const_cast<ContactMatrix<I1> *>(&cmatrix)};
  std::string chr_name_{chr_name.data(), chr_name.size()};
  write_or_append_cmatrices_to_file(
      absl::MakeConstSpan(cmatrices), absl::MakeConstSpan(&chr_name_, 1),
      absl::MakeConstSpan(&chr_start, 1), absl::MakeConstSpan(&chr_end, 1),
      absl::MakeConstSpan(&chr_length, 1), quiet);
}

template <typename I1, typename I2>
void Cooler::write_or_append_cmatrices_to_file(const std::vector<ContactMatrix<I1>> &cmatrices,
                                               const std::vector<std::string> &chr_names,
                                               const std::vector<I2> &chr_starts,
                                               const std::vector<I2> &chr_ends,
                                               const std::vector<I2> &chr_sizes, bool quiet) {
  std::vector<ContactMatrix<I1> *> v(cmatrices.size());
  // Convert references to ptrs
  std::transform(cmatrices.begin(), cmatrices.end(), v.begin(), [](const auto &m) { return &m; });
  write_or_append_cmatrices_to_file(v, chr_names, chr_starts, chr_ends, chr_sizes, quiet);
}

template <typename I1, typename I2>
void Cooler::write_or_append_cmatrices_to_file(const std::vector<ContactMatrix<I1> *> &cmatrices,
                                               const std::vector<std::string> &chr_names,
                                               const std::vector<I2> &chr_starts,
                                               const std::vector<I2> &chr_ends,
                                               const std::vector<I2> &chr_sizes, bool quiet) {
  Cooler::write_or_append_cmatrices_to_file(
      absl::MakeConstSpan(cmatrices), absl::MakeConstSpan(chr_names),
      absl::MakeConstSpan(chr_starts), absl::MakeConstSpan(chr_ends),
      absl::MakeConstSpan(chr_sizes), quiet);
}

template <typename I1, typename I2>
void Cooler::write_or_append_cmatrices_to_file(absl::Span<ContactMatrix<I1> *const> cmatrices,
                                               absl::Span<const std::string> chr_names,
                                               absl::Span<const I2> chr_starts,
                                               absl::Span<const I2> chr_ends,
                                               absl::Span<const I2> chr_sizes, bool quiet) {
  static_assert(std::is_integral_v<I1>, "I1 should be an integral type.");
  static_assert(std::is_integral_v<I1>, "I2 should be an integral type.");

  assert(this->_bin_size != 0);

  this->_nchroms += static_cast<int64_t>(cmatrices.size());
  if (const auto n = chr_names.size();
      n != chr_starts.size() || n != chr_ends.size() || n != chr_sizes.size()) {
    utils::throw_with_trace(std::runtime_error(fmt::format(
        FMT_STRING("Message for the DEVS: The vectors passed to function "
                   "cooler::write_modle_cmatrices_to_cooler() should all have "
                   "the same size!\n  cmatrices.size()={}\n  chr_names.size()={}\n  "
                   "chr_starts={}\n  chr_ends={}\n  chr_sizes={}\n"),
        cmatrices.size(), chr_names.size(), chr_starts.size(), chr_ends.size(), chr_sizes.size())));
  }

  auto t0 = absl::Now();

  // Declare one buffer for each dataset
  constexpr std::size_t BUFF_SIZE = 1024 * 1024 / sizeof(int64_t);  // 1 MB
  std::vector<int64_t> bin_pos_buff;
  std::vector<int32_t> bin_chrom_buff;
  std::vector<int64_t> pixel_b1_idx_buff;
  std::vector<int64_t> pixel_b2_idx_buff;
  std::vector<int64_t> pixel_count_buff;
  std::vector<int64_t> idx_bin1_offset_buff;
  std::vector<int64_t> idx_chrom_offset_buff;

  // Reserve memory for max buff size
  bin_pos_buff.reserve(BUFF_SIZE);
  bin_chrom_buff.reserve(BUFF_SIZE);
  pixel_b1_idx_buff.reserve(BUFF_SIZE);
  pixel_b2_idx_buff.reserve(BUFF_SIZE);
  pixel_count_buff.reserve(BUFF_SIZE);
  idx_bin1_offset_buff.reserve(BUFF_SIZE);
  idx_chrom_offset_buff.reserve(cmatrices.size() + 1);

  // Declare aliases to file offsets for HDF5 datasets
  auto &chr_name_h5_foffset = this->_dataset_file_offsets[CHR_NAME];
  auto &chr_length_h5_foffset = this->_dataset_file_offsets[CHR_LEN];
  auto &bin_chr_name_h5_foffset = this->_dataset_file_offsets[BIN_CHROM];
  auto &bin_start_h5_foffset = this->_dataset_file_offsets[BIN_START];
  auto &bin_end_h5_foffset = this->_dataset_file_offsets[BIN_END];
  auto &pixel_b1_id_h5_foffset = this->_dataset_file_offsets[PXL_B1];
  auto &pixel_b2_id_h5_foffset = this->_dataset_file_offsets[PXL_B2];
  auto &pixel_count_h5_foffset = this->_dataset_file_offsets[PXL_COUNT];
  auto &idx_bin1_offset_h5_foffset = this->_dataset_file_offsets[IDX_BIN1];
  auto &idx_chr_offset_h5_foffset = this->_dataset_file_offsets[IDX_CHR];

  auto &d = this->_datasets;

  // Number of non-zero (values)
  this->_nnz = hdf5::has_attribute(*this->_fp, "nnz", this->_root_path)
                   ? hdf5::read_attribute_int(*this->_fp, "nnz", this->_root_path)
                   : 0;
  // Idx of the first non-zero pixel of the current chr
  auto pxl_offset =
      static_cast<hsize_t>(hdf5::has_attribute(*this->_fp, "nbins", this->_root_path)
                               ? hdf5::read_attribute_int(*this->_fp, "nbins", this->_root_path) - 1
                               : 0);
  auto chr_offset =
      static_cast<hsize_t>(hdf5::has_attribute(*this->_fp, "nchroms", this->_root_path)
                               ? hdf5::read_attribute_int(*this->_fp, "nchroms", this->_root_path)
                               : 0);

  auto write_pixels_to_file =
      [&]() {  // Lambda used to write pixel data to file. Mostly useful to reduce code bloat
        pixel_b1_id_h5_foffset =
            hdf5::write_numbers(pixel_b1_idx_buff, d[PXL_B1], pixel_b1_id_h5_foffset);
        pixel_b2_id_h5_foffset =
            hdf5::write_numbers(pixel_b2_idx_buff, d[PXL_B2], pixel_b2_id_h5_foffset);
        pixel_count_h5_foffset =
            hdf5::write_numbers(pixel_count_buff, d[PXL_COUNT], pixel_count_h5_foffset);

        pixel_b1_idx_buff.clear();
        pixel_b2_idx_buff.clear();
        pixel_count_buff.clear();
      };

  for (auto chr_idx = 0UL; chr_offset + chr_idx < static_cast<std::size_t>(this->_nchroms);
       ++chr_idx) {
    // Declare several aliases/variables to improve code readability in later sections
    const auto &cmatrix = cmatrices[chr_idx];
    const auto &chr_name = chr_names[chr_idx];
    const auto chr_total_len = static_cast<int64_t>(chr_sizes[chr_idx]);
    const auto chr_start = static_cast<uint64_t>(chr_starts[chr_idx]);
    const auto chr_end = static_cast<uint64_t>(chr_ends[chr_idx]);
    const auto chr_simulated_len = chr_end - chr_start;

    if (!quiet) {
      fmt::print(stderr, FMT_STRING("Writing contacts for '{}' ({:.2f} Mbp)..."), chr_name,
                 static_cast<double>(chr_simulated_len) / 1.0e6);  // NOLINT
      if (const auto n = cmatrix->get_n_of_missed_updates(); n != 0) {
        fmt::print(
            stderr, FMT_STRING(" WARNING: Detected {} missed updates ({:.4f}% of the total)."), n,
            static_cast<double>(100U * n) / static_cast<double>(cmatrix->get_tot_contacts()));
      }
    }
    const auto t1 = absl::Now();

    // Write chr name and size
    chr_name_h5_foffset = hdf5::write_str(chr_name, d[CHR_NAME], STR_TYPE, chr_name_h5_foffset);
    chr_length_h5_foffset = hdf5::write_number(chr_total_len, d[CHR_LEN], chr_length_h5_foffset);
    // Add idx of the first bin belonging to the chromosome that is being processed
    idx_chrom_offset_buff.emplace_back(this->_nbins);

    // Write all fixed-size bins for the current chromosome. Maybe in the future we can switch to
    // use variable-size bins
    this->_nbins = this->write_bins(chr_offset + chr_idx, chr_total_len, this->_bin_size,
                                    bin_chrom_buff, bin_pos_buff, this->_nbins);
    // Update file offsets of bin_* datasets
    bin_chr_name_h5_foffset = this->_nbins;
    bin_start_h5_foffset = this->_nbins;
    bin_end_h5_foffset = this->_nbins;

    // In case we are simulating a subset of a chromosome (i.e. end - start != chr size), write
    // the index for all the bins corresponding to genomic coordinates before the start position
    // Example: suppose we are writing contacts for a chromosome "C" of size 10 Mbp. Suppose we
    // also know that there are no contacts in the first and last 2 Mbps. In this case chr_start
    // will be 2 Mbp and chr_end will be 8 Mbp. This for loop write the index for all the bins
    // corresponding to the genomic region 0-2 Mbp. A more elegant solution would be to use
    // variable bin-size and write a single 2 Mbp bin.
    idx_bin1_offset_buff.resize(idx_bin1_offset_buff.size() + (chr_start / this->_bin_size),
                                this->_nnz);
    idx_bin1_offset_h5_foffset =
        hdf5::write_numbers(idx_bin1_offset_buff, d[IDX_BIN1], idx_bin1_offset_h5_foffset);
    idx_bin1_offset_buff.clear();

    if (cmatrix->ncols() == 0) {
      DISABLE_WARNING_PUSH
      DISABLE_WARNING_SIGN_CONVERSION
      idx_bin1_offset_buff.resize(
          (chr_total_len / this->_bin_size) + (chr_total_len % this->_bin_size != 0), this->_nnz);
      DISABLE_WARNING_POP
      idx_bin1_offset_h5_foffset =
          hdf5::write_numbers(idx_bin1_offset_buff, d[IDX_BIN1], idx_bin1_offset_h5_foffset);
      idx_bin1_offset_buff.clear();
    }

    pxl_offset += chr_start / this->_bin_size;
    for (auto i = 0UL; i < cmatrix->ncols(); ++i) {  // Iterate over columns in the cmatrix
      // Write first pixel that refers to a given bin1 to the index
      idx_bin1_offset_buff.push_back(this->_nnz);

      // Iterate over rows of the cmatrix. The first condition serves the purpose to avoid reading
      // data from regions that are full of zeros by design (because cmatrix only stores contacts
      // for a certain width along the diagonal). The second condition makes sure we are not
      // reading past the array end
      for (auto j = i; j < i + cmatrix->nrows() && j < cmatrix->ncols(); ++j) {
        if (const auto m = cmatrix->get(i, j); m != 0) {  // Only write nnz pixels
#ifndef NDEBUG
          // Make sure we are always reading from the upper-triangle of the underlying square
          // contact matrix
          if (pxl_offset + i > pxl_offset + j) {
            utils::throw_with_trace(std::runtime_error(
                fmt::format(FMT_STRING("Cooler::write_or_append_cmatrix_to_file(): b1 > b2: b1={}; "
                                       "b2={}; offset={}; m={}\n"),
                            pixel_b1_idx_buff.back(), pixel_b2_idx_buff.back(), pxl_offset, m)));
          }
#endif
          pixel_b1_idx_buff.push_back(static_cast<int64_t>(pxl_offset + i));
          pixel_b2_idx_buff.push_back(static_cast<int64_t>(pxl_offset + j));
          pixel_count_buff.push_back(m);
          ++this->_nnz;

          if (pixel_b1_idx_buff.size() == BUFF_SIZE) {  // Write pixels to disk when buffer is full
            write_pixels_to_file();
          }
        }
      }

      if (idx_bin1_offset_buff.size() == BUFF_SIZE) {  // Write bin idx when buffer is full
        idx_bin1_offset_h5_foffset =
            hdf5::write_numbers(idx_bin1_offset_buff, d[IDX_BIN1], idx_bin1_offset_h5_foffset);
        idx_bin1_offset_buff.clear();
      }
    }

    // In case we are simulating a subset of a chromosome (i.e. end - start != chr size), write
    // the index for all the bins corresponding to genomic coordinates after the end position. See
    // previous comment for an example
    idx_bin1_offset_buff.resize(
        idx_bin1_offset_buff.size() +
            ((static_cast<std::size_t>(chr_total_len) - chr_end) / this->_bin_size) +
            ((static_cast<std::size_t>(chr_total_len) - chr_end) % this->_bin_size != 0),
        this->_nnz);
    idx_bin1_offset_h5_foffset =
        hdf5::write_numbers(idx_bin1_offset_buff, d[IDX_BIN1], idx_bin1_offset_h5_foffset);
    idx_bin1_offset_buff.clear();

    pxl_offset = this->_nbins;
    if (!quiet) {
      fmt::print(stderr, FMT_STRING(" DONE in {}!\n"), absl::FormatDuration(absl::Now() - t1));
    }
  }

  // Write all non-empty buffers to disk
  if (!pixel_b1_idx_buff.empty()) {
    write_pixels_to_file();
  }

  // Writing the chr_index only at the end should be ok, given that most of the time we are
  // processing 20-100 chr
  idx_chr_offset_h5_foffset =
      hdf5::write_numbers(idx_chrom_offset_buff, d[IDX_CHR], idx_chr_offset_h5_foffset);
  idx_bin1_offset_h5_foffset =
      hdf5::write_numbers(idx_bin1_offset_buff, d[IDX_BIN1], idx_bin1_offset_h5_foffset);

  auto &f = *this->_fp;
  hdf5::write_or_create_attribute(f, "nchroms", this->_nchroms);
  const auto nbins = static_cast<int64_t>(this->_nbins);
  hdf5::write_or_create_attribute(f, "nbins", nbins);
  hdf5::write_or_create_attribute(f, "nnz", this->_nnz);

  if (!quiet) {
    fmt::print(stderr, "All contacts have been written to file {}. Saved {:.2f}M pixels in {}.\n",
               this->_path_to_file, static_cast<double>(this->_nnz) / 1.0e6,
               absl::FormatDuration(absl::Now() - t0));
  }
}

ContactMatrix<uint32_t> Cooler::cooler_to_cmatrix(
    std::string_view chr_name, std::size_t nrows,
    std::pair<std::size_t, std::size_t> chr_boundaries, bool try_common_chr_prefixes,
    bool prefer_using_balanced_counts) {
  assert(this->_fp);
  assert(!this->_datasets.empty());
  assert(!this->_idx_bin1_offset.empty());
  assert(!this->_idx_chrom_offset.empty());

  const auto chr_idx = this->get_chr_idx(chr_name, try_common_chr_prefixes);
  const auto chr_size = static_cast<std::size_t>(this->get_chr_sizes()[chr_idx]);
  if (chr_boundaries.second > chr_size) {
    chr_boundaries.second = chr_size;
  }
  const auto bin_range = std::make_pair(static_cast<std::size_t>(this->_idx_chrom_offset[chr_idx]) +
                                            (chr_boundaries.first / this->_bin_size),
                                        static_cast<std::size_t>(this->_idx_chrom_offset[chr_idx]) +
                                            (chr_boundaries.second / this->_bin_size));
  const auto bin1_offset_idx = this->get_bin1_offset_idx_for_chr(chr_idx, chr_boundaries);
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
      pxl_count_scaling_factor = buff[chr_idx];
    } else {  // standard or --trans-only balancing produces a single scale factor
      hdf5::read_attribute(d, "scale", pxl_count_scaling_factor);
    }
  }

  return cooler_to_cmatrix(bin_range, bin1_offset_idx, nrows, pxl_count_scaling_factor,
                           prefer_using_balanced_counts);
}

ContactMatrix<uint32_t> Cooler::cooler_to_cmatrix(
    std::string_view chr_name, std::size_t diagonal_width, std::size_t bin_size,
    std::pair<std::size_t, std::size_t> chr_boundaries, bool try_common_chr_prefixes,
    bool prefer_using_balanced_counts) {
  assert(this->_bin_size != 0);
  if (bin_size != 0 && this->_bin_size != bin_size) {
    throw std::runtime_error(fmt::format(
        FMT_STRING(
            "Unable to read a Cooler file with bin size {} in a contact matrix with bin size {}"),
        this->_bin_size, bin_size));
  }
  const auto nrows = (diagonal_width / this->_bin_size) +
                     static_cast<std::size_t>(diagonal_width % this->_bin_size != 0);
  return cooler_to_cmatrix(chr_name, nrows, chr_boundaries, try_common_chr_prefixes,
                           prefer_using_balanced_counts);
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
  (void)fill_value;
  (void)type;

  H5::DSetCreatPropList prop{};
  prop.setChunk(1, &chunk_size);
  prop.setDeflate(compression_lvl);
  if constexpr (!std::is_constructible_v<H5std_string, T2>) {
    prop.setFillValue(type, &fill_value);
  }

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
    prop.setChunkCache(default_multiplier * (cache_size / chunk_size), cache_size, 0.01);
    (void)type;
  } else {
    if (type == H5::PredType::NATIVE_INT64) {  // int64_t
      prop.setChunkCache(default_multiplier * (cache_size / chunk_size), cache_size, rdcc_w0);
    } else if (type == H5::PredType::NATIVE_INT32) {  // int32_t
      prop.setChunkCache(default_multiplier * (cache_size / (chunk_size + chunk_size)), cache_size,
                         rdcc_w0);
    } else if (type == H5::PredType::NATIVE_DOUBLE) {
      prop.setChunkCache(default_multiplier * (cache_size / chunk_size), cache_size, rdcc_w0);
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
    throw std::runtime_error(hdf5::construct_error_stack());
  }
}

std::vector<H5::Group> Cooler::open_groups(H5::H5File &f, bool create_if_not_exist,
                                           std::size_t bin_size) {
  std::vector<H5::Group> groups(4);  // NOLINT

  auto open_or_create_group = [&](Cooler::Groups g, const std::string &group_name) {
    groups[g] = create_if_not_exist && !f.nameExists(group_name) ? f.createGroup(group_name)
                                                                 : f.openGroup(group_name);
  };
  const std::string root_path = bin_size != 0 && hdf5::has_group(f, "/resolutions")
                                    ? absl::StrCat("/resolutions/", bin_size, "/")
                                    : "";
  try {
    open_or_create_group(CHR, absl::StrCat(root_path, "chroms"));
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
    dset[CHR_LEN] =
        f.createDataSet(absl::StrCat(r, "/chroms/length"), INT64, mem_space, cp64, ap64);
    dset[CHR_NAME] = f.createDataSet(absl::StrCat(r, "/chroms/name"), STR, mem_space, cps, aps);

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
    d[CHR_LEN] = f.openDataSet(absl::StrCat(this->_root_path, "chroms/length"), ai64);
    d[CHR_NAME] = f.openDataSet(absl::StrCat(this->_root_path, "chroms/name"), as);

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

template <typename I1, typename I2, typename I3>
hsize_t Cooler::write_bins(I1 chrom_, I2 length_, I3 bin_size_, std::vector<int32_t> &buff32,
                           std::vector<int64_t> &buff64, hsize_t file_offset, hsize_t buff_size) {
  assert(bin_size_ != 0);
  buff32.resize(buff_size);
  buff64.resize(buff_size);
  const auto chrom = static_cast<int32_t>(chrom_);
  const auto length = static_cast<int64_t>(length_);
  const auto bin_size = static_cast<int64_t>(bin_size_);

  std::fill(buff32.begin(), buff32.end(), chrom);

  int64_t start = 0;
  int64_t end = std::min(length, bin_size);
  hsize_t bins_processed = 0;
  const int64_t nbins = (length / bin_size) + (length % bin_size != 0);
  const int64_t nchunks =
      (nbins / static_cast<int64_t>(buff_size)) + (nbins % static_cast<int64_t>(buff_size) != 0);
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

absl::Span<const int64_t> Cooler::get_bin1_offset_idx_for_chr(
    std::size_t chr_idx, std::pair<std::size_t, std::size_t> chr_subrange) {
  assert(!this->_idx_bin1_offset.empty());           // NOLINT
  assert(!this->_idx_chrom_offset.empty());          // NOLINT
  assert(chr_idx < this->_idx_chrom_offset.size());  // NOLINT
  const auto chr_start_bin = static_cast<std::size_t>(this->_idx_chrom_offset[chr_idx]) +
                             (chr_subrange.first / this->_bin_size);
  const auto chr_end_bin =
      chr_subrange.second == std::numeric_limits<decltype(chr_subrange.second)>::max()
          ? static_cast<std::size_t>(this->_idx_chrom_offset[chr_idx + 1])
          : static_cast<std::size_t>(this->_idx_chrom_offset[chr_idx]) +
                (chr_subrange.second / this->_bin_size);
  DISABLE_WARNING_PUSH
  DISABLE_WARNING_SIGN_COMPARE
  assert(chr_end_bin <= this->_idx_chrom_offset[chr_idx + 1]);
  assert(chr_end_bin >= chr_start_bin);  // NOLINT
  DISABLE_WARNING_POP

  return absl::MakeConstSpan(this->_idx_bin1_offset)
      .subspan(chr_start_bin, chr_end_bin - chr_start_bin);
}

absl::Span<const int64_t> Cooler::get_bin1_offset_idx_for_chr(
    std::string_view chr_name, std::pair<std::size_t, std::size_t> chr_subrange) {
  const auto chr_idx = get_chr_idx(chr_name);
  return get_bin1_offset_idx_for_chr(chr_idx, chr_subrange);
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

ContactMatrix<uint32_t> Cooler::cooler_to_cmatrix(std::pair<hsize_t, hsize_t> bin_range,
                                                  absl::Span<const int64_t> bin1_offset_idx,
                                                  std::size_t nrows, double bias_scaling_factor,
                                                  bool prefer_using_balanced_counts) {
  if (this->_datasets.empty()) {
    this->open_default_datasets();
  }

  const auto &[first_bin, last_bin] = bin_range;
  assert(first_bin < last_bin);
  assert(last_bin <= first_bin + bin1_offset_idx.size());
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
        std::min(static_cast<std::size_t>(bin1_offset_idx[i] - bin1_offset_idx[i - 1]), nrows);
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
                                                  std::size_t nrows, double scaling_factor,
                                                  bool prefer_using_balanced_counts) {
  return this->cooler_to_cmatrix(bin_range, absl::MakeConstSpan(bin1_offset_idx), nrows,
                                 scaling_factor, prefer_using_balanced_counts);
}

std::size_t Cooler::get_chr_idx(std::string_view query_chr_name, bool try_common_chr_prefixes) {
  // Here's the issue: for a given genome assembly (say hg19), some tools name chromosome as
  // chrN, where N is the chromosome number, while others only use the number N. The purpose
  // of this lambda is to try to guess few reasonably common prefixes, look them up in the
  // Cooler file, then return the chromosome name variant that produced a hit together with
  // the chr index

  try {
    const auto chr_names = hdf5::read_strings(this->_datasets[CHR_NAME], 0);
    auto match = std::find(chr_names.begin(), chr_names.end(), query_chr_name);
    if (match != chr_names.end()) {
      return static_cast<std::size_t>(std::distance(chr_names.begin(), match));
    }
    if (try_common_chr_prefixes) {
      std::array<std::string_view, 3> queries;
      constexpr std::array<std::string_view, 3> prefixes = {"chr", "CHR", "Chr"};
      for (auto i = 0U; i < prefixes.size(); ++i) {
        queries[i] = absl::StartsWith(query_chr_name, prefixes[i])
                         ? absl::StripPrefix(query_chr_name, prefixes[i])
                         : absl::StrCat(prefixes[i], query_chr_name);
      }
      for (const auto &q : queries) {
        match = std::find(chr_names.begin(), chr_names.end(), q);
        if (match != chr_names.end()) {
          return static_cast<std::size_t>(std::distance(chr_names.begin(), match));
        }
      }
      throw std::runtime_error(
          fmt::format(FMT_STRING("Unable to find a chromosome named '{}'. The following chromosome "
                                 "name variants were searched: '{}'"),
                      query_chr_name, absl::StrJoin(queries, "', '")));
    }
    throw std::runtime_error(
        fmt::format(FMT_STRING("Unable to find a chromosome named '{}'"), query_chr_name));
  } catch (const std::exception &e) {
    if (absl::StartsWith(e.what(), "Unable to find a chromosome")) {
      throw;
    }
    throw std::runtime_error(fmt::format(FMT_STRING("The following error occurred while looking up "
                                                    "'{}' in dataset chroms/name of file '{}': {}"),
                                         query_chr_name, this->_path_to_file, e.what()));
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

bool Cooler::validate_cool_flavor(H5::H5File &f, std::size_t bin_size, std::string_view root_path,
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

bool Cooler::validate_multires_cool_flavor(H5::H5File &f, std::size_t bin_size,
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
