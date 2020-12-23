#pragma once

#include <H5Cpp.h>
#include <absl/time/clock.h>

#include <algorithm>
#include <range/v3/action/remove_if.hpp>
#include <range/v3/action/transform.hpp>
#include <range/v3/algorithm/transform.hpp>
#include <range/v3/view/chunk.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/filter.hpp>
#include <range/v3/view/generate.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/zip.hpp>
#include <string_view>
#include <type_traits>
#include <vector>

#include "modle/contacts.hpp"
#include "modle/suppress_compiler_warnings.hpp"
#include "modle/utils.hpp"

namespace modle::cooler {

template <typename DataType>
H5::PredType getH5_type() {
  // Taken from github.com/DavidAce/h5pp:
  // https://github.com/DavidAce/h5pp/blob/master/h5pp/include/h5pp/details/h5ppUtils.h
  // clang-format off
  using DecayType = typename std::decay<DataType>::type;
  if constexpr (std::is_same_v<DecayType, short>) return H5::PredType::NATIVE_SHORT;               // NOLINT
  if constexpr (std::is_same_v<DecayType, int>) return H5::PredType::NATIVE_INT;                   // NOLINT
  if constexpr (std::is_same_v<DecayType, long>) return H5::PredType::NATIVE_LONG;                 // NOLINT
  if constexpr (std::is_same_v<DecayType, long long>) return H5::PredType::NATIVE_LLONG;           // NOLINT
  if constexpr (std::is_same_v<DecayType, unsigned short>) return H5::PredType::NATIVE_USHORT;     // NOLINT
  if constexpr (std::is_same_v<DecayType, unsigned int>) return H5::PredType::NATIVE_UINT;         // NOLINT
  if constexpr (std::is_same_v<DecayType, unsigned long>) return H5::PredType::NATIVE_ULONG;       // NOLINT
  if constexpr (std::is_same_v<DecayType, unsigned long long>) return H5::PredType::NATIVE_ULLONG; // NOLINT
  if constexpr (std::is_same_v<DecayType, double>) return H5::PredType::NATIVE_DOUBLE;             // NOLINT
  if constexpr (std::is_same_v<DecayType, long double>) return H5::PredType::NATIVE_LDOUBLE;       // NOLINT
  if constexpr (std::is_same_v<DecayType, float>) return H5::PredType::NATIVE_FLOAT;               // NOLINT
  if constexpr (std::is_same_v<DecayType, bool>) return H5::PredType::NATIVE_UINT8;                // NOLINT
  if constexpr (std::is_same_v<DecayType, std::string>) return H5::PredType::C_S1;                 // NOLINT
  if constexpr (std::is_same_v<DecayType, char>) return H5::PredType::C_S1;                        // NOLINT
  // clang-format on

  throw std::logic_error("getH5_type(): Unable to map C++ type to a H5T.");
}

template <typename S>
hsize_t write_str(const S &str, H5::H5File &f, std::string_view dataset_name,
                  hsize_t MAX_STR_LENGTH, hsize_t file_offset, hsize_t CHUNK_DIMS,
                  uint8_t COMPRESSION_LEVEL) {
  static_assert(std::is_convertible<S, H5std_string>::value,
                "S should be a type convertible to std::string.");

  const H5std_string d{dataset_name};  // Probably unnecessary
  constexpr const hsize_t RANK{1};     // i.e. number of dimensions
  constexpr const hsize_t BUFF_SIZE{1};
  constexpr const hsize_t MAXDIMS{H5S_UNLIMITED};  // extensible dataset
  const H5std_string FILL_VAL_STR{'\0'};

  try {
    auto mem_space{H5::DataSpace(RANK, &BUFF_SIZE, &MAXDIMS)};

    // Create the datatype as follows
    H5::StrType STD_STR(H5::PredType::C_S1, MAX_STR_LENGTH);
    STD_STR.setStrpad(H5T_STR_NULLPAD);  // This seems to be the most robust way of terminating
    // strings. Using NULLPAD sometimes causes some garbage to
    // appear at end of string
    STD_STR.setCset(H5T_CSET_ASCII);

    H5::DSetCreatPropList cparms{};
    cparms.setChunk(RANK, &CHUNK_DIMS);
    cparms.setFillValue(STD_STR, &FILL_VAL_STR);
    cparms.setDeflate(COMPRESSION_LEVEL);

    auto dataset =
        f.nameExists(d) ? f.openDataSet(d) : f.createDataSet(d, STD_STR, mem_space, cparms);

    H5::DataSpace file_space;  // File space

    // This is probably not very efficient, but it should be fine, given that we don't expect to
    // write that many strings
    hsize_t file_size{file_offset + 1};
    dataset.extend(&file_size);
    file_space = dataset.getSpace();
    file_space.selectHyperslab(H5S_SELECT_SET, &BUFF_SIZE, &file_offset);
    dataset.write(std::string{str}, STD_STR, mem_space, file_space);

    return file_size;
  } catch (const H5::FileIException &e) {
    throw std::runtime_error(e.getDetailMsg());
  } catch (const H5::GroupIException &e) {
    throw std::runtime_error(e.getDetailMsg());
  } catch (const H5::DataSetIException &e) {
    throw std::runtime_error(e.getDetailMsg());
  } catch (const H5::DataSpaceIException &e) {
    throw std::runtime_error(e.getDetailMsg());
  }
}

template <typename S>
hsize_t write_vect_of_str(std::vector<S> &data, H5::H5File &f, std::string_view dataset_name,
                          hsize_t file_offset, hsize_t CHUNK_DIMS, uint8_t COMPRESSION_LEVEL) {
  static_assert(std::is_convertible<S, H5std_string>::value,
                "S should be a type convertible to std::string.");

  const H5std_string d{dataset_name};  // Probably unnecessary
  constexpr const hsize_t RANK{1};     // i.e. number of dimensions
  constexpr const hsize_t BUFF_SIZE{1};
  constexpr const hsize_t MAXDIMS{H5S_UNLIMITED};  // extensible dataset
  const H5std_string FILL_VAL_STR{'\0'};

  try {
    auto mem_space{H5::DataSpace(RANK, &BUFF_SIZE, &MAXDIMS)};

    const hsize_t MAX_STR_LENGTH =
        std::max_element(data.begin(), data.end(),
                         [](const auto &s1, const auto &s2) { return s1.size() < s2.size(); })
            ->size() +
        1;
    // Create the datatype as follows
    H5::StrType STD_STR(H5::PredType::C_S1, MAX_STR_LENGTH);
    STD_STR.setStrpad(H5T_STR_NULLTERM);  // This seems to be the most robust way of terminating
                                          // strings. Using NULLPAD sometimes causes some garbage to
                                          // appear at end of string
    STD_STR.setCset(H5T_CSET_ASCII);

    H5::DSetCreatPropList cparms{};
    cparms.setChunk(RANK, &CHUNK_DIMS);
    cparms.setFillValue(STD_STR, &FILL_VAL_STR);
    cparms.setDeflate(COMPRESSION_LEVEL);

    auto dataset =
        f.nameExists(d) ? f.openDataSet(d) : f.createDataSet(d, STD_STR, mem_space, cparms);

    H5::DataSpace file_space;  // File space

    // This is probably not very efficient, but it should be fine, given that we don't expect to
    // write that many strings
    hsize_t file_size{file_offset + data.size()};
    dataset.extend(&file_size);
    for (const auto &s : data) {
      file_space = dataset.getSpace();
      file_space.selectHyperslab(H5S_SELECT_SET, &BUFF_SIZE, &file_offset);
      dataset.write(s.c_str(), STD_STR, mem_space, file_space);
      ++file_offset;
    }
    return file_offset;
  } catch (const H5::FileIException &e) {
    throw std::runtime_error(e.getDetailMsg());
  } catch (const H5::GroupIException &e) {
    throw std::runtime_error(e.getDetailMsg());
  } catch (const H5::DataSetIException &e) {
    throw std::runtime_error(e.getDetailMsg());
  } catch (const H5::DataSpaceIException &e) {
    throw std::runtime_error(e.getDetailMsg());
  }
}

template <typename I>
hsize_t write_int(I num, H5::H5File &f, std::string_view dataset_name, hsize_t file_offset,
                  hsize_t CHUNK_DIMS, uint8_t COMPRESSION_LEVEL) {
  static_assert(std::is_integral<I>::value, "I should be an integer type.");

  const H5std_string d{dataset_name};  // Probably unnecessary
  constexpr const hsize_t RANK{1};     // i.e. number of dimensions
  const hsize_t BUFF_SIZE{1};
  constexpr const hsize_t MAXDIMS{H5S_UNLIMITED};  // extensible dataset
  constexpr const I FILL_VAL{0};

  try {
    const auto mem_space{H5::DataSpace(RANK, &BUFF_SIZE, &MAXDIMS)};

    H5::DSetCreatPropList cparms{};
    cparms.setChunk(RANK, &CHUNK_DIMS);
    cparms.setFillValue(getH5_type<I>(), &FILL_VAL);
    cparms.setDeflate(COMPRESSION_LEVEL);

    auto dataset =
        f.nameExists(d) ? f.openDataSet(d) : f.createDataSet(d, getH5_type<I>(), mem_space, cparms);

    H5::DataSpace file_space;

    hsize_t file_size{file_offset + BUFF_SIZE};
    dataset.extend(&file_size);
    file_space = dataset.getSpace();
    file_space.selectHyperslab(H5S_SELECT_SET, &BUFF_SIZE, &file_offset);
    dataset.write(&num, getH5_type<I>(), mem_space, file_space);
    return file_size;

  } catch (const H5::FileIException &e) {
    throw std::runtime_error(e.getDetailMsg());
  } catch (const H5::GroupIException &e) {
    throw std::runtime_error(e.getDetailMsg());
  } catch (const H5::DataSetIException &e) {
    throw std::runtime_error(e.getDetailMsg());
  } catch (const H5::DataSpaceIException &e) {
    throw std::runtime_error(e.getDetailMsg());
  }
}

template <typename I>
hsize_t write_vect_of_int(std::vector<I> &data, H5::H5File &f, std::string_view dataset_name,
                          hsize_t file_offset, hsize_t CHUNK_DIMS, uint8_t COMPRESSION_LEVEL) {
  static_assert(std::is_integral<I>::value, "I should be an integer type.");

  const H5std_string d{dataset_name};  // Probably unnecessary
  constexpr const hsize_t RANK{1};     // i.e. number of dimensions
  const hsize_t BUFF_SIZE{data.size()};
  constexpr const hsize_t MAXDIMS{H5S_UNLIMITED};  // extensible dataset
  constexpr const I FILL_VAL{0};

  try {
    const auto mem_space{H5::DataSpace(RANK, &BUFF_SIZE, &MAXDIMS)};

    H5::DSetCreatPropList cparms{};
    cparms.setChunk(RANK, &CHUNK_DIMS);
    cparms.setFillValue(getH5_type<I>(), &FILL_VAL);
    cparms.setDeflate(COMPRESSION_LEVEL);

    auto dataset =
        f.nameExists(d) ? f.openDataSet(d) : f.createDataSet(d, getH5_type<I>(), mem_space, cparms);

    H5::DataSpace file_space;

    hsize_t file_size{file_offset + BUFF_SIZE};
    dataset.extend(&file_size);
    DISABLE_WARNING_PUSH
    DISABLE_WARNING_C_VLA
    const auto *cbuff = reinterpret_cast<I(*)[BUFF_SIZE]>(data.data());  // NOLINT
    DISABLE_WARNING_POP
    file_space = dataset.getSpace();
    file_space.selectHyperslab(H5S_SELECT_SET, &BUFF_SIZE, &file_offset);
    dataset.write(cbuff, getH5_type<I>(), mem_space, file_space);
    return file_size;

  } catch (const H5::FileIException &e) {
    throw std::runtime_error(e.getDetailMsg());
  } catch (const H5::GroupIException &e) {
    throw std::runtime_error(e.getDetailMsg());
  } catch (const H5::DataSetIException &e) {
    throw std::runtime_error(e.getDetailMsg());
  } catch (const H5::DataSpaceIException &e) {
    throw std::runtime_error(e.getDetailMsg());
  }
}

template <typename I>
hsize_t read_int(H5::H5File &f, std::string_view dataset_name, I &BUFF, hsize_t file_offset) {
  static_assert(std::is_integral<I>::value, "I should be an integer type.");
  try {
    H5::Exception::dontPrint();

    const H5std_string d{dataset_name};
    constexpr const hsize_t DIMS = 1;
    constexpr const hsize_t RANK = 1;

    auto dataset = f.openDataSet(d);
    const auto DTYPE = dataset.getDataType();
    auto file_space = dataset.getSpace();

    auto mem_space{H5::DataSpace(RANK, &DIMS)};

    file_space = dataset.getSpace();
    file_space.selectHyperslab(H5S_SELECT_SET, &RANK, &file_offset);
    dataset.read(&BUFF, DTYPE, mem_space, file_space);

    return file_offset + 1;

  } catch (const H5::Exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("Failed to read data into an integer variable (dataset = {}; file offset = "
                   "{}): Function '{}' failed with error '{}'"),
        dataset_name, file_offset, e.getFuncName(), e.getDetailMsg()));
  }
}

template <typename I>
hsize_t read_vect_of_int(H5::H5File &f, std::string_view dataset_name, std::vector<I> &BUFF,
                         hsize_t file_offset) {
  static_assert(std::is_integral<I>::value, "I should be an integer type.");

  const H5std_string d{dataset_name};  // Probably unnecessary
  constexpr const hsize_t RANK{1};     // i.e. number of dimensions
  hsize_t BUFF_SIZE{BUFF.size()};

  try {
    auto dataset = f.openDataSet(d);
    auto file_space = dataset.getSpace();

    if (const auto n_items = static_cast<hsize_t>(file_space.getSimpleExtentNpoints());
        n_items < BUFF_SIZE || BUFF_SIZE == 0) {
      BUFF_SIZE = n_items;
      BUFF.resize(BUFF_SIZE);
    }

    const auto mem_space{H5::DataSpace(RANK, &BUFF_SIZE)};

    DISABLE_WARNING_PUSH
    DISABLE_WARNING_C_VLA
    auto *cbuff = reinterpret_cast<I(*)[BUFF_SIZE]>(BUFF.data());  // NOLINT
    DISABLE_WARNING_POP

    file_space.selectHyperslab(H5S_SELECT_SET, &BUFF_SIZE, &file_offset);
    dataset.read(cbuff, getH5_type<I>(), mem_space, file_space);

    return file_offset + BUFF_SIZE;
  } catch (const H5::Exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("Failed to read data into a vector of integers (dataset = {}; file offset = "
                   "{}): Function '{}' failed with error '{}'"),
        dataset_name, file_offset, e.getFuncName(), e.getDetailMsg()));
  }
}

/*
template <typename I>
void write_modle_cmatrix_to_cooler(const ContactMatrix<I> &cmatrix, std::string_view chr_name,
                                   uint64_t chr_start, uint64_t chr_end, uint64_t chr_length,
                                   std::string_view output_file, bool force_overwrite) {
  const std::vector<ContactMatrix<I>> cmatrices{cmatrix};
  write_modle_cmatrices_to_cooler(cmatrices, std::vector < std::string_view{chr_name},
                                std::vector<uint64_t>{chr_start}, std::vector<uint64_t>{chr_end},
                                std::vector<uint64_t>{chr_length}, output_file, force_overwrite);
}
 */

template <typename I>
void write_modle_cmatrices_to_cooler(const std::vector<ContactMatrix<I>> &cmatrices,
                                     const std::vector<std::string> &chr_names,
                                     const std::vector<uint64_t> &chr_starts,
                                     const std::vector<uint64_t> &chr_ends,
                                     const std::vector<uint64_t> &chr_sizes, uint64_t bin_size,
                                     std::string_view output_file, bool force_overwrite) {
  static_assert(std::is_integral<I>::value, "I should be an integral type.");

  fmt::print(stderr, FMT_STRING("Writing {} contact matrices for to file '{}'...\n"),
             cmatrices.size(), output_file);

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

  if (force_overwrite && std::filesystem::exists(output_file)) {
    std::filesystem::remove(output_file);
  }

  std::filesystem::create_directories(std::filesystem::path(output_file).parent_path());

  auto f = cooler::init_file(output_file, force_overwrite);

  std::size_t chr_name_foffset = 0;
  std::size_t chr_length_foffset = 0;

  constexpr const std::size_t BUFF_SIZE = 1024 * 1024 / sizeof(int64_t);  // 1 MB
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

  auto write_pixels_to_file = [&]() {
    pixel_b1_idx_h5df_offset =
        cooler::write_vect_of_int(pixel_b1_idx_buff, f, "pixels/bin1_id", pixel_b1_idx_h5df_offset);
    pixel_b2_idx_h5df_offset =
        cooler::write_vect_of_int(pixel_b2_idx_buff, f, "pixels/bin2_id", pixel_b2_idx_h5df_offset);
    pixel_count_h5df_offset =
        cooler::write_vect_of_int(pixel_count_buff, f, "pixels/count", pixel_count_h5df_offset);

    pixel_b1_idx_buff.clear();
    pixel_b2_idx_buff.clear();
    pixel_count_buff.clear();
  };

  int64_t nnz = 0;
  int64_t offset = 0;
  int64_t nbins = 0;

  cooler::write_metadata(f, static_cast<int32_t>(bin_size));

  const auto max_chrom_name_length =
      std::max_element(chr_names.begin(), chr_names.end(),
                       [](const auto &s1, const auto &s2) { return s1.size() < s2.size(); })
          ->size() +
      1;
  for (auto chr_idx = 0UL; chr_idx < n_chromosomes; ++chr_idx) {
    const auto &cmatrix = cmatrices[chr_idx];
    const auto &chr_name = chr_names[chr_idx];
    const auto &chr_total_len = chr_sizes[chr_idx];
    const auto &chr_start = chr_starts[chr_idx];
    const auto &chr_end = chr_ends[chr_idx];
    const auto chr_simulated_len = chr_end - chr_start;
    chr_name_foffset =
        cooler::write_str(chr_name, f, "chroms/name", max_chrom_name_length, chr_name_foffset);
    chr_length_foffset = cooler::write_int(chr_total_len, f, "chroms/length", chr_length_foffset);
    DISABLE_WARNING_PUSH
    DISABLE_WARNING_CONVERSION
    DISABLE_WARNING_SIGN_CONVERSION
    fmt::print(stderr, FMT_STRING("Processing '{}' ({:.2f} Mbp)..."), chr_name,
               chr_simulated_len / 1.0e6);  // NOLINT
    const auto t1 = absl::Now();
    idx_chrom_offset_buff.push_back(nbins);
    nbins = cooler::write_bins(f, chr_idx++, chr_simulated_len, bin_size,  // NOLINT
                               bin_chrom_buff, bin_pos_buff, nbins);
    DISABLE_WARNING_POP

    for (auto i = 0UL; i < chr_start / bin_size; ++i) {
      idx_bin1_offset_buff.push_back(nnz);
      if (idx_bin1_offset_buff.size() == BUFF_SIZE) {
        idx_bin1_offset_h5df_offset = cooler::write_vect_of_int(
            idx_bin1_offset_buff, f, "indexes/bin1_offset", idx_bin1_offset_h5df_offset);
        idx_bin1_offset_buff.clear();
      }
    }

    offset += static_cast<int64_t>(chr_start / bin_size);
    for (auto i = 0UL; i < cmatrix.n_cols(); ++i) {
      idx_bin1_offset_buff.push_back(nnz);
      for (auto j = i; j < i + cmatrix.n_rows() && j < cmatrix.n_cols(); ++j) {
        if (const auto m = cmatrix.get(i, j); m != 0) {
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
        idx_bin1_offset_h5df_offset = cooler::write_vect_of_int(
            idx_bin1_offset_buff, f, "indexes/bin1_offset", idx_bin1_offset_h5df_offset);
        idx_bin1_offset_buff.clear();
      }
    }
    for (auto i = 0UL; i < (chr_total_len - chr_end) / bin_size; ++i) {
      idx_bin1_offset_buff.push_back(nnz);
      if (idx_bin1_offset_buff.size() == BUFF_SIZE) {
        idx_bin1_offset_h5df_offset = cooler::write_vect_of_int(
            idx_bin1_offset_buff, f, "indexes/bin1_offset", idx_bin1_offset_h5df_offset);
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

  cooler::write_vect_of_int(idx_chrom_offset_buff, f, "indexes/chrom_offset");
  cooler::write_vect_of_int(idx_bin1_offset_buff, f, "indexes/bin1_offset",
                            idx_bin1_offset_h5df_offset);

  H5::DataSpace att_space(H5S_SCALAR);
  auto att = f.createAttribute("nbins", H5::PredType::NATIVE_INT, att_space);
  att.write(H5::PredType::NATIVE_INT, &nbins);
  att = f.createAttribute("nnz", H5::PredType::NATIVE_INT, att_space);
  att.write(H5::PredType::NATIVE_INT, &nnz);

  fmt::print(stderr, "DONE! Saved {} contacts in {}.\n", nnz,
             absl::FormatDuration(absl::Now() - t0));
}

}  // namespace modle::cooler