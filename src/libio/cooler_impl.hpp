
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
    STD_STR.setStrpad(H5T_STR_NULLTERM);  // This seems to be the most roboust way of terminating
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

hsize_t write_vect_of_enums(std::vector<int32_t> &data, const H5::EnumType &ENUM, H5::H5File &f,
                            std::string_view dataset_name, hsize_t file_offset, hsize_t CHUNK_DIMS,
                            uint8_t COMPRESSION_LEVEL) {
  const H5std_string d{dataset_name};  // Probably unnecessary
  constexpr const hsize_t RANK{1};     // i.e. number of dimensions
  const hsize_t BUFF_SIZE{data.size()};
  constexpr const hsize_t MAXDIMS{H5S_UNLIMITED};  // extensible dataset

  try {
    const auto mem_space{H5::DataSpace(RANK, &BUFF_SIZE, &MAXDIMS)};

    H5::DSetCreatPropList cparms{};
    cparms.setChunk(RANK, &CHUNK_DIMS);
    cparms.setDeflate(COMPRESSION_LEVEL);

    auto dataset = f.nameExists(d) ? f.openDataSet(d) : f.createDataSet(d, ENUM, mem_space, cparms);

    hsize_t file_size{file_offset + BUFF_SIZE};
    dataset.extend(&file_size);
    auto file_space = dataset.getSpace();
    const auto *cbuff = reinterpret_cast<int32_t(*)[BUFF_SIZE]>(data.data());  // NOLINT
    file_space.selectHyperslab(H5S_SELECT_SET, &BUFF_SIZE, &file_offset);
    dataset.write(cbuff, ENUM, mem_space, file_space);
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
    const auto *cbuff = reinterpret_cast<I(*)[BUFF_SIZE]>(data.data());  // NOLINT
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

[[nodiscard]] H5::H5File init_file(std::string_view path_to_output_file, bool force_overwrite) {
  H5::H5File f(std::string{path_to_output_file}, force_overwrite ? H5F_ACC_TRUNC : H5F_ACC_CREAT);

  for (const auto &g : {"chroms", "bins", "pixels", "indexes"}) {
    f.createGroup(g);
  }
  return f;
}

H5::EnumType init_enum_from_strs(const std::vector<std::string> &data, int32_t offset) {
  H5::EnumType ENUM{getH5_type<int32_t>()};
  DISABLE_WARNING_PUSH
  DISABLE_WARNING_SIGN_COMPARE
  DISABLE_WARNING_SIGN_CONVERSION
  for (int32_t i = offset; i < data.size(); ++i) {
    ENUM.insert(data[i], &i);
  }
  DISABLE_WARNING_POP
  return ENUM;
}

H5::EnumType write_chroms(H5::H5File &f,
                          const typename std::vector<ContactMatrix<int64_t>::Header> &headers) {
  std::vector<std::string> names(headers.size());
  std::vector<int64_t> lengths(headers.size());
  for (const auto &[i, h] : headers | ranges::views::enumerate) {
    names[i] = std::string{h.chr_name};
    lengths[i] = static_cast<int64_t>(h.end);
  }
  write_vect_of_int(lengths, f, "chroms/length");
  write_vect_of_str(names, f, "chroms/name");

  H5::DataSpace att_space(H5S_SCALAR);
  auto att = f.createAttribute("nchroms", H5::PredType::NATIVE_INT, att_space);
  const auto nchroms = headers.size();
  att.write(H5::PredType::NATIVE_INT, &nchroms);

  return init_enum_from_strs(names);
}

hsize_t write_bins(H5::H5File &f, H5::EnumType &ENUM, const std::string &chrom, int64_t length,
                   int64_t bin_size, std::vector<int32_t> &buff32, std::vector<int64_t> &buff64,
                   hsize_t file_offset, hsize_t BUFF_SIZE) {
  return write_bins(f, ENUM.getMemberIndex(chrom), length, bin_size, buff32, buff64, file_offset,
                    BUFF_SIZE);
}

hsize_t write_bins(H5::H5File &f, int32_t chrom, int64_t length, int64_t bin_size,
                   std::vector<int32_t> &buff32, std::vector<int64_t> &buff64, hsize_t file_offset,
                   hsize_t BUFF_SIZE) {
  buff32.resize(BUFF_SIZE);
  buff64.resize(BUFF_SIZE);

  std::fill(buff32.begin(), buff32.end(), chrom);

  int64_t start = 0;
  int64_t end = std::min(length, bin_size);
  int64_t bins_processed = 0;
  const int64_t nbins = (length / bin_size) + 1;
  const int64_t nchunks = (nbins / static_cast<int64_t>(BUFF_SIZE)) + 1;

  for (auto i = 0; i < nchunks; ++i) {
    const auto chunk_size = std::min(BUFF_SIZE, static_cast<hsize_t>(nbins - bins_processed));
    if (chunk_size != BUFF_SIZE) {
      buff32.resize(chunk_size);
      buff64.resize(chunk_size);
    }
    write_vect_of_int(buff32, f, "bins/chrom", file_offset);

    std::generate(buff64.begin(), buff64.end(), [&]() { return (start += bin_size) - bin_size; });
    write_vect_of_int(buff64, f, "bins/start", file_offset);

    std::generate(buff64.begin(), buff64.end(), [&]() { return (end += bin_size) - bin_size; });
    if (chunk_size != BUFF_SIZE) {
      buff64.back() = length;
    }
    file_offset = write_vect_of_int(buff64, f, "bins/end", file_offset);
    bins_processed += chunk_size;
  }
  assert(buff64.back() == length);

  return file_offset;
}

void write_contacts(H5::H5File &f, const std::vector<std::string_view> &path_to_cmatrices,
                    H5::EnumType &ENUM) {
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
  idx_chrom_offset_buff.reserve(path_to_cmatrices.size());

  hsize_t pixel_b1_idx_h5df_offset = 0;
  hsize_t pixel_b2_idx_h5df_offset = 0;
  hsize_t pixel_count_h5df_offset = 0;
  hsize_t idx_bin1_offset_h5df_offset = 0;

  auto write_pixels_to_file = [&]() {
    pixel_b1_idx_h5df_offset =
        write_vect_of_int(pixel_b1_idx_buff, f, "pixels/bin1_id", pixel_b1_idx_h5df_offset);
    pixel_b2_idx_h5df_offset =
        write_vect_of_int(pixel_b2_idx_buff, f, "pixels/bin2_id", pixel_b2_idx_h5df_offset);
    pixel_count_h5df_offset =
        write_vect_of_int(pixel_count_buff, f, "pixels/count", pixel_count_h5df_offset);

    pixel_b1_idx_buff.clear();
    pixel_b2_idx_buff.clear();
    pixel_count_buff.clear();
  };

  int64_t n = 0;
  int64_t offset = 0;
  int32_t chr_idx{};
  int64_t nbins = 0;

  for (const auto &path : path_to_cmatrices) {
    const auto header = ContactMatrix<uint32_t>::parse_header(path);
    fmt::print(stderr, FMT_STRING("Processing '{}' ({:.2f} Mbp)..."), header.chr_name,
               header.end / 1.0e6);  // NOLINT
    const auto t0 = absl::Now();
    const auto cmatrix = ContactMatrix<uint32_t>(path);
    const auto bs = header.bin_size;
    chr_idx = ENUM.getMemberIndex(header.chr_name);
    idx_chrom_offset_buff.push_back(nbins);
    nbins = write_bins(f, chr_idx, header.end, bs, bin_chrom_buff, bin_pos_buff, nbins);  // NOLINT

    for (auto i = 0UL; i < header.start / bs; ++i) {
      idx_bin1_offset_buff.push_back(n);
      if (idx_bin1_offset_buff.size() == BUFF_SIZE) {
        idx_bin1_offset_h5df_offset = write_vect_of_int(
            idx_bin1_offset_buff, f, "indexes/bin1_offset", idx_bin1_offset_h5df_offset);
        idx_bin1_offset_buff.clear();
      }
    }

    for (auto i = 0UL; i < cmatrix.n_cols(); ++i) {
      idx_bin1_offset_buff.push_back(n);
      for (auto j = i; j < i + cmatrix.n_rows() && j < cmatrix.n_cols(); ++j) {
        if (const auto m = cmatrix.get(i, j); m != 0) {
          pixel_b1_idx_buff.push_back(offset + static_cast<int64_t>(i));
          pixel_b2_idx_buff.push_back(offset + static_cast<int64_t>(j));
#ifndef NDEBUG
          if (offset + static_cast<int64_t>(i) > offset + static_cast<int64_t>(j)) {
            throw std::runtime_error(
                fmt::format(FMT_STRING("b1 > b2: b1={}; b2={}; offset={}; m={}\n"),
                            pixel_b1_idx_buff.back(), pixel_b2_idx_buff.back(), offset, m));
          }
#endif
          pixel_count_buff.push_back(m);
          ++n;
          if (pixel_b1_idx_buff.size() == BUFF_SIZE) {
            write_pixels_to_file();
          }
        }
      }
      if (idx_bin1_offset_buff.size() == BUFF_SIZE) {
        idx_bin1_offset_h5df_offset = write_vect_of_int(
            idx_bin1_offset_buff, f, "indexes/bin1_offset", idx_bin1_offset_h5df_offset);
        idx_bin1_offset_buff.clear();
      }
    }
    offset = nbins;
    fmt::print(stderr, FMT_STRING(" DONE in {}!\n"), absl::FormatDuration(absl::Now() - t0));
  }

  if (!pixel_b1_idx_buff.empty()) {
    write_pixels_to_file();
  }

  idx_chrom_offset_buff.push_back(nbins);
  idx_bin1_offset_buff.push_back(n);

  write_vect_of_int(idx_chrom_offset_buff, f, "indexes/chrom_offset");
  write_vect_of_int(idx_bin1_offset_buff, f, "indexes/bin1_offset", idx_bin1_offset_h5df_offset);

  H5::DataSpace att_space(H5S_SCALAR);
  auto att = f.createAttribute("nbins", H5::PredType::NATIVE_INT, att_space);
  att.write(H5::PredType::NATIVE_INT, &nbins);
  att = f.createAttribute("nnz", H5::PredType::NATIVE_INT, att_space);
  att.write(H5::PredType::NATIVE_INT, &n);
}

void write_metadata(H5::H5File &f, int32_t bin_size, std::string_view assembly_name) {
  H5::StrType STR_TYPE(H5::PredType::C_S1, H5T_VARIABLE);
  STR_TYPE.setCset(H5T_CSET_UTF8);
  H5::DataSpace att_space(H5S_SCALAR);
  int32_t int_buff{};
  std::string str_buff{};

  auto att = f.createAttribute("format", STR_TYPE, att_space);
  str_buff = "HDF5::Cooler";
  att.write(STR_TYPE, str_buff);

  att = f.createAttribute("format-version", H5::PredType::NATIVE_INT, att_space);
  int_buff = 3;
  att.write(H5::PredType::NATIVE_INT, &int_buff);

  att = f.createAttribute("bin-type", STR_TYPE, att_space);
  str_buff = "fixed";
  att.write(STR_TYPE, str_buff);

  att = f.createAttribute("bin-size", H5::PredType::NATIVE_INT, att_space);
  att.write(H5::PredType::NATIVE_INT, &bin_size);

  att = f.createAttribute("storage-mode", STR_TYPE, att_space);
  str_buff = "symmetric-upper";
  att.write(STR_TYPE, str_buff);

  if (!assembly_name.empty()) {
    str_buff = std::string{assembly_name};
    att = f.createAttribute("assembly-name", STR_TYPE, att_space);
    att.write(STR_TYPE, &str_buff);
  }

  att = f.createAttribute("generated-by", STR_TYPE, att_space);
  str_buff = fmt::format("ModLE-v{}", "0.0.1");  // TODO make ModLE ver a tunable
  att.write(STR_TYPE, str_buff);

  att = f.createAttribute("creation-date", STR_TYPE, att_space);
  str_buff = absl::FormatTime(absl::Now(), absl::UTCTimeZone());
  att.write(STR_TYPE, str_buff);
}

void modle_to_cooler(const std::vector<std::string_view> &path_to_cmatrices,
                     std::string_view path_to_output) {
  using CMatrix = ContactMatrix<int64_t>;
  std::vector<CMatrix::Header> headers;
  headers.reserve(path_to_cmatrices.size());

  std::transform(path_to_cmatrices.begin(), path_to_cmatrices.end(), std::back_inserter(headers),
                 [](auto path_to_cmatrix) { return CMatrix::parse_header(path_to_cmatrix); });

  auto cooler_file = init_file(path_to_output, true);
  write_metadata(cooler_file,
                 headers[0].bin_size);  // TODO: figure out a cleaner way to provide bin sizes
  auto CHR_ENUM = write_chroms(cooler_file, headers);

  write_contacts(cooler_file, path_to_cmatrices, CHR_ENUM);
}
}  // namespace modle::cooler