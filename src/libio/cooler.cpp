#include "modle/cooler.hpp"

#include <H5Cpp.h>
#include <absl/container/flat_hash_map.h>

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

#include "modle/contacts.hpp"
#include "modle/suppress_compiler_warnings.hpp"

namespace modle::cooler {

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
    DISABLE_WARNING_PUSH
    DISABLE_WARNING_C_VLA
    const auto *cbuff = reinterpret_cast<int32_t(*)[BUFF_SIZE]>(data.data());  // NOLINT
    DISABLE_WARNING_POP
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

H5::H5File init_file(std::string_view path_to_output_file, bool force_overwrite) {
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

hsize_t write_bins(H5::H5File &f, int32_t chrom, int64_t length, int64_t bin_size,
                   std::vector<int32_t> &buff32, std::vector<int64_t> &buff64, hsize_t file_offset,
                   hsize_t BUFF_SIZE) {
  buff32.resize(BUFF_SIZE);
  buff64.resize(BUFF_SIZE);

  std::fill(buff32.begin(), buff32.end(), chrom);

  int64_t start = 0;
  int64_t end = std::min(length, bin_size);
  hsize_t bins_processed = 0;
  const int64_t nbins = (length / bin_size) + 1;
  const int64_t nchunks = (nbins / static_cast<int64_t>(BUFF_SIZE)) + 1;

  for (auto i = 0; i < nchunks; ++i) {
    const auto chunk_size = std::min(BUFF_SIZE, static_cast<hsize_t>(nbins) - bins_processed);
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

absl::flat_hash_map<std::string, H5O_info_t> parse_file_structure(H5::H5File &f) {
  DISABLE_WARNING_PUSH
  DISABLE_WARNING_UNUSED_PARAMETER
  DISABLE_WARNING_SHADOW  // I believe the -Wshadow raised by GCC is a false positive
      auto op_func_L =
          [](hid_t loc_id, const char *name, const H5L_info_t *info, void *operator_data) {
            auto op_func = [](hid_t loc_id, const char *name, const H5O_info_t *info,
                              void *operator_data) {
              auto *structure =
                  reinterpret_cast<absl::flat_hash_map<std::string, H5O_info_t> *>(  // NOLINT
                      operator_data);
              assert(!structure->contains(name));  // NOLINT
              structure->emplace(std::string{name}, *info);
              return 0;
            };

            H5O_info_t infobuf;
            if (const auto status =
                    H5Oget_info_by_name(loc_id, name, &infobuf, H5O_INFO_BASIC, H5P_DEFAULT);
                status != 0) {
              throw fmt::system_error(status, "An error occurred during traversal");
            }

            return op_func(loc_id, name, &infobuf, operator_data);
          };
  DISABLE_WARNING_POP

  absl::flat_hash_map<std::string, H5O_info_t> obj_info;
  if (const auto status = H5Lvisit(f.getId(), H5_INDEX_NAME, H5_ITER_NATIVE,
                                   static_cast<H5L_iterate2_t>(op_func_L), &obj_info);
      status != 0) {
    throw fmt::system_error(status, "An error occurred during traversal");
  }
  std::vector<std::string> attributes;
  for (const auto &[name, info] : obj_info) {
    if (info.type == H5O_TYPE_GROUP) {
      auto g = f.openGroup(name);
      for (auto i = 0; i < g.getNumAttrs(); ++i) {
        auto attr = g.openAttribute(static_cast<uint32_t>(i));
        attributes.emplace_back(attr.getName());
      }
    }
  }
  for (const auto &attr : attributes) {
    H5O_info_t info;
    info.type = H5O_TYPE_UNKNOWN;
    obj_info.emplace(attr, info);
  }
  std::vector<std::string> attrs;
  return obj_info;
}

void validate_cooler(H5::H5File &f) {
  H5::Exception::dontPrint();
  H5::StrType STR_TYPE(H5::PredType::C_S1, H5T_VARIABLE);

  const auto file_structure = parse_file_structure(f);
  const auto is_multires = file_structure.contains("resolutions");
  const auto is_singlecell = file_structure.contains("cells");

  if (!is_multires && !is_singlecell) {
    for (const auto &group : {"chroms", "bins", "pixels", "indexes"}) {
      if (!file_structure.contains(group)) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("Mandatory group '{}' was not found"), group));
      }
    }
  }

  if (is_multires) {
    throw std::runtime_error("Multires cooler files are not yet supported");
  }
  if (is_singlecell) {
    throw std::runtime_error("Single-cell cooler files are not yet supported");
  }
}

H5::H5File open_for_reading(std::string_view path_to_file, bool validate) {
  try {
    H5::Exception::dontPrint();
    auto f = H5::H5File(std::string{path_to_file}, H5F_ACC_RDONLY);
    if (validate) {
      try {
        validate_cooler(f);
      } catch (const std::runtime_error &e) {
        throw std::runtime_error(fmt::format(FMT_STRING("File validation failed: {}"), e.what()));
      }
    }
    return f;

  } catch (const H5::Exception &e) {
    if (e.getDetailMsg() == "H5Aopen failed") {
      throw std::runtime_error(
          fmt::format(FMT_STRING("An error occurred while opening file '{}' for reading: Function "
                                 "'{}' failed with error '{}'"),
                      path_to_file, e.getFuncName(), e.getDetailMsg()));
    }
    throw std::runtime_error(
        fmt::format(FMT_STRING("An error occurred while opening file '{}' for reading: function "
                               "'{}' failed with error '{}'"),
                    path_to_file, e.getFuncName(), e.getDetailMsg()));
  } catch (const std::runtime_error &e) {
    if (absl::StartsWith(e.what(), "An error occurred during traversal")) {
      throw std::runtime_error(fmt::format(FMT_STRING("{} of file '{}'"), e.what(), path_to_file));
    }
    throw std::runtime_error(fmt::format(
        FMT_STRING("An error occurred while reading file '{}': {}"), path_to_file, e.what()));
  }
}

hsize_t read_str(H5::H5File &f, std::string_view dataset_name, std::string &BUFF,
                 hsize_t file_offset) {
  try {
    H5::Exception::dontPrint();

    const H5std_string d{dataset_name};
    constexpr const hsize_t DIMS = 1;
    constexpr const hsize_t RANK = 1;

    auto dataset = f.openDataSet(d);
    auto STR_TYPE = dataset.getDataType();
    auto file_space = dataset.getSpace();

    auto mem_space{H5::DataSpace(RANK, &DIMS)};

    file_space = dataset.getSpace();
    file_space.selectHyperslab(H5S_SELECT_SET, &RANK, &file_offset);
    dataset.read(BUFF, STR_TYPE, mem_space, file_space);

    return file_offset + 1;

  } catch (const H5::Exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Failed to read data into a string (dataset = {}; file offset = "
                               "{}): Function '{}' failed with error '{}'"),
                    dataset_name, file_offset, e.getFuncName(), e.getDetailMsg()));
  }
}

std::vector<std::string> read_vect_of_str(H5::H5File &f, std::string_view dataset_name,
                                          hsize_t file_offset) {
  std::vector<std::string> buff;
  read_vect_of_str(f, dataset_name, buff, file_offset);
  return buff;
}

std::size_t read_vect_of_str(H5::H5File &f, std::string_view dataset_name,
                             std::vector<std::string> &BUFF, hsize_t file_offset) {
  try {
    H5::Exception::dontPrint();

    const H5std_string d{dataset_name};
    constexpr const hsize_t DIMS = 1;
    constexpr const hsize_t RANK = 1;

    auto dataset = f.openDataSet(d);
    auto STR_TYPE = dataset.getDataType();
    auto file_space = dataset.getSpace();

    const auto N_STRINGS = static_cast<hsize_t>(file_space.getSimpleExtentNpoints());
    BUFF.resize(N_STRINGS);

    auto mem_space{H5::DataSpace(RANK, &DIMS)};

    for (hsize_t i = 0U; i < N_STRINGS; ++i) {
      file_space = dataset.getSpace();
      file_space.selectHyperslab(H5S_SELECT_SET, &RANK, &i);
      dataset.read(BUFF[i], STR_TYPE, mem_space, file_space);
    }

    return N_STRINGS;

  } catch (const H5::Exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("Failed to read data into a vector of strings (dataset = {}; file offset = "
                   "{}): Function '{}' failed with error '{}'"),
        dataset_name, file_offset, e.getFuncName(), e.getDetailMsg()));
  }
}

std::size_t read_chr_offset_idx(H5::H5File &f, std::string_view chr_name) {
  try {
    H5::Exception::dontPrint();

    constexpr const hsize_t DIMS = 1;
    constexpr const hsize_t RANK = 1;
    std::string BUFF;

    auto dataset = f.openDataSet("chroms/name");
    auto STR_TYPE = dataset.getDataType();
    auto file_space = dataset.getSpace();

    const auto nchroms = static_cast<hsize_t>(file_space.getSimpleExtentNpoints());

    auto mem_space{H5::DataSpace(RANK, &DIMS)};

    for (hsize_t i = 0U; i < nchroms; ++i) {
      file_space = dataset.getSpace();
      file_space.selectHyperslab(H5S_SELECT_SET, &RANK, &i);
      dataset.read(BUFF, STR_TYPE, mem_space, file_space);
      if (BUFF == chr_name) {
        return i;
      }
    }

    throw std::runtime_error(
        fmt::format(FMT_STRING("Unable to find a chromosome named '{}'"), chr_name));
  } catch (const H5::Exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("An error occurred while looking for chromosome "
                               "'{}': Function '{}' failed with error '{}'"),
                    chr_name, e.getFuncName(), e.getDetailMsg()));
  }
}

std::vector<int64_t> read_chr_offset_idx(H5::H5File &f) {
  std::vector<int64_t> idx;
  read_vect_of_int(f, "indexes/chrom_offset", idx, 0);
  return idx;
}

std::vector<int64_t> read_bin1_offset_idx(H5::H5File &f) {
  std::vector<int64_t> idx;
  read_vect_of_int(f, "indexes/bin1_offset", idx, 0);
  return idx;
}

std::vector<int64_t> read_bin1_offset_idx(H5::H5File &f, std::string_view chr_name) {
  const auto chr_idx = read_chr_offset_idx(f, chr_name);
  return read_bin1_offset_idx(f, chr_idx);
}

std::vector<int64_t> read_bin1_offset_idx(H5::H5File &f, std::size_t chr_idx) {
  const auto chr_offset_idx = read_chr_offset_idx(f);
  assert(chr_idx < chr_offset_idx.size());  // NOLINT
  const auto chr_start_bin = static_cast<std::size_t>(chr_offset_idx[chr_idx]);
  const auto chr_end_bin = static_cast<std::size_t>(chr_offset_idx[chr_idx + 1]);
  assert(chr_end_bin >= chr_start_bin);  // NOLINT

  std::vector<int64_t> bin1_offset_idx(chr_end_bin - chr_start_bin);
  read_vect_of_int(f, "indexes/bin1_offset", bin1_offset_idx, chr_start_bin);

  return bin1_offset_idx;
}

std::pair<int64_t, int64_t> read_chrom_pixels_boundaries(H5::H5File &f, std::string_view chr_name) {
  const auto chr_idx = read_chr_offset_idx(f, chr_name);
  return read_chrom_pixels_boundaries(f, chr_idx);
}

std::pair<int64_t, int64_t> read_chrom_pixels_boundaries(H5::H5File &f, std::size_t chr_idx) {
  const auto chr_offset_idx = read_chr_offset_idx(f);
  assert(chr_idx < chr_offset_idx.size());  // NOLINT
  const auto chr_start_bin = static_cast<std::size_t>(chr_offset_idx[chr_idx]);
  const auto chr_end_bin = static_cast<std::size_t>(chr_offset_idx[chr_idx + 1]);
  assert(chr_end_bin >= chr_start_bin);  // NOLINT

  std::pair<int64_t, int64_t> pixel_boundaries{};
  read_int(f, "indexes/bin1_offset", pixel_boundaries.first, chr_start_bin);
  read_int(f, "indexes/bin1_offset", pixel_boundaries.second, chr_end_bin);

  return pixel_boundaries;
}

ContactMatrix<uint32_t> cooler_to_cmatrix(H5::H5File &f, int64_t bin_offset,
                                          const std::vector<int64_t> &bin1_offset_idx,
                                          std::size_t diagonal_width, std::size_t bin_size) {
  // TODO: Add check for Cooler bin size
  const auto nrows = diagonal_width / bin_size;
  return cooler_to_cmatrix(f, bin_offset, bin1_offset_idx, nrows);
}

ContactMatrix<uint32_t> cooler_to_cmatrix(H5::H5File &f, int64_t bin_offset,
                                          const std::vector<int64_t> &bin1_offset_idx,
                                          std::size_t nrows) {
  ContactMatrix<uint32_t> cmatrix(nrows, bin1_offset_idx.size());
  std::vector<int64_t> bin1_BUFF;
  std::vector<int64_t> bin2_BUFF;
  std::vector<int64_t> count_BUFF;

  for (auto i = 1UL; i < bin1_offset_idx.size(); ++i) {
    const auto file_offset = static_cast<hsize_t>(bin1_offset_idx[i - 1]);
    const auto buff_size =
        std::min(nrows, static_cast<std::size_t>(bin1_offset_idx[i] - bin1_offset_idx[i - 1])) + 1;
    if (buff_size == 0) {
      continue;
    }

    bin1_BUFF.resize(buff_size);
    bin2_BUFF.resize(buff_size);
    count_BUFF.resize(buff_size);

    read_vect_of_int(f, "pixels/bin1_id", bin1_BUFF, file_offset);
    read_vect_of_int(f, "pixels/bin2_id", bin2_BUFF, file_offset);
    read_vect_of_int(f, "pixels/count", count_BUFF, file_offset);

    assert(bin1_BUFF.size() == buff_size);   // NOLINT
    assert(bin2_BUFF.size() == buff_size);   // NOLINT
    assert(count_BUFF.size() == buff_size);  // NOLINT

    for (auto j = 0UL; j < buff_size; ++j) {
      assert(count_BUFF[j] != 0);  // NOLINT
      DISABLE_WARNING_PUSH
      DISABLE_WARNING_SIGN_CONVERSION
      DISABLE_WARNING_SIGN_COMPARE
      if (bin2_BUFF[j] - bin_offset >= i + nrows - 1) {
        break;
      }
      // fmt::print(stderr, "m[{}][{}]={} (nrows={}; ncols={})\n", bin2_BUFF[j] - bin_offset,
      //           bin1_BUFF[j] - bin_offset, count_BUFF[j], cmatrix.n_rows(), cmatrix.n_cols());
      cmatrix.set(bin2_BUFF[j] - bin_offset, bin1_BUFF[j] - bin_offset, count_BUFF[j]);

      DISABLE_WARNING_POP
    }
  }
  return cmatrix;
}

ContactMatrix<uint32_t> cooler_to_cmatrix(std::string_view path_to_file, std::string_view chr_name,
                                          std::size_t diagonal_width, std::size_t bin_size) {
  // TODO: Add check for Cooler bin size
  const auto nrows = diagonal_width / bin_size;
  return cooler_to_cmatrix(path_to_file, chr_name, nrows);
}

ContactMatrix<uint32_t> cooler_to_cmatrix(std::string_view path_to_file, std::string_view chr_name,
                                          std::size_t nrows) {
  auto f = open_for_reading(path_to_file);
  const auto bin1_offset_idx = cooler::read_bin1_offset_idx(f, chr_name);
  const auto chr_idx = cooler::read_chr_offset_idx(f, chr_name);
  const auto bin_offset = cooler::read_chr_offset_idx(f)[chr_idx];

  return cooler_to_cmatrix(f, bin_offset, bin1_offset_idx, nrows);
}

}  // namespace modle::cooler