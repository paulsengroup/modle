#include "modle/cooler.hpp"

#include <H5Cpp.h>

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

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

}  // namespace modle::cooler