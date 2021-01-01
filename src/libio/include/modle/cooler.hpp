#pragma once

#include <H5Cpp.h>

#include <filesystem>
#include <memory>
#include <string_view>
#include <vector>

#include "modle/contacts.hpp"

namespace modle::cooler {
class Cooler {
 private:
  inline static H5::StrType generate_default_str_type();

 public:
  enum IO_MODE { READ_ONLY, WRITE_ONLY };
  inline static H5::StrType STR_TYPE{generate_default_str_type()};    // NOLINT
  inline static H5::PredType INT64_TYPE{H5::PredType::NATIVE_INT64};  // NOLINT
  inline static H5::PredType INT32_TYPE{H5::PredType::NATIVE_INT32};  // NOLINT
  enum Flavor {
    UNK = 0,
    AUTO = 1,
    COOL = 2,
    MCOOL = 3,
    SCOOL = 4  // For the time being, we do not support SCOOL
  };
  enum Groups { CHR = 0, BIN = 1, PXL = 2, IDX = 3 };
  enum Datasets {
    CHR_LEN = 0,
    CHR_NAME = 1,
    BIN_CHROM = 2,
    BIN_START = 3,
    BIN_END = 4,
    PXL_B1 = 5,
    PXL_B2 = 6,
    PXL_COUNT = 7,
    IDX_BIN1 = 8,
    IDX_CHR = 9
  };

  Cooler() = delete;
  inline explicit Cooler(std::string_view path_to_file, IO_MODE mode = READ_ONLY,
                         std::size_t bin_size = 0, Flavor flavor = AUTO, bool validate = true,
                         uint8_t compression_lvl = 9,                    // NOLINT
                         std::size_t chunk_size = 1024 * 1024ULL,        // 1 MB NOLINT
                         std::size_t cache_size = 16 * 1024 * 1024ULL);  // 16 MB NOLINT

  [[nodiscard]] inline static Flavor detect_file_flavor(H5::H5File& f);

  [[nodiscard]] inline static bool validate_file_format(H5::H5File& f, Flavor expected_flavor,
                                                        IO_MODE mode, std::size_t bin_size = 0,
                                                        bool throw_on_failure = true);

  [[nodiscard]] inline static bool validate_cool_flavor(H5::H5File& f, std::size_t bin_size,
                                                        std::string_view root_path = "/",
                                                        bool throw_on_failure = true);

  [[nodiscard]] inline static bool validate_multires_cool_flavor(H5::H5File& f,
                                                                 std::size_t bin_size,
                                                                 std::string_view root_path = "/",
                                                                 bool throw_on_failure = true);
  template <typename I>
  inline void write_cmatrix_to_file(const ContactMatrix<I>& cmatrix,
                                    bool wipe_before_writing = false);

  inline void write_metadata();
  [[nodiscard]] inline bool is_read_only() const;

 private:
  std::filesystem::path _path_to_file;
  IO_MODE _mode;
  std::size_t _bin_size;
  Flavor _flavor;
  std::unique_ptr<H5::H5File> _fp{nullptr};
  std::vector<H5::Group> _groups{};
  std::vector<H5::DataSet> _datasets{};
  std::unique_ptr<H5::DataSpace> _mem_space{nullptr};
  std::vector<std::pair<H5::DataSpace, hsize_t>> _fspaces{};  // fspaces + file offset

  uint8_t _compression_lvl;
  hsize_t _chunk_size;
  hsize_t _cache_size;

  std::unique_ptr<H5::DSetCreatPropList> _cprop_str{nullptr};
  std::unique_ptr<H5::DSetCreatPropList> _cprop_int32{nullptr};
  std::unique_ptr<H5::DSetCreatPropList> _cprop_int64{nullptr};

  std::unique_ptr<H5::DSetAccPropList> _aprop_str{nullptr};
  std::unique_ptr<H5::DSetAccPropList> _aprop_int32{nullptr};
  std::unique_ptr<H5::DSetAccPropList> _aprop_int64{nullptr};

  [[nodiscard]] inline static std::unique_ptr<H5::H5File> open_file(
      const std::filesystem::path& path, IO_MODE mode, std::size_t bin_size = 0,
      Flavor flavor = AUTO, bool validate = true);
  [[nodiscard]] inline static std::vector<H5::Group> open_groups(H5::H5File& f,
                                                                 bool create_if_not_exist = false);
  [[nodiscard]] inline static std::string flavor_to_string(Flavor f);
  [[nodiscard]] inline static bool check_version(int64_t min_ver = 2, int64_t max_ver = 3,
                                                 bool throw_on_failure = false);

  template <typename T1, typename T2>
  std::unique_ptr<H5::DSetCreatPropList> generate_default_cprop(hsize_t chunk_size,
                                                                uint8_t compression_lvl, T1 type,
                                                                T2 fill_value);
  template <typename T>
  inline static std::unique_ptr<H5::DSetAccPropList> generate_default_aprop(T type,
                                                                            hsize_t chunk_size,
                                                                            hsize_t cache_size);
};

}  // namespace modle::cooler

#include "../../cooler_impl.hpp"