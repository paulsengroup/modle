#include <H5Cpp.h>

#include <catch2/catch.hpp>
#include <filesystem>
#include <limits>
#include <string>
#include <vector>

#include "modle/hdf5.hpp"

namespace modle::test::hdf5 {
using namespace modle::hdf5;

inline const std::filesystem::path test_dir{"/tmp/modle/unit_tests"};  // NOLINT
[[nodiscard]] inline H5::StrType init_str_type() {
  auto st = H5::StrType(H5::PredType::C_S1, H5T_VARIABLE);
  st.setStrpad(H5T_STR_NULLPAD);
  st.setCset(H5T_CSET_ASCII);

  return st;
}

inline H5::DataSet init_test_str_dataset(H5::H5File& f, std::string_view path = "/test") {
  constexpr hsize_t rank{1};
  constexpr hsize_t chunk_dims{1024};
  const H5std_string fill_var{"\0"};
  constexpr hsize_t max_dims{H5S_UNLIMITED};
  const auto str_type{init_str_type()};

  H5::DSetCreatPropList cprop{};
  cprop.setChunk(rank, &chunk_dims);
  cprop.setFillValue(str_type, &fill_var);
  cprop.setDeflate(1);

  auto mem_space{H5::DataSpace(rank, &rank, &max_dims)};

  return f.createDataSet(std::string{path.data(), path.size()}, str_type, mem_space, cprop);
}

inline H5::DataSet init_test_int64_dataset(H5::H5File& f, std::string_view path = "/test") {
  constexpr hsize_t rank{1};
  constexpr hsize_t chunk_dims{1024};
  constexpr int64_t fill_var{0};
  constexpr hsize_t max_dims{H5S_UNLIMITED};

  H5::DSetCreatPropList cprop{};
  cprop.setChunk(rank, &chunk_dims);
  cprop.setFillValue(H5::PredType::NATIVE_INT64, &fill_var);
  cprop.setDeflate(1);

  auto mem_space{H5::DataSpace(rank, &rank, &max_dims)};

  return f.createDataSet(std::string{path.data(), path.size()}, H5::PredType::NATIVE_INT64,
                         mem_space, cprop);
}

TEST_CASE("read_write_strings HDF5", "[io][hdf5][short]") {
  const auto test_file = test_dir / "rw_strings.hdf5";
  std::filesystem::create_directories(test_dir);

  std::vector<std::string> v{"0", "QWERTY", "ABCDE", "0123456789", "!\"#¤%&/()=?"};

  H5::H5File f(test_file.string(), H5F_ACC_TRUNC);
  auto dataset = init_test_str_dataset(f);

  CHECK(write_strings(v, dataset, init_str_type(), 0) == v.size());
  for (auto i = 0UL; i < v.size(); ++i) {
    const auto buff = read_str(dataset, i);
    CHECK(buff == v[i]);
  }

  const auto buffv = read_strings(dataset, 0);
  REQUIRE(buffv.size() == v.size());
  for (auto i = 0UL; i < v.size(); ++i) {
    CHECK(v[i] == buffv[i]);
  }
  std::filesystem::remove_all(test_dir);
  if (const auto& p = test_dir.parent_path(); std::filesystem::is_empty(p)) {
    std::filesystem::remove(p);
  }
}

TEST_CASE("read_write_ints HDF5", "[io][hdf5][short]") {
  const auto test_file = test_dir / "rw_ints.hdf5";
  std::filesystem::create_directories(test_dir);
  std::vector<int64_t> v{std::numeric_limits<int64_t>::min(), -10, 0, 10,  // NOLINT
                         std::numeric_limits<int64_t>::max()};
  H5::H5File f(test_file.string(), H5F_ACC_TRUNC);
  auto dataset = init_test_int64_dataset(f);
  int64_t buff{};
  for (auto i = 0UL; i < v.size(); ++i) {
    CHECK(write_number(v[i], dataset, i) == i + 1);
    CHECK(read_number(dataset, buff, i) == i + 1);
    CHECK(buff == v[i]);
  }

  std::vector<int64_t> buffv(v.size());
  CHECK(write_numbers(v, dataset, 0) == v.size());
  CHECK(read_numbers(dataset, buffv, 0));
  REQUIRE(buffv.size() == v.size());
  for (auto i = 0UL; i < v.size(); ++i) {
    CHECK(v[i] == buffv[i]);
  }
  std::filesystem::remove_all(test_dir);
  if (const auto& p = test_dir.parent_path(); std::filesystem::is_empty(p)) {
    std::filesystem::remove(p);
  }
}

TEST_CASE("group_exists HDF5", "[io][hdf5][short]") {
  const auto test_file = test_dir / "group_exists.hdf5";
  std::filesystem::create_directories(test_dir);
  H5::H5File f(test_file.string(), H5F_ACC_TRUNC);
  (void)init_test_int64_dataset(f);
  f.createGroup("/g1");
  f.createGroup("/g1/1");

  CHECK(group_exists(f, "g1"));
  CHECK(group_exists(f, "1", "g1"));
  CHECK(group_exists(f, "1", "/g1"));
  CHECK(group_exists(f, "1", "/g1/"));
  CHECK(!group_exists(f, "g2"));
  CHECK(!group_exists(f, "2", "g2"));
  CHECK_THROWS_WITH(group_exists(f, "test"),
                    Catch::Matchers::EndsWith("exists but is not a group"));

  std::filesystem::remove_all(test_dir);
  if (const auto& p = test_dir.parent_path(); std::filesystem::is_empty(p)) {
    std::filesystem::remove(p);
  }
}

TEST_CASE("check_dataset_type HDF5", "[io][hdf5][short]") {
  const auto test_file = test_dir / "dset_exisrs.hdf5";
  std::filesystem::create_directories(test_dir);
  H5::H5File f(test_file.string(), H5F_ACC_TRUNC);
  f.createGroup("/g1");
  auto int_dataset = init_test_int64_dataset(f);
  auto str_dataset = init_test_str_dataset(f, "/g1/test");
  CHECK(dataset_exists(f, "/test"));
  CHECK(!dataset_exists(f, "/test2"));

  CHECK(check_dataset_type(int_dataset, H5::PredType::NATIVE_INT64));
  CHECK(check_dataset_type(str_dataset, init_str_type()));
  CHECK(!check_dataset_type(str_dataset, H5::PredType::NATIVE_INT64, false));

  CHECK_THROWS_WITH(check_dataset_type(int_dataset, H5::PredType::NATIVE_INT),
                    Catch::Matchers::Contains("incorrect datasize"));
  CHECK_THROWS_WITH(check_dataset_type(int_dataset, H5::PredType::NATIVE_FLOAT),
                    Catch::Matchers::EndsWith("incorrect datatype"));
  CHECK_THROWS_WITH(check_dataset_type(int_dataset, H5::PredType::NATIVE_UINT64),
                    Catch::Matchers::Contains("incorrect signedness"));

  std::filesystem::remove_all(test_dir);
  if (const auto& p = test_dir.parent_path(); std::filesystem::is_empty(p)) {
    std::filesystem::remove(p);
  }
}

/*
TEST_CASE("read_chr_offset_idx 001", "[parsers][cooler][short]") {
  auto f = cooler::open_for_reading(
      "test/data/cooler/Dixon2012-H1hESC-HindIII-allreps-filtered.100kb.cool");
  const auto chr1_idx = cooler::read_chr_offset_idx(f, "chr1");
  const auto chr5_idx = cooler::read_chr_offset_idx(f, "chr5");
  const auto chrX_idx = cooler::read_chr_offset_idx(f, "chrX");

  CHECK(chr1_idx == 0);
  CHECK(chr5_idx == 4);
  CHECK(chrX_idx == 22);

  CHECK_THROWS_WITH(cooler::read_chr_offset_idx(f, "chr0"),
                    "Unable to find a chromosome named 'chr0'");
}

TEST_CASE("read_bin1_offset 001", "[parsers][cooler][short]") {
  auto f = cooler::open_for_reading(
      "test/data/cooler/Dixon2012-H1hESC-HindIII-allreps-filtered.100kb.cool");
  const auto pixels_idx = cooler::read_bin1_offset_idx(f, "chr5");

  CHECK(pixels_idx.front() == 22082852);
  CHECK(pixels_idx.back() == 25968678);
}

TEST_CASE("read_pixel_boundaries 001", "[parsers][cooler][short]") {
  auto f = cooler::open_for_reading(
      "test/data/cooler/Dixon2012-H1hESC-HindIII-allreps-filtered.100kb.cool");
  const auto pixel_boundaries = cooler::read_chrom_pixels_boundaries(f, "chr5");

  CHECK(pixel_boundaries.first == 22082852);
  CHECK(pixel_boundaries.second == 25968678);
}

TEST_CASE("cooler_to_cmatrix 001", "[parsers][cooler][short]") {
  const auto test_file = test_dir / "cooler_to_cmatrix_002.cool";
  std::filesystem::remove(test_file);
  const std::string chr_name = "chr_test";

  std::filesystem::create_directories(test_file.parent_path());

  const std::size_t nrows = 5;
  const std::size_t ncols = 20;
  const std::size_t bin_size = 1'000;
  const auto diagonal_width = nrows * bin_size;

  std::vector<ContactMatrix<uint32_t>> cmatrices;
  const auto ref_cmatrix = cmatrices.emplace_back(nrows, ncols, true);

  cooler::write_modle_cmatrices_to_cooler(cmatrices, {chr_name}, {0}, {(ncols * bin_size) - 1},
                                          {(ncols * bin_size) - 1}, bin_size, test_file.string(),
                                          true);

  const auto cmatrix =
      cooler::cooler_to_cmatrix(test_file.string(), chr_name, diagonal_width, bin_size);

  CHECK(ref_cmatrix.n_rows() == cmatrix.n_rows());
  CHECK(ref_cmatrix.n_cols() == cmatrix.n_cols());
  CHECK(ref_cmatrix.get_tot_contacts() == cmatrix.get_tot_contacts());
  CHECK(ref_cmatrix.get_n_of_missed_updates() == cmatrix.get_n_of_missed_updates());

  for (auto i = 0UL; i < ref_cmatrix.n_cols(); ++i) {
    for (auto j = 0UL; j < ref_cmatrix.n_cols(); ++j) {
      CHECK(ref_cmatrix.get(i, j) == cmatrix.get(i, j));
    }
  }
  std::filesystem::remove_all(test_file.parent_path());
  if (const auto& p = test_file.parent_path().parent_path(); std::filesystem::is_empty(p)) {
    std::filesystem::remove(p);
  }
}
 */

}  // namespace modle::test::hdf5
