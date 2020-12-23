#include "modle/cooler.hpp"

#include <catch2/catch.hpp>
#include <filesystem>
#include <limits>
#include <string>

namespace modle::cooler::test {

inline const std::filesystem::path test_dir{"/tmp/modle/unit_tests"};  // NOLINT

TEST_CASE("read_write_strings HDF5", "[parsers][hdf5][short]") {
  const auto test_file = test_dir / "rw_strings.hdf5";
  std::filesystem::create_directories(test_dir);
  std::vector<std::string> v{"0", "QWERTY", "ABCDE", "0123456789", "!\"#¤%&/()=?"};
  H5::H5File f(test_file.string(), H5F_ACC_TRUNC);
  std::string buff;
  for (auto i = 0UL; i < v.size(); ++i) {
    cooler::write_str(v[i], f, "/test", v.back().size() + 1, i);
    cooler::read_str(f, "/test", buff, i);
    CHECK(buff == v[i]);
  }

  const auto buffv = cooler::read_vect_of_str(f, "/test", 0);
  REQUIRE(buffv.size() == v.size());
  for (auto i = 0UL; i < v.size(); ++i) {
    CHECK(v[i] == buffv[i]);
  }
  std::filesystem::remove_all(test_dir);
  if (const auto& p = test_dir.parent_path(); std::filesystem::is_empty(p)) {
    std::filesystem::remove(p);
  }
}

TEST_CASE("read_write_vect_of_strings HDF5", "[parsers][hdf5][short]") {
  const auto test_file = test_dir / "rw_strings.hdf5";
  std::filesystem::create_directories(test_dir);
  std::vector<std::string> v{"0", "QWERTY", "ABCDE", "0123456789", "!\"#¤%&/()=?"};
  H5::H5File f(test_file.string(), H5F_ACC_TRUNC);
  cooler::write_vect_of_str(v, f, "/test");

  auto buff = cooler::read_vect_of_str(f, "/test", 0);
  REQUIRE(v.size() == buff.size());
  for (auto i = 0UL; i < v.size(); ++i) {
    CHECK(v[i] == buff[i]);
  }

  buff.clear();
  cooler::read_vect_of_str(f, "/test", buff, 0);

  REQUIRE(v.size() == buff.size());
  for (auto i = 0UL; i < v.size(); ++i) {
    CHECK(v[i] == buff[i]);
  }
  std::filesystem::remove_all(test_dir);
  if (const auto& p = test_dir.parent_path(); std::filesystem::is_empty(p)) {
    std::filesystem::remove(p);
  }
}

TEST_CASE("read_write_ints HDF5", "[parsers][hdf5][short]") {
  const auto test_file = test_dir / "rw_ints.hdf5";
  std::filesystem::create_directories(test_dir);
  std::vector<int64_t> v{std::numeric_limits<int64_t>::min(), -10, 0, 10,  // NOLINT
                         std::numeric_limits<int64_t>::max()};
  H5::H5File f(test_file.string(), H5F_ACC_TRUNC);
  int64_t buff{};
  for (auto i = 0UL; i < v.size(); ++i) {
    cooler::write_int(v[i], f, "/test", i);
    cooler::read_int(f, "/test", buff, i);
    CHECK(buff == v[i]);
  }

  std::vector<int64_t> buffv(v.size());
  cooler::read_vect_of_int(f, "/test", buffv, 0);
  REQUIRE(buffv.size() == v.size());
  for (auto i = 0UL; i < v.size(); ++i) {
    CHECK(v[i] == buffv[i]);
  }
  std::filesystem::remove_all(test_dir);
  if (const auto& p = test_dir.parent_path(); std::filesystem::is_empty(p)) {
    std::filesystem::remove(p);
  }
}

TEST_CASE("read_write_vect_of_ints HDF5", "[parsers][hdf5][short]") {
  const auto test_file = test_dir / "rw_ints.hdf5";
  std::filesystem::create_directories(test_dir);
  std::vector<int64_t> v{std::numeric_limits<int64_t>::min(), -10, 0, 10,  // NOLINT
                         std::numeric_limits<int64_t>::max()};
  H5::H5File f(test_file.string(), H5F_ACC_TRUNC);
  cooler::write_vect_of_int(v, f, "/test");

  std::vector<int64_t> buffv(v.size());
  cooler::read_vect_of_int(f, "/test", buffv, 0);
  REQUIRE(v.size() == buffv.size());
  for (auto i = 0UL; i < v.size(); ++i) {
    CHECK(v[i] == buffv[i]);
  }

  std::filesystem::remove_all(test_dir);
  if (const auto& p = test_dir.parent_path(); std::filesystem::is_empty(p)) {
    std::filesystem::remove(p);
  }
}

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

}  // namespace modle::cooler::test
