#include "modle/cooler.hpp"

#include <catch2/catch.hpp>
#include <filesystem>
#include <string>
#include <string_view>

#include "modle/contacts.hpp"

namespace modle::cooler::test {

inline const std::filesystem::path test_dir{"/tmp/modle/unit_tests"};  // NOLINT
inline const std::filesystem::path data_dir{"test/data/cooler"};       // NOLINT

TEST_CASE("cooler ctor", "[io][cooler][short]") {
  const auto test_file = data_dir / "Dixon2012-H1hESC-HindIII-allreps-filtered.100kb.cool";

  {
    auto c = Cooler(test_file.string(), Cooler::READ_ONLY, 100'000);  // NOLINT
    CHECK(c.is_read_only());
  }

  H5::H5File f(test_file.string(), H5F_ACC_RDONLY);
  CHECK(Cooler::detect_file_flavor(f) == Cooler::COOL);
  CHECK(Cooler::validate_file_format(f, Cooler::COOL));
  CHECK_THROWS_WITH(!Cooler::validate_file_format(f, Cooler::MCOOL, Cooler::READ_ONLY, 1000),
                    Catch::Matchers::Contains("Expected format flavor MCOOL, found COOL"));
}

TEST_CASE("CMatrix to cooler", "[io][cooler][short]") {
  const auto test_file = test_dir / "cmatrix_to_cooler.cool";
  std::filesystem::create_directories(test_dir);

  constexpr uint64_t start = 0;
  constexpr uint64_t end = 5'000'000;
  constexpr uint64_t bin_size = 1000;
  constexpr uint64_t nrows = 25;
  constexpr uint64_t ncols = end / bin_size;

  auto c = Cooler(test_file.string(), Cooler::WRITE_ONLY, bin_size);

  ContactMatrix<int32_t> cmatrix(nrows, ncols, true);
  c.write_cmatrix_to_file(cmatrix, "chr0", start, end, end);

  /*
  std::filesystem::remove_all(test_dir);
  if (const auto& p = test_dir.parent_path(); std::filesystem::is_empty(p)) {
    std::filesystem::remove(p);
  }
   */
}

}  // namespace modle::cooler::test