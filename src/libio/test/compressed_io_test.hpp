#pragma once

#include <catch2/catch.hpp>  // for AssertionHandler, operator""_catch_sr, SourceLineInfo
#include <filesystem>        // for exists
#include <fstream>           // for ifstrea, getline
#include <iostream>
#include <string>  // for string, basic_string, operator==, char_traits, stoull

#include "modle/compressed_io.hpp"  // for Reader

namespace modle::test::libarchivexx {

[[maybe_unused]] void print_buffers(const std::string& buff1, const std::string& buff2) {
  std::cerr << "buff1='" << buff1 << "'\n"
            << "buff2='" << buff2 << "'\n";
}

TEST_CASE("Reader plain", "[reader][short]") {
  const std::filesystem::path ptext_file = "test/data/unit_tests/sample.bed9";

  modle::libarchivexx::Reader r(ptext_file);
  std::ifstream fp(ptext_file);

  std::string buff1;
  std::string buff2;
  while (r.getline(buff1) && std::getline(fp, buff2)) {
    CHECK(!!fp);
    CHECK(buff1 == buff2);
  }
  CHECK(!r.getline(buff2));
  CHECK(!std::getline(fp, buff2));
}

TEST_CASE("Reader plain sv", "[reader][short]") {
  const std::filesystem::path ptext_file = "test/data/unit_tests/sample.bed9";

  modle::libarchivexx::Reader r(ptext_file);
  std::ifstream fp(ptext_file);

  std::string_view buff1;
  std::string buff2;
  while (!(buff1 = r.getline()).empty() && std::getline(fp, buff2)) {
    CHECK(!!fp);
    CHECK(buff1 == buff2);
  }
  CHECK(!r.getline(buff2));
  CHECK(!std::getline(fp, buff2));
}

TEST_CASE("Reader gzip", "[reader][short]") {
  const std::filesystem::path compressed_file = "test/data/unit_tests/sample.bed9.gz";
  const std::filesystem::path ptext_file = "test/data/unit_tests/sample.bed9";

  modle::libarchivexx::Reader r1(compressed_file);
  modle::libarchivexx::Reader r2(ptext_file);

  std::string buff1;
  std::string buff2;
  while (r1.getline(buff1) && r2.getline(buff2)) {
    // print_buffers(buff1, buff2);
    CHECK(buff1 == buff2);
  }
  CHECK(!r1.getline(buff1));
  CHECK(!r2.getline(buff1));
}

TEST_CASE("Reader bzip2", "[reader][short]") {
  const std::filesystem::path compressed_file = "test/data/unit_tests/sample.bed9.bz2";
  const std::filesystem::path ptext_file = "test/data/unit_tests/sample.bed9";

  modle::libarchivexx::Reader r1(compressed_file);
  modle::libarchivexx::Reader r2(ptext_file);

  std::string buff1;
  std::string buff2;
  while (r1.getline(buff1) && r2.getline(buff2)) {
    CHECK(buff1 == buff2);
  }
  CHECK(!r1.getline(buff1));
  CHECK(!r2.getline(buff1));
}

TEST_CASE("Reader lz4", "[reader][short]") {
  const std::filesystem::path compressed_file = "test/data/unit_tests/sample.bed9.lz4";
  const std::filesystem::path ptext_file = "test/data/unit_tests/sample.bed9";

  modle::libarchivexx::Reader r1(compressed_file);
  modle::libarchivexx::Reader r2(ptext_file);

  std::string buff1;
  std::string buff2;
  while (r1.getline(buff1) && r2.getline(buff2)) {
    CHECK(buff1 == buff2);
  }
  CHECK(!r1.getline(buff1));
  CHECK(!r2.getline(buff1));
}

TEST_CASE("Reader zstd", "[reader][short]") {
  const std::filesystem::path compressed_file = "test/data/unit_tests/sample.bed9.zst";
  const std::filesystem::path ptext_file = "test/data/unit_tests/sample.bed9";

  modle::libarchivexx::Reader r1(compressed_file);
  modle::libarchivexx::Reader r2(ptext_file);

  std::string buff1;
  std::string buff2;
  while (r1.getline(buff1) && r2.getline(buff2)) {
    CHECK(buff1 == buff2);
  }
  CHECK(!r1.getline(buff1));
  CHECK(!r2.getline(buff1));
}

TEST_CASE("Reader lzma", "[reader][short]") {
  const std::filesystem::path compressed_file = "test/data/unit_tests/sample.bed9.xz";
  const std::filesystem::path ptext_file = "test/data/unit_tests/sample.bed9";

  modle::libarchivexx::Reader r1(compressed_file);
  modle::libarchivexx::Reader r2(ptext_file);

  std::string buff1;
  std::string buff2;
  while (r1.getline(buff1) && r2.getline(buff2)) {
    CHECK(buff1 == buff2);
  }
  CHECK(!r1.getline(buff1));
  CHECK(!r2.getline(buff1));
}

}  // namespace modle::test::libarchivexx
