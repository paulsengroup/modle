// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "modle/compressed_io/compressed_io.hpp"  // for Reader

#include <boost/filesystem/operations.hpp>  // for remove
#include <boost/filesystem/path.hpp>        // for path, operator/
#include <catch2/catch.hpp>                 // for AssertionHandler, operator""_catch_sr, Source...
#include <fstream>                          // for ifstream, basic_ios, basic_istream, operator<<
#include <iostream>                         // for cerr
#include <string>                           // for operator==, string, basic_string, getline

#include "modle/common/smartdir.hpp"  // for SmartDir

namespace modle::test {
inline const SmartDir testdir{true};  // NOLINT(cert-err58-cpp)
}  // namespace modle::test

namespace modle::test::compressed_io {

[[maybe_unused]] inline void print_buffers(const std::string& buff1, const std::string& buff2) {
  std::cerr << "buff1='" << buff1 << "'\n"
            << "buff2='" << buff2 << "'\n";
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Reader plain", "[io][reader][short]") {
  const boost::filesystem::path ptext_file = "test/data/unit_tests/sample.bed9";

  modle::compressed_io::Reader r(ptext_file);
  std::ifstream fp(ptext_file.string());

  std::string buff1;
  std::string buff2;
  while (r.getline(buff1) && std::getline(fp, buff2)) {
    CHECK(!!fp);
    CHECK(buff1 == buff2);
  }
  CHECK(!r.getline(buff2));
  CHECK(!std::getline(fp, buff2));
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Reader plain sv", "[io][reader][short]") {
  const boost::filesystem::path ptext_file = "test/data/unit_tests/sample.bed9";

  modle::compressed_io::Reader r(ptext_file);
  std::ifstream fp(ptext_file.string());

  std::string_view buff1;
  std::string buff2;
  while (!(buff1 = r.getline()).empty() && std::getline(fp, buff2)) {
    CHECK(!!fp);
    CHECK(buff1 == buff2);
  }
  CHECK(!r.getline(buff2));
  CHECK(!std::getline(fp, buff2));
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Reader gzip", "[io][reader][short]") {
  const boost::filesystem::path compressed_file = "test/data/unit_tests/sample.bed9.gz";
  const boost::filesystem::path ptext_file = "test/data/unit_tests/sample.bed9";

  modle::compressed_io::Reader r1(compressed_file);
  modle::compressed_io::Reader r2(ptext_file);

  std::string buff1;
  std::string buff2;
  while (r1.getline(buff1) && r2.getline(buff2)) {
    // print_buffers(buff1, buff2);
    CHECK(buff1 == buff2);
  }
  CHECK(!r1.getline(buff1));
  CHECK(!r2.getline(buff1));
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Reader bzip2", "[io][reader][short]") {
  const boost::filesystem::path compressed_file = "test/data/unit_tests/sample.bed9.bz2";
  const boost::filesystem::path ptext_file = "test/data/unit_tests/sample.bed9";

  modle::compressed_io::Reader r1(compressed_file);
  modle::compressed_io::Reader r2(ptext_file);

  std::string buff1;
  std::string buff2;
  while (r1.getline(buff1) && r2.getline(buff2)) {
    CHECK(buff1 == buff2);
  }
  CHECK(!r1.getline(buff1));
  CHECK(!r2.getline(buff1));
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Reader lz4", "[io][reader][short]") {
  const boost::filesystem::path compressed_file = "test/data/unit_tests/sample.bed9.lz4";
  const boost::filesystem::path ptext_file = "test/data/unit_tests/sample.bed9";

  modle::compressed_io::Reader r1(compressed_file);
  modle::compressed_io::Reader r2(ptext_file);

  std::string buff1;
  std::string buff2;
  while (r1.getline(buff1) && r2.getline(buff2)) {
    CHECK(buff1 == buff2);
  }
  CHECK(!r1.getline(buff1));
  CHECK(!r2.getline(buff1));
}
/*
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Reader zstd", "[io][reader][short]") {
  const boost::filesystem::path compressed_file = "test/data/unit_tests/sample.bed9.zst";
  const boost::filesystem::path ptext_file = "test/data/unit_tests/sample.bed9";

  modle::compressed_io::Reader r1(compressed_file);
  modle::compressed_io::Reader r2(ptext_file);

  std::string buff1;
  std::string buff2;
  while (r1.getline(buff1) && r2.getline(buff2)) {
    CHECK(buff1 == buff2);
  }
  CHECK(!r1.getline(buff1));
  CHECK(!r2.getline(buff1));
}
 */

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Reader lzma", "[io][reader][short]") {
  const boost::filesystem::path compressed_file = "test/data/unit_tests/sample.bed9.xz";
  const boost::filesystem::path ptext_file = "test/data/unit_tests/sample.bed9";

  modle::compressed_io::Reader r1(compressed_file);
  modle::compressed_io::Reader r2(ptext_file);

  std::string buff1;
  std::string buff2;
  while (r1.getline(buff1) && r2.getline(buff2)) {
    CHECK(buff1 == buff2);
  }
  CHECK(!r1.getline(buff1));
  CHECK(!r2.getline(buff1));
}

inline void test_writer(modle::compressed_io::Reader& r1, modle::compressed_io::Writer& w) {
  std::string buff1;
  std::string buff2;
  while (r1.getline(buff1)) {
    buff1.append("\n");
    w.write(buff1);
  }
  CHECK(!r1.getline(buff1));
  r1.reset();
  w.close();

  modle::compressed_io::Reader r2(w.path());
  while (r1.getline(buff1) && r2.getline(buff2)) {
    CHECK(buff1 == buff2);
  }
  CHECK(!r1.getline(buff1));
  CHECK(!r2.getline(buff2));
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Reader plain - readall", "[io][reader][long]") {
  const boost::filesystem::path ptext_file = "test/data/unit_tests/sample.bed9";

  const auto buff1 = modle::compressed_io::Reader(ptext_file).readall();

  const auto buff2 = [&]() {
    std::ifstream fp(ptext_file.string(), std::ios::ate);
    const auto size = fp.tellg();
    fp.seekg(0, std::ios::beg);
    std::string buff(static_cast<usize>(size), '\0');
    REQUIRE(fp.read(buff.data(), size));
    return buff;
  }();

  CHECK(buff1 == buff2);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Reader plain - empty file", "[io][reader][short]") {
  const auto test_file = testdir() / "empty_file.txt";
  { std::ofstream fp(test_file.string()); }

  modle::compressed_io::Reader r(test_file);
  CHECK(r.eof());
  std::string buff;
  CHECK(!r.getline(buff));
  boost::filesystem::remove(test_file);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Reader plain - one newline", "[io][reader][short]") {
  const auto test_file = testdir() / "one_newline_file.txt";
  {
    std::ofstream fp(test_file.string());
    const auto c = '\n';
    fp.write(&c, 1);
  }

  std::string buff;
  modle::compressed_io::Reader r(test_file);

  CHECK(!r.eof());
  CHECK(r.getline(buff));
  CHECK(buff.empty());
  CHECK(!r.eof());
  CHECK(!r.getline(buff));
  CHECK(r.eof());
  boost::filesystem::remove(test_file);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Reader plain - truncated file", "[io][reader][short]") {
  const auto test_file = testdir() / "truncated_file.txt";
  const std::string buff1{"test"};
  {
    std::ofstream fp(test_file.string());
    fp.write(buff1.data(), static_cast<std::streamsize>(buff1.size()));
  }

  std::string buff2;
  modle::compressed_io::Reader r(test_file);
  CHECK(r.getline(buff2));
  CHECK(buff1 == buff2);
  CHECK(r.eof());
  CHECK(!r.getline(buff2));
  boost::filesystem::remove(test_file);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Writer plain", "[io][writer][short]") {
  const boost::filesystem::path ptext_file = "test/data/unit_tests/sample.bed9";
  const auto tmpout = testdir() / ptext_file.filename();

  modle::compressed_io::Reader r1(ptext_file);
  modle::compressed_io::Writer w(tmpout);

  test_writer(r1, w);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Writer gzip", "[io][writer][short]") {
  const boost::filesystem::path ptext_file = "test/data/unit_tests/sample.bed9";
  const auto tmpout = testdir() / ptext_file.filename();

  modle::compressed_io::Reader r1(ptext_file);
  modle::compressed_io::Writer w(tmpout, modle::compressed_io::Writer::GZIP);

  test_writer(r1, w);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Writer bzip2", "[io][writer][short]") {
  const boost::filesystem::path ptext_file = "test/data/unit_tests/sample.bed9";
  const auto tmpout = testdir() / ptext_file.filename();

  modle::compressed_io::Reader r1(ptext_file);
  modle::compressed_io::Writer w(tmpout, modle::compressed_io::Writer::BZIP2);

  test_writer(r1, w);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Writer lzma", "[io][writer][short]") {
  const boost::filesystem::path ptext_file = "test/data/unit_tests/sample.bed9";
  const auto tmpout = testdir() / ptext_file.filename();

  modle::compressed_io::Reader r1(ptext_file);
  modle::compressed_io::Writer w(tmpout, modle::compressed_io::Writer::LZMA);

  test_writer(r1, w);
}
/*
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Writer zstd", "[io][writer][short]") {
  const boost::filesystem::path ptext_file = "test/data/unit_tests/sample.bed9";
  const auto tmpout = testdir() / ptext_file.filename();

  modle::compressed_io::Reader r1(ptext_file);
  modle::compressed_io::Writer w(tmpout, modle::compressed_io::Writer::ZSTD);

  test_writer(r1, w);
}
 */
}  // namespace modle::test::compressed_io
