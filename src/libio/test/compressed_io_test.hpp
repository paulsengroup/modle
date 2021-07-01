#pragma once

#include <catch2/catch.hpp>  // for AssertionHandler, operator""_catch_sr, SourceLineInfo
#include <boost/filesystem/path.hpp>        // for exists
#include <fstream>           // for ifstrea, getline
#include <iostream>
#include <string>  // for string, basic_string, operator==, char_traits, stoull

#include "modle/common/smartdir.hpp"  // IWYU pragma: keep
#include "modle/compressed_io.hpp"    // for Reader

namespace modle::test {
const extern SmartDir testdir;  // NOLINT
}  // namespace modle::test

namespace modle::test::compressed_io {

[[maybe_unused]] inline void print_buffers(const std::string& buff1, const std::string& buff2) {
  std::cerr << "buff1='" << buff1 << "'\n"
            << "buff2='" << buff2 << "'\n";
}

TEST_CASE("Reader plain", "[io][reader][short]") {
  const boost::filesystem::path ptext_file = "test/data/unit_tests/sample.bed9";

  modle::compressed_io::Reader r(ptext_file);
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

TEST_CASE("Reader plain sv", "[io][reader][short]") {
  const boost::filesystem::path ptext_file = "test/data/unit_tests/sample.bed9";

  modle::compressed_io::Reader r(ptext_file);
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

TEST_CASE("Writer plain", "[io][writer][short]") {
  const boost::filesystem::path ptext_file = "test/data/unit_tests/sample.bed9";
  const auto tmpout = testdir() / ptext_file.filename();

  modle::compressed_io::Reader r1(ptext_file);
  modle::compressed_io::Writer w(tmpout);

  test_writer(r1, w);
}

TEST_CASE("Writer gzip", "[io][writer][short]") {
  const boost::filesystem::path ptext_file = "test/data/unit_tests/sample.bed9";
  const auto tmpout = testdir() / ptext_file.filename();

  modle::compressed_io::Reader r1(ptext_file);
  modle::compressed_io::Writer w(tmpout, modle::compressed_io::Writer::GZIP);

  test_writer(r1, w);
}

TEST_CASE("Writer bzip2", "[io][writer][short]") {
  const boost::filesystem::path ptext_file = "test/data/unit_tests/sample.bed9";
  const auto tmpout = testdir() / ptext_file.filename();

  modle::compressed_io::Reader r1(ptext_file);
  modle::compressed_io::Writer w(tmpout, modle::compressed_io::Writer::BZIP2);

  test_writer(r1, w);
}

TEST_CASE("Writer lzma", "[io][writer][short]") {
  const boost::filesystem::path ptext_file = "test/data/unit_tests/sample.bed9";
  const auto tmpout = testdir() / ptext_file.filename();

  modle::compressed_io::Reader r1(ptext_file);
  modle::compressed_io::Writer w(tmpout, modle::compressed_io::Writer::LZMA);

  test_writer(r1, w);
}

TEST_CASE("Writer zstd", "[io][writer][short]") {
  const boost::filesystem::path ptext_file = "test/data/unit_tests/sample.bed9";
  const auto tmpout = testdir() / ptext_file.filename();

  modle::compressed_io::Reader r1(ptext_file);
  modle::compressed_io::Writer w(tmpout, modle::compressed_io::Writer::ZSTD);

  test_writer(r1, w);
}
}  // namespace modle::test::compressed_io
