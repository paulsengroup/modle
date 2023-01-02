// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "modle/compressed_io/compressed_io.hpp"  // for Reader

#include <catch2/catch_test_macros.hpp>
#include <filesystem>  // for path, operator/
#include <fstream>     // for ifstream, basic_ios, basic_istream, operator<<
#include <iostream>    // for cerr
#include <string>      // for operator==, string, basic_string, getline

#include "modle/test/self_deleting_folder.hpp"  // for SelfDeletingFolder

namespace modle::test {
inline const SelfDeletingFolder testdir{true};  // NOLINT(cert-err58-cpp)
}  // namespace modle::test

namespace modle::compressed_io::test {

constexpr auto& testdir = modle::test::testdir;

[[maybe_unused]] static const std::filesystem::path& data_dir() {
  static const std::filesystem::path data_dir{"test/data/unit_tests"};
  return data_dir;
}

[[maybe_unused]] static void print_buffers(const std::string& buff1, const std::string& buff2) {
  std::cerr << "buff1='" << buff1 << "'\n"
            << "buff2='" << buff2 << "'\n";
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Reader plain", "[io][reader][short]") {
  const auto ptext_file = data_dir() / "genomic_intervals" / "intervals.bed9";

  Reader r(ptext_file);
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
  const auto ptext_file = data_dir() / "genomic_intervals" / "intervals.bed9";

  Reader r(ptext_file);
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
  const auto compressed_file = data_dir() / "genomic_intervals" / "intervals.bed9.gz";
  const auto ptext_file = data_dir() / "genomic_intervals" / "intervals.bed9";

  Reader r1(compressed_file);
  Reader r2(ptext_file);

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
  const auto compressed_file = data_dir() / "genomic_intervals" / "intervals.bed9.bz2";
  const auto ptext_file = data_dir() / "genomic_intervals" / "intervals.bed9";

  Reader r1(compressed_file);
  Reader r2(ptext_file);

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
  const auto compressed_file = data_dir() / "genomic_intervals" / "intervals.bed9.lz4";
  const auto ptext_file = data_dir() / "genomic_intervals" / "intervals.bed9";

  Reader r1(compressed_file);
  Reader r2(ptext_file);

  std::string buff1;
  std::string buff2;
  while (r1.getline(buff1) && r2.getline(buff2)) {
    CHECK(buff1 == buff2);
  }
  CHECK(!r1.getline(buff1));
  CHECK(!r2.getline(buff1));
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Reader zstd", "[io][reader][short]") {
  const auto compressed_file = data_dir() / "genomic_intervals" / "intervals.bed9.zst";
  const auto ptext_file = data_dir() / "genomic_intervals" / "intervals.bed9";

  Reader r1(compressed_file);
  Reader r2(ptext_file);

  std::string buff1;
  std::string buff2;
  while (r1.getline(buff1) && r2.getline(buff2)) {
    CHECK(buff1 == buff2);
  }
  CHECK(!r1.getline(buff1));
  CHECK(!r2.getline(buff1));
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Reader lzma", "[io][reader][short]") {
  const auto compressed_file = data_dir() / "genomic_intervals" / "intervals.bed9.xz";
  const auto ptext_file = data_dir() / "genomic_intervals" / "intervals.bed9";

  Reader r1(compressed_file);
  Reader r2(ptext_file);

  std::string buff1;
  std::string buff2;
  while (r1.getline(buff1) && r2.getline(buff2)) {
    CHECK(buff1 == buff2);
  }
  CHECK(!r1.getline(buff1));
  CHECK(!r2.getline(buff1));
}

static void compare_compressed_writer_with_reader(Reader& r1, Writer& w) {
  std::string buff1;
  std::string buff2;
  while (r1.getline(buff1)) {
    buff1.append("\n");
    w.write(buff1);
  }
  CHECK(!r1.getline(buff1));
  r1.reset();
  w.close();

  Reader r2(w.path());
  while (r1.getline(buff1) && r2.getline(buff2)) {
    CHECK(buff1 == buff2);
  }
  CHECK(!r1.getline(buff1));
  CHECK(!r2.getline(buff2));
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Reader plain - readall", "[io][reader][long]") {
  const auto ptext_file = data_dir() / "genomic_intervals" / "intervals.bed9";

  const auto buff1 = Reader(ptext_file).readall();

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

  Reader r(test_file);
  CHECK(r.eof());
  std::string buff;
  CHECK(!r.getline(buff));
  std::filesystem::remove(test_file);
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
  Reader r(test_file);

  CHECK(!r.eof());
  CHECK(r.getline(buff));
  CHECK(buff.empty());
  CHECK(!r.eof());
  CHECK(!r.getline(buff));
  CHECK(r.eof());
  std::filesystem::remove(test_file);
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
  Reader r(test_file);
  CHECK(r.getline(buff2));
  CHECK(buff1 == buff2);
  CHECK(r.eof());
  CHECK(!r.getline(buff2));
  std::filesystem::remove(test_file);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Writer plain", "[io][writer][short]") {
  const auto ptext_file = data_dir() / "genomic_intervals" / "intervals.bed9";
  const auto tmpout = testdir() / ptext_file.filename();

  Reader r1(ptext_file);
  Writer w(tmpout);

  compare_compressed_writer_with_reader(r1, w);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Writer gzip", "[io][writer][short]") {
  const auto ptext_file = data_dir() / "genomic_intervals" / "intervals.bed9";
  const auto tmpout = testdir() / ptext_file.filename();

  Reader r1(ptext_file);
  Writer w(tmpout, Writer::GZIP);

  compare_compressed_writer_with_reader(r1, w);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Writer bzip2", "[io][writer][short]") {
  const auto ptext_file = data_dir() / "genomic_intervals" / "intervals.bed9";
  const auto tmpout = testdir() / ptext_file.filename();

  Reader r1(ptext_file);
  Writer w(tmpout, Writer::BZIP2);

  compare_compressed_writer_with_reader(r1, w);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Writer lzma", "[io][writer][short]") {
  const auto ptext_file = data_dir() / "genomic_intervals" / "intervals.bed9";
  const auto tmpout = testdir() / ptext_file.filename();

  Reader r1(ptext_file);
  Writer w(tmpout, Writer::LZMA);

  compare_compressed_writer_with_reader(r1, w);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Writer zstd", "[io][writer][short]") {
  const auto ptext_file = data_dir() / "genomic_intervals" / "intervals.bed9";
  const auto tmpout = testdir() / ptext_file.filename();

  Reader r1(ptext_file);
  Writer w(tmpout, Writer::ZSTD);

  compare_compressed_writer_with_reader(r1, w);
}
}  // namespace modle::compressed_io::test
