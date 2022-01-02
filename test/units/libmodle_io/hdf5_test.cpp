// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "modle/hdf5/hdf5.hpp"

#include <H5Cpp.h>
#include <H5DataSpace.h>   // for DataSpace
#include <H5DcreatProp.h>  // for DSetCreatPropList
#include <H5Exception.h>   // for Exception, H5std_string
#include <H5Spublic.h>     // for H5S_UNLIMITED
#include <fmt/format.h>    // for format

#include <algorithm>                        // for max
#include <boost/filesystem/operations.hpp>  // for create_directories, is_empty, remove, remove_all
#include <boost/filesystem/path.hpp>        // for operator/, path
#include <catch2/catch.hpp>                 // for operator""_catch_sr, AssertionHandler, Source...
#include <limits>                           // for numeric_limits
#include <string>                           // for string, basic_string, allocator, operator==
#include <string_view>                      // for string_view
#include <tuple>                            // for ignore
#include <vector>                           // for vector

#include "H5Ppublic.h"                // for H5F_ACC_TRUNC, H5T_CSET_ASCII, H5T_STR_NULLPAD
#include "modle/common/common.hpp"    // for i64, usize
#include "modle/common/smartdir.hpp"  // for SmartDir

namespace modle::test {
[[maybe_unused]] static const boost::filesystem::path& testdir(bool cleanup_on_exit = true) {
  static const SmartDir dir{cleanup_on_exit};
  return dir();
}
}  // namespace modle::test

namespace modle::test::hdf5 {
using namespace modle::hdf5;

static const usize MAX_STR_LENGTH = 32;

[[nodiscard]] inline H5::StrType init_str_type() {
  auto st = H5::StrType(H5::PredType::C_S1, MAX_STR_LENGTH);
  st.setStrpad(H5T_STR_NULLPAD);
  st.setCset(H5T_CSET_ASCII);

  return st;
}

inline H5::DataSet init_test_str_dataset(H5::H5File& f, std::string_view path = "/test") {
  const hsize_t rank{1};
  const hsize_t chunk_dims{1024};
  const H5std_string fill_var{"\0"};  // NOLINT bugprone-string-literal-with-embedded-nul
  const hsize_t max_dims{H5S_UNLIMITED};
  const auto str_type{init_str_type()};

  H5::DSetCreatPropList cprop{};
  cprop.setChunk(rank, &chunk_dims);
  cprop.setFillValue(str_type, &fill_var);
  cprop.setDeflate(1);

  auto mem_space{H5::DataSpace(rank, &rank, &max_dims)};

  return f.createDataSet(std::string{path.data(), path.size()}, str_type, mem_space, cprop);
}

inline H5::DataSet init_test_int64_dataset(H5::H5File& f, std::string_view path = "/test") {
  const hsize_t rank{1};
  const hsize_t chunk_dims{1024};
  constexpr i64 fill_var{0};
  const hsize_t max_dims{H5S_UNLIMITED};

  H5::DSetCreatPropList cprop{};
  cprop.setChunk(rank, &chunk_dims);
  cprop.setFillValue(H5::PredType::NATIVE_INT64, &fill_var);
  cprop.setDeflate(1);

  auto mem_space{H5::DataSpace(rank, &rank, &max_dims)};

  return f.createDataSet(std::string{path.data(), path.size()}, H5::PredType::NATIVE_INT64,
                         mem_space, cprop);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("read_write_strings HDF5", "[io][hdf5][short]") {
  const auto test_file = testdir() / "rw_strings.hdf5";
  boost::filesystem::create_directories(testdir());

  std::vector<std::string> v{"0", "QWERTY", "ABCDE", "0123456789", "!\"#¤%&/()=?"};

  H5::H5File f(test_file.string(), H5F_ACC_TRUNC);
  auto dataset = init_test_str_dataset(f);

  CHECK(write_strings(v, dataset, init_str_type(), 0) == v.size());
  for (usize i = 0; i < v.size(); ++i) {
    const auto buff = read_str(dataset, i);
    CHECK(buff == v[i]);
  }

  const auto buffv = read_strings(dataset, 0);
  REQUIRE(buffv.size() == v.size());
  for (usize i = 0; i < v.size(); ++i) {
    CHECK(v[i] == buffv[i]);
  }
  boost::filesystem::remove_all(testdir());
  if (const auto& p = testdir().parent_path(); boost::filesystem::is_empty(p)) {
    boost::filesystem::remove(p);
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("read_write_ints HDF5", "[io][hdf5][short]") {
  const auto test_file = testdir() / "rw_ints.hdf5";
  boost::filesystem::create_directories(testdir());
  std::vector<i64> v{(std::numeric_limits<i64>::min)(), -10, 0, 10,
                     (std::numeric_limits<i64>::max)()};
  H5::H5File f(test_file.string(), H5F_ACC_TRUNC);
  auto dataset = init_test_int64_dataset(f);
  i64 buff{};
  for (usize i = 0; i < v.size(); ++i) {
    CHECK(write_number(v[i], dataset, i) == i + 1);
    CHECK(read_number(dataset, buff, i) == i + 1);
    CHECK(buff == v[i]);
  }

  std::vector<i64> buffv(v.size());
  CHECK(write_numbers(v, dataset, 0) == v.size());
  CHECK(read_numbers(dataset, buffv, 0));
  REQUIRE(buffv.size() == v.size());
  for (usize i = 0; i < v.size(); ++i) {
    CHECK(v[i] == buffv[i]);
  }
  boost::filesystem::remove_all(testdir());
  if (const auto& p = testdir().parent_path(); boost::filesystem::is_empty(p)) {
    boost::filesystem::remove(p);
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("has_group HDF5", "[io][hdf5][short]") {
  const auto test_file = testdir() / "has_group.hdf5";
  boost::filesystem::create_directories(testdir());
  H5::H5File f(test_file.string(), H5F_ACC_TRUNC);
  std::ignore = init_test_int64_dataset(f);
  f.createGroup("/g1");
  f.createGroup("/g1/1");

  CHECK(has_group(f, "g1"));
  CHECK(has_group(f, "1", "g1"));
  CHECK(has_group(f, "1", "/g1"));
  CHECK(has_group(f, "1", "/g1/"));
  CHECK(!has_group(f, "g2"));
  CHECK(!has_group(f, "2", "g2"));
  CHECK_THROWS_WITH(has_group(f, "test"), Catch::Matchers::EndsWith("exists but is not a group"));

  boost::filesystem::remove_all(testdir());
  if (const auto& p = testdir().parent_path(); boost::filesystem::is_empty(p)) {
    boost::filesystem::remove(p);
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("check_dataset_type HDF5", "[io][hdf5][short]") {
  const auto test_file = testdir() / "dset_type.hdf5";
  boost::filesystem::create_directories(testdir());
  H5::H5File f(test_file.string(), H5F_ACC_TRUNC);
  f.createGroup("/g1");
  auto int_dataset = init_test_int64_dataset(f);
  auto str_dataset = init_test_str_dataset(f, "/g1/test");
  CHECK(has_dataset(f, "/test"));
  CHECK(!has_dataset(f, "/test2"));

  CHECK(check_dataset_type(int_dataset, H5::PredType::NATIVE_INT64));
  CHECK(check_dataset_type(str_dataset, init_str_type()));
  CHECK(!check_dataset_type(str_dataset, H5::PredType::NATIVE_INT64, false));

  CHECK_THROWS_WITH(check_dataset_type(int_dataset, H5::PredType::NATIVE_INT),
                    Catch::Matchers::Contains("incorrect datasize"));
  CHECK_THROWS_WITH(check_dataset_type(int_dataset, H5::PredType::NATIVE_FLOAT),
                    Catch::Matchers::EndsWith("incorrect datatype"));
  CHECK_THROWS_WITH(check_dataset_type(int_dataset, H5::PredType::NATIVE_UINT64),
                    Catch::Matchers::Contains("incorrect signedness"));

  boost::filesystem::remove_all(testdir());
  if (const auto& p = testdir().parent_path(); boost::filesystem::is_empty(p)) {
    boost::filesystem::remove(p);
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("read_write_string HDF5 - long string", "[io][hdf5][short]") {
  const auto test_file = testdir() / "rw_strings.hdf5";
  boost::filesystem::create_directories(testdir());

  std::string s(MAX_STR_LENGTH + 1, '#');

  H5::H5File f(test_file.string(), H5F_ACC_TRUNC);
  auto dataset = init_test_str_dataset(f);

#ifndef NDEBUG
  CHECK_THROWS_WITH(write_str(s, dataset, init_str_type(), 0),
                    Catch::Contains("string does not fit in the receiving dataset"));
#else
  CHECK(write_str(s, dataset, init_str_type(), 0) == 1);
#endif
  boost::filesystem::remove_all(testdir());
  if (const auto& p = testdir().parent_path(); boost::filesystem::is_empty(p)) {
    boost::filesystem::remove(p);
  }
}

}  // namespace modle::test::hdf5
