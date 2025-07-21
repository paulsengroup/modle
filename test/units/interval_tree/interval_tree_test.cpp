// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "modle/interval_tree.hpp"

#include <absl/strings/str_split.h>

#include <catch2/catch_test_macros.hpp>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

#include "modle/common/numeric_utils.hpp"
#include "modle/compressed_io/compressed_io.hpp"

namespace modle::interval_tree::test {

[[maybe_unused]] static const std::filesystem::path& data_dir() {
  static const std::filesystem::path data_dir{"test/data/unit_tests"};
  return data_dir;
}

struct Record {
  usize start{};
  usize end{};
  std::string data{};

  Record() = default;

  explicit Record(std::string_view buff) {
    const std::vector<std::string_view> toks = absl::StrSplit(buff, absl::ByAnyChar("\t "));
    REQUIRE(toks.size() == 3UL);
    utils::parse_numeric_or_throw(toks[1], this->start);
    utils::parse_numeric_or_throw(toks[2], this->end);
    this->data = toks[0];
  }
};

template <class IITreeT>
static void test_find_overlaps(const IITreeT& tree, usize start, usize end,
                               usize num_expected_overlaps) {
  const auto [overlap_begin, overlap_end] = tree.find_overlaps(start, end);
  CHECK(static_cast<usize>(overlap_end - overlap_begin) == num_expected_overlaps);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Interval tree simple", "[interval-tree][short]") {
  IITree<usize, Record> tree{};

  tree.insert(0, 10, Record{});
  tree.insert(5, 15, Record{});

  tree.make_BST();

  REQUIRE(tree.size() == 2);

  CHECK(tree.overlaps_with(0, 4));
  CHECK(tree.overlaps_with(0, 10));
  CHECK(tree.overlaps_with(5, 15));
  CHECK(tree.overlaps_with(0, 15));
  CHECK(tree.overlaps_with(0, 25));
  CHECK(!tree.overlaps_with(16, 25));

  test_find_overlaps(tree, 0, 4, 1);
  test_find_overlaps(tree, 0, 10, 2);
  test_find_overlaps(tree, 5, 15, 2);
  test_find_overlaps(tree, 0, 15, 2);
  test_find_overlaps(tree, 0, 25, 2);
  test_find_overlaps(tree, 16, 25, 0);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Interval tree chrX", "[interval-tree][short]") {
  const auto parent = data_dir() / "genomic_intervals";
  const auto all_intervals = parent / "test_intervals.bed3.xz";
  const auto interval_queries_hits = parent / "test_interval_queries_001.bed.xz";
  const auto interval_queries_miss = parent / "test_interval_queries_002.bed.xz";

  std::string buff;

  // Construct interval tree
  const auto tree = [&]() {
    IITree<usize, Record> tree_{};

    compressed_io::Reader r(all_intervals);
    REQUIRE(r.is_open());
    while (r.getline(buff)) {
      const Record record(buff);
      tree_.insert(record.start, record.end, record);
    }

    REQUIRE(tree_.size() == 2118);
    tree_.make_BST();
    return tree_;
  }();

  SECTION("positive queries") {
    compressed_io::Reader r(interval_queries_hits);
    REQUIRE(r.is_open());

    while (r.getline(buff)) {
      const Record query(buff);
      CHECK(tree.overlaps_with(query.start, query.end));
      CHECK(tree.count(query.start, query.end) == 1);
      const auto overlap_begin = tree.find_overlaps(query.start, query.end).first;

      CHECK(overlap_begin->start == query.start);
      CHECK(overlap_begin->end == query.end);
      CHECK(overlap_begin->data == query.data);
    }
  }

  SECTION("negative queries") {
    compressed_io::Reader r(interval_queries_miss);
    REQUIRE(r.is_open());

    while (r.getline(buff)) {
      const Record query(buff);

      CHECK(!tree.overlaps_with(query.start, query.end));
      CHECK(tree.count(query.start, query.end) == 0);
      test_find_overlaps(tree, query.start, query.end, 0);
    }

    for (auto it1 = tree.starts_begin(), it2 = tree.starts_begin() + 1; it2 != tree.starts_end();
         ++it1, ++it2) {
      CHECK(*it1 < *it2);
    }
  }
}

}  // namespace modle::interval_tree::test
