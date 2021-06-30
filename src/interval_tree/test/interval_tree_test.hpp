#pragma once

#include <absl/strings/str_split.h>

#include <catch2/catch.hpp>  // for AssertionHandler, operator""_catch_sr
#include <fstream>
#include <string_view>
#include <vector>

#include "modle/common/utils.hpp"
#include "modle/interval_tree.hpp"

namespace modle::test::interval_tree {

struct Record {
  size_t start;
  size_t end;
  std::string data;
};

using IITree_t = IITree<size_t, Record>;

void test_find_overlaps(const IITree_t& tree, size_t start, size_t end,
                        size_t num_expected_overlaps) {
  const auto [overlap_begin, overlap_end] = tree.find_overlaps(start, end);
  CHECK(overlap_end - overlap_begin == num_expected_overlaps);
}

TEST_CASE("Interval tree simple", "[interval-tree][short]") {
  IITree_t tree{};

  tree.insert(0, 10, Record{});  // NOLINT
  tree.insert(5, 15, Record{});  // NOLINT

  tree.make_BST();

  REQUIRE(tree.size() == 2);

  CHECK(tree.overlaps_with(0, 4));
  CHECK(tree.overlaps_with(0, 10));
  CHECK(tree.overlaps_with(5, 15));
  CHECK(tree.overlaps_with(0, 15));
  CHECK(tree.overlaps_with(0, 25));
  CHECK(!tree.overlaps_with(16, 25));

  test_find_overlaps(tree, 0, 4, 1);    // NOLINT
  test_find_overlaps(tree, 0, 10, 2);   // NOLINT
  test_find_overlaps(tree, 5, 15, 2);   // NOLINT
  test_find_overlaps(tree, 0, 15, 2);   // NOLINT
  test_find_overlaps(tree, 0, 25, 2);   // NOLINT
  test_find_overlaps(tree, 16, 25, 0);  // NOLINT
}

TEST_CASE("Interval tree chrX", "[interval-tree][short]") {
  const std::string all_intervals = "test/data/unit_tests/interval_tree_all.bed.gz";
  const std::string subset_intervals = "test/data/unit_tests/interval_tree_subset.bed.gz";
  const std::string complement_intervals = "test/data/unit_tests/interval_tree_complement.bed.gz";

  std::string buff;
  std::vector<std::string_view> toks;
  Record record{};
  IITree_t tree{};

  auto toks_to_record = [&]() {
    REQUIRE(toks.size() == 3UL);
    utils::parse_numeric_or_throw(toks[1], record.start);
    utils::parse_numeric_or_throw(toks[2], record.end);
    record.data = toks[0];
  };

  compressed_io::Reader r(all_intervals);
  REQUIRE(r.is_open());
  while (r.getline(buff)) {
    toks = absl::StrSplit(buff, absl::ByAnyChar("\t "));
    toks_to_record();

    tree.insert(record.start, record.end, record);
  }

  REQUIRE(tree.size() == 2118);  // NOLINT
  tree.make_BST();

  r.open(subset_intervals);
  REQUIRE(r.is_open());

  while (r.getline(buff)) {
    toks = absl::StrSplit(buff, absl::ByAnyChar("\t "));
    toks_to_record();

    CHECK(tree.overlaps_with(record.start, record.end));
    CHECK(tree.count(record.start, record.end) == 1);
    const auto [overlap_begin, overlap_end] = tree.find_overlaps(record.start, record.end);

    CHECK(overlap_begin->start == record.start);
    CHECK(overlap_begin->end == record.end);
    CHECK(overlap_begin->data == record.data);
  }

  r.open(complement_intervals);
  REQUIRE(r.is_open());

  while (r.getline(buff)) {
    toks = absl::StrSplit(buff, absl::ByAnyChar("\t "));
    toks_to_record();

    CHECK(!tree.overlaps_with(record.start, record.end));
    CHECK(tree.count(record.start, record.end) == 0);
    test_find_overlaps(tree, record.start, record.end, 0);
  }

  for (auto it1 = tree.starts_begin(), it2 = tree.starts_begin() + 1; it2 != tree.starts_end();
       ++it1, ++it2) {
    CHECK(*it1 < *it2);
  }
}

}  // namespace modle::test::interval_tree
