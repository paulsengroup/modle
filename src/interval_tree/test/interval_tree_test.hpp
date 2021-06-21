#pragma once

#include <absl/strings/str_split.h>

#include <catch2/catch.hpp>  // for AssertionHandler, operator""_catch_sr
#include <fstream>
#include <string_view>
#include <vector>

#include "modle/common/utils.hpp"
#include "modle/interval_tree.hpp"

namespace modle::test::interval_tree {
using IITree_t = IITree<size_t, std::string>;

struct Record {
  size_t start;
  size_t end;
  std::string data;
};

TEST_CASE("Interval tree simple", "[interval-tree][short]") {
  IITree_t tree{};

  tree.insert(0, 10, "");
  tree.insert(5, 15, "");

  tree.make_BST();

  REQUIRE(tree.size() == 2);

  std::vector<size_t> buff;

  CHECK(tree.overlaps_with(0, 4));
  CHECK(tree.overlaps_with(0, 10));
  CHECK(tree.overlaps_with(5, 15));
  CHECK(tree.overlaps_with(0, 15));
  CHECK(tree.overlaps_with(0, 25));
  CHECK(!tree.overlaps_with(16, 25));

  CHECK(tree.find_overlaps(0, 4, buff));
  CHECK(buff.size() == 1);
  CHECK(tree.find_overlaps(0, 10, buff));
  CHECK(buff.size() == 2);
  CHECK(tree.find_overlaps(5, 15, buff));
  CHECK(buff.size() == 2);
  CHECK(tree.find_overlaps(0, 15, buff));
  CHECK(buff.size() == 2);
  CHECK(tree.find_overlaps(0, 25, buff));
  CHECK(buff.size() == 2);
  CHECK(!tree.find_overlaps(16, 25, buff));
  CHECK(buff.empty());
}

TEST_CASE("Interval tree chrX", "[interval-tree][short]") {
  const std::string all_intervals = "test/data/unit_tests/interval_tree_all.bed";
  const std::string subset_intervals = "test/data/unit_tests/interval_tree_subset.bed";
  const std::string complement_intervals = "test/data/unit_tests/interval_tree_complement.bed";

  std::string buff;
  std::vector<std::string_view> toks;
  Record record{};
  IITree_t tree{};
  std::vector<size_t> idx_buff;

  auto toks_to_record = [&]() {
    REQUIRE(toks.size() == 3UL);
    utils::parse_numeric_or_throw(toks[1], record.start);
    utils::parse_numeric_or_throw(toks[2], record.end);
    record.data = toks[0];
  };

  {
    std::ifstream f(all_intervals);
    REQUIRE((f.is_open() && f.good()));
    while (std::getline(f, buff)) {
      toks = absl::StrSplit(buff, absl::ByAnyChar("\t "));
      toks_to_record();

      tree.insert(record.start, record.end, record.data);
    }

    REQUIRE(f.eof());
    REQUIRE(tree.size() == 2118);  // NOLINT
  }
  tree.make_BST();

  {
    std::ifstream f(subset_intervals);
    REQUIRE((f.is_open() && f.good()));

    while (std::getline(f, buff)) {
      toks = absl::StrSplit(buff, absl::ByAnyChar("\t "));
      toks_to_record();

      CHECK(tree.overlaps_with(record.start, record.end));
      CHECK(tree.find_overlaps(record.start, record.end, idx_buff));
      CHECK(idx_buff.size() == 1);

      CHECK(tree.get_overlap_start(idx_buff.front()) == record.start);
      CHECK(tree.get_overlap_end(idx_buff.front()) == record.end);
      CHECK(tree.get_overlap_data(idx_buff.front()) == record.data);
    }

    REQUIRE(f.eof());
  }

  std::ifstream f(complement_intervals);
  REQUIRE((f.is_open() && f.good()));

  while (std::getline(f, buff)) {
    toks = absl::StrSplit(buff, absl::ByAnyChar("\t "));
    toks_to_record();

    CHECK(!tree.overlaps_with(record.start, record.end));
    CHECK(!tree.find_overlaps(record.start, record.end, idx_buff));
    CHECK(idx_buff.empty());
  }

  REQUIRE(f.eof());

  for (auto it1 = tree.starts_begin(), it2 = tree.starts_begin() + 1; it2 != tree.starts_end();
       ++it1, ++it2) {
    CHECK(*it1 < *it2);
  }
}

}  // namespace modle::test::interval_tree
