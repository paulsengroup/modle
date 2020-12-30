#include "modle_tools/eval.hpp"

#include <absl/algorithm/container.h>
#include <absl/container/btree_set.h>
#include <absl/strings/match.h>

#include "modle/cooler.hpp"
#include "modle/utils.hpp"

namespace modle::tools {

std::vector<std::pair<std::string, int64_t>> select_chromosomes_for_eval(
    std::string_view path_to_cooler1, std::string_view path_to_cooler2) {
  std::vector<std::string> str_buff;
  std::vector<int64_t> int_buff;

  auto build_chr_set = [&](std::string_view path_to_cooler) {
    try {
      auto f = cooler::open_for_reading(path_to_cooler);
      const auto nchroms = static_cast<std::size_t>(cooler::read_attribute_int(f, "nchroms"));
      str_buff.resize(nchroms);
      int_buff.resize(nchroms);

      cooler::read_vect_of_str(f, "chroms/name", str_buff, 0);
      cooler::read_vect_of_int(f, "chroms/length", int_buff, 0);

      assert(str_buff.size() == nchroms);
      assert(int_buff.size() == nchroms);

      absl::btree_set<std::pair<std::string, int64_t>> chr_set;

      for (auto i = 0UL; i < nchroms; ++i) {
        chr_set.emplace(str_buff[i], int_buff[i]);
      }
      return chr_set;

    } catch (const std::runtime_error& e) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("An error occurred while reading file '{}': {}"), path_to_cooler, e.what()));
    }
  };

  std::vector<std::pair<std::string, int64_t>> chr_intersection;
  const auto chr_set1 = build_chr_set(path_to_cooler1);
  const auto chr_set2 = build_chr_set(path_to_cooler2);

  absl::c_set_intersection(
      chr_set1, chr_set2, std::back_inserter(chr_intersection),
      [](const auto& c1, const auto& c2) { return utils::chr_less_than_operator(c1, c2); });

  return chr_intersection;
}

}  // namespace modle::tools
