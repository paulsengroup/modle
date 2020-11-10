#pragma once

#include <algorithm>
#include <filesystem>
#include <random>
#include <string>
#include <string_view>

#include "absl/strings/str_cat.h"

namespace modle::tools::utils {

std::string generate_random_path_name(std::string_view base_name = "") {
  static constexpr auto chars =
      "0123456789"
      "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
      "abcdefghijklmnopqrstuvwxyz";
  std::random_device rndev;
  std::mt19937 rnd_eng{rndev()};
  auto dist = std::uniform_int_distribution{{}, std::strlen(chars) - 1};
  auto rand_name = base_name.data() + std::string(16, '\0');
  do {
    std::generate(rand_name.begin() + base_name.size(), rand_name.end(),
                  [&]() { return chars[dist(rnd_eng)]; });
  } while (base_name.empty() || std::filesystem::exists(absl::StrCat(base_name, rand_name)));

  return rand_name;
}
}  // namespace modle::tools::utils
