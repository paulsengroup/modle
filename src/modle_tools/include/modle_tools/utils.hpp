#pragma once

#include <absl/strings/match.h>
#include <absl/strings/str_cat.h>

#include <algorithm>
#include <boost/process.hpp>
#include <filesystem>
#include <modle_tools/cli.hpp>
#include <modle_tools/utils.hpp>
#include <random>
#include <string>
#include <string_view>

namespace modle::tools::utils {

inline std::string generate_random_path_name(std::string_view base_name) {
  static constexpr auto chars =
      "0123456789"
      "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
      "abcdefghijklmnopqrstuvwxyz";
  std::random_device rndev;
  std::mt19937 rnd_eng{rndev()};
  auto dist = std::uniform_int_distribution{{}, std::strlen(chars) - 1};
  auto rand_name = base_name.data() + std::string(16, '\0');
  do {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-conversion"
    std::generate(rand_name.begin() + base_name.size(), rand_name.end(),
                  [&]() { return chars[dist(rnd_eng)]; });
#pragma GCC diagnostic pop
  } while (base_name.empty() || std::filesystem::exists(absl::StrCat(base_name, rand_name)));

  return rand_name;
}

}  // namespace modle::tools::utils
