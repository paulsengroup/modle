#pragma once

#include <algorithm>
#include <boost/process.hpp>
#include <filesystem>
#include <random>
#include <string>
#include <string_view>

#include "absl/strings/match.h"
#include "absl/strings/str_cat.h"
#include "modle_tools/cli.hpp"
#include "modle_tools/utils.hpp"

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
    std::generate(rand_name.begin() + base_name.size(), rand_name.end(),
                  [&]() { return chars[dist(rnd_eng)]; });
  } while (base_name.empty() || std::filesystem::exists(absl::StrCat(base_name, rand_name)));

  return rand_name;
}

inline std::string init_juicer_tools_argv(const modle::tools::config& c) {
  std::string argv = c.path_to_juicer_tools;
  if (absl::EndsWith(argv, ".jar")) {
    // TODO: Check that java >= 1.7
    auto java = boost::process::search_path("java").string();
    if (java.empty())
      throw std::runtime_error(
          "--path-to-juicer-tools points to a jar file, but we were unable to find java in "
          "your "
          "path");
    argv = absl::StrFormat("java -Xms512m -Xmx%.0fm -jar %s", c.juicer_tools_mem / 1e6,
                           c.path_to_juicer_tools);
  }
  return argv;
}

}  // namespace modle::tools::utils
