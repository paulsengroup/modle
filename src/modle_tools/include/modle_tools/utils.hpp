#pragma once

#include <absl/strings/match.h>
#include <absl/strings/str_cat.h>
#include <absl/strings/str_split.h>
#include <absl/strings/strip.h>
#include <fmt/format.h>

#include <algorithm>
#include <boost/process/child.hpp>
#include <boost/process/pipe.hpp>
#include <filesystem>
#include <modle_tools/utils.hpp>
#include <random>
#include <string>
#include <string_view>

#include "modle/suppress_compiler_warnings.hpp"
#include "modle/utils.hpp"

namespace modle::tools::utils {

inline std::string generate_random_path(std::string_view base_name = "",
                                        std::string_view suffix = "") {
  static constexpr auto chars =
      "0123456789"
      "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
      "abcdefghijklmnopqrstuvwxyz";
  std::random_device rndev;
  std::mt19937 rnd_eng{rndev()};
  auto dist = std::uniform_int_distribution{{}, std::strlen(chars) - 1};
  std::string rand_path = absl::StrCat(base_name, std::string(16, '\0'), suffix);
  do {
    DISABLE_WARNING_PUSH
    DISABLE_WARNING_SIGN_CONVERSION
    std::generate(rand_path.begin() + base_name.size(), rand_path.end() - suffix.size(),
                  [&]() { return chars[dist(rnd_eng)]; });
    DISABLE_WARNING_POP
  } while (base_name.empty() || std::filesystem::exists(rand_path));

  return rand_path;
}

[[nodiscard]] inline std::vector<uint8_t> parse_version(std::vector<std::string_view> toks) {
  std::vector<uint8_t> ver;
  for (const auto& tok : toks) {
    uint8_t& n = ver.emplace_back();
    modle::utils::parse_numeric_or_throw(tok, n);
  }

  return ver;
}

[[nodiscard]] inline std::vector<uint8_t> detect_juicer_tools_version(
    std::string_view template_argv) {
  boost::process::ipstream juicer_tools_stdout;

  try {
    boost::process::child juicer_tools(
        fmt::format("{} -V", template_argv), boost::process::std_in.close(),
        boost::process::std_err.close(), boost::process::std_out > juicer_tools_stdout);
    std::string buff;
    while (juicer_tools.running() && std::getline(juicer_tools_stdout, buff)) {
      if (absl::StartsWithIgnoreCase(buff, "Juicer Tools Version")) {
        buff = absl::StripPrefix(buff, "Juicer Tools Version ");
        return parse_version(absl::StrSplit(buff, '.'));
      }
    }
  } catch (const std::runtime_error& err) {
    throw std::runtime_error(fmt::format(
        "An error occurred while trying to detect Juicer Tools version: {}", err.what()));
  }
  fmt::print(
      stderr,
      "WARNING: Unable to detect Juicer Tools version. Assuming the latest version is used\n.");
  return {255, 255, 255};
}

[[nodiscard]] inline bool juicer_tools_version_is_greater(const std::vector<uint8_t>& v1,
                                                          const std::vector<uint8_t>& v2) {
  for (auto i = 0UL; i < std::min(v1.size(), v2.size()); ++i) {
    if (v1[i] > v2[i]) {
      return true;
    }
  }
  return v2.size() > v1.size();
}

[[nodiscard]] inline bool juicer_tools_version_is_greater_or_equal(const std::vector<uint8_t>& v1,
                                                                   const std::vector<uint8_t>& v2) {
  for (auto i = 0UL; i < std::min(v1.size(), v2.size()); ++i) {
    if (v1[i] >= v2[i]) {
      return true;
    }
  }
  return v2.size() > v1.size();
}

[[nodiscard]] inline bool juicer_tools_version_is_equal(const std::vector<uint8_t>& v1,
                                                        const std::vector<uint8_t>& v2) {
  if (v1.size() != v2.size()) {
    return false;
  }

  for (auto i = 0UL; i < v1.size(); ++i) {
    if (v1[i] != v2[i]) {
      return false;
    }
  }
  return true;
}

}  // namespace modle::tools::utils
