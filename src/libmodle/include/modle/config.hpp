#pragma once

#include <absl/container/flat_hash_map.h>  // for flat_hash_map, operator!=, raw_hash_set<>::con...
#include <absl/strings/str_cat.h>          // for StrAppend
#include <fmt/format.h>                    // for format, FMT_STRING, print
#include <fmt/ostream.h>                   // for formatbuf<>::int_type

#include <cmath>             // for round
#include <cstdint>           // for uint32_t, uint64_t
#include <cstdio>            // for stderr
#include <filesystem>        // for weakly_canonical, operator<<, path
#include <initializer_list>  // for initializer_list
#include <sstream>           // for basic_stringbuf<>::int_type, basic_stringbuf<>...
#include <string>            // for string, allocator, basic_string
#include <string_view>       // for operator""sv, basic_string_view, string_view
#include <thread>            // for thread
#include <utility>           // for tuple_element<>::type
#include <vector>            // for vector

#include "modle/common.hpp"  // for Bp

namespace modle {
using namespace std::literals::string_view_literals;
class Cli;
class Simulation;

struct Config {
  friend Cli;
  friend Simulation;

 public:
  // clang-format off
  // IO
  std::filesystem::path path_to_chrom_sizes;
  std::filesystem::path path_to_chrom_subranges;
  std::filesystem::path path_to_output_file;
  std::filesystem::path path_to_log_file;
  std::filesystem::path path_to_output_file_w_noise{};
  std::filesystem::path path_to_extr_barriers;
  bool force{false};
  bool write_contacts_for_ko_chroms{false};
  bool write_raw_contacts{true};
  bool write_contacts_with_noise{false};

  // General settings
  bp_t bin_size{1'000};  // NOLINT(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  bp_t fwd_extrusion_speed{std::numeric_limits<bp_t>::max()};
  bp_t rev_extrusion_speed{std::numeric_limits<bp_t>::max()};
  double fwd_extrusion_speed_std{0.05};
  double rev_extrusion_speed_std{fwd_extrusion_speed_std};
  std::size_t ncells{5'000};    // NOLINT(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  std::size_t nthreads{std::thread::hardware_concurrency()};
  bp_t diagonal_width{3'000'000 /* 3 Mbp */};  // NOLINT(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  bp_t simulation_iterations{200};      // NOLINT(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  double target_contact_density{0.0};
  double number_of_lefs_per_mbp;
  bp_t average_lef_lifetime{100'000};  // NOLINT(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  double hard_stall_multiplier{2.0};       // NOLINT(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  double soft_stall_multiplier{0.5};       // NOLINT(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  bool allow_lef_lifetime_extension{true};
  bool skip_burnin{false};
  uint64_t seed{0};

  // Misc probabilities
  double probability_of_extrusion_unit_bypass{0.25};  // NOLINT(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  double lef_fraction_contact_sampling{0.025};  // NOLINT(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  double ctcf_occupied_self_prob{0.875};
  double ctcf_not_occupied_self_prob{0.70};

  // Misc
  bool exclude_chr_wo_extr_barriers{true};
  bool skip_output{false};
  double random_noise_mean{0};
  bp_t random_noise_std{7'500 /* 7.5 Kbp */}; // NOLINT(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)

  int argc;
  char** argv;
  // clang-format on

  [[nodiscard]] inline std::string to_string() const {
    struct cli_tokens {
      std::string name;
      std::string value;
    };

    // clang-format off
    // TODO Update with new params
    const absl::flat_hash_map<std::string_view, std::vector<cli_tokens>> tokens{
    {"Input/Output"sv,
    std::vector<cli_tokens>{
     {"Path to chr. sizes", fmt::format(FMT_STRING("{}"), std::filesystem::weakly_canonical(this->path_to_chrom_sizes))},
     {"Path to extr. barriers", fmt::format(FMT_STRING("{}"), std::filesystem::weakly_canonical(this->path_to_extr_barriers))},
     {"Output directory", fmt::format(FMT_STRING("{}"), std::filesystem::weakly_canonical(this->path_to_output_file))},
     {"Skip output", this->skip_output ? "Yes" : "No"},
     {"Force overwrite", this->force ? "Yes" : "No"}}},
    {"General settings"sv,
    std::vector<cli_tokens>{
     {"Bin size (bp)", fmt::format(FMT_STRING("{}"), this->bin_size)},
     {"Diagonal width (bp)", fmt::format(FMT_STRING("{}"), this->diagonal_width)},
     {"# of iterations", fmt::format(FMT_STRING("{}"), this->simulation_iterations != 0 ? "Disabled"sv : fmt::format("{}", this->simulation_iterations))},
     {"Target contact density", fmt::format(FMT_STRING("{}"), this->target_contact_density != 0 ? fmt::format("{}", this->target_contact_density) : "Disabled"sv)},
     {"Seed", fmt::format(FMT_STRING("{}"), this->seed)},
     {"Contact sampling interval", fmt::format(FMT_STRING("{}"), this->lef_fraction_contact_sampling)}}},
    {"LEF settings"sv,
    std::vector<cli_tokens>{
     {"# of LEFs per Mbp", fmt::format(FMT_STRING("{}"), this->number_of_lefs_per_mbp)},
     {"Avg. LEF lifetime", fmt::format(FMT_STRING("{}"), this->average_lef_lifetime)},
     {"Prob. of LEF bypass", fmt::format(FMT_STRING("{:.4G}"), this->probability_of_extrusion_unit_bypass)},
     {"LEF unloader strength", fmt::format(FMT_STRING("{:.4G}"), this->hard_stall_multiplier)}}},
    {"Extr. barrier settings"sv,
    std::vector<cli_tokens>{
     {"CTCF occupied -> Occupied transition prob.", fmt::format(FMT_STRING("{}"), this->ctcf_occupied_self_prob)},
     {"CTCF not-occupied -> not-occupied transition prob.", fmt::format(FMT_STRING("{}"), this->ctcf_not_occupied_self_prob)}}},
    {"Burn-in phase"sv,
    std::vector<cli_tokens>{
     {"Skip burn-in", this->skip_burnin ? "Yes" : "No"}}}
    };
    // clang-format on

    std::size_t max_column1_length = 0;
    std::size_t max_column2_length = 0;
    for (const auto& [title, options] : tokens) {
      for (const auto& opt : options) {
        max_column1_length =
            opt.name.size() > max_column1_length ? opt.name.size() : max_column1_length;
        max_column2_length =
            opt.value.size() > max_column2_length ? opt.value.size() : max_column2_length;
      }
    }

    constexpr auto col1_constant_padding_length = 9U;
    constexpr auto col2_constant_padding_length = 6U;
    constexpr auto title_constant_padding_length = 11U;

    max_column1_length += col1_constant_padding_length;
    max_column2_length += col2_constant_padding_length;
    const auto max_column_length = max_column1_length + max_column2_length;

    std::string buff = fmt::format(FMT_STRING("{:#^{}}\n##{:^{}}##\n"), "   CONFIG SUMMARY   ",
                                   max_column_length, "", max_column_length - 4);  // NOLINT
    for (const auto& title : {"Input/Output"sv, "General settings"sv, "LEF settings"sv,
                              "Extr. barrier settings"sv, "Burn-in phase"sv}) {
      absl::StrAppend(
          &buff, fmt::format(FMT_STRING("########## {}{:#<{}}\n"), title, " ",
                             max_column_length - title_constant_padding_length - title.size()));
      for (const auto& opt : tokens.at(title)) {
        absl::StrAppend(&buff,
                        fmt::format(FMT_STRING("##    {:<{}}  #  {:<{}}  ##\n"), opt.name,
                                    max_column1_length - col1_constant_padding_length, opt.value,
                                    max_column2_length - col2_constant_padding_length));
      }
    }
    absl::StrAppend(&buff, fmt::format(FMT_STRING("{:#^{}}\n\n"), "   END OF CONFIG SUMMARY   ",
                                       max_column_length));
    return buff;
  }

  inline void print() const { fmt::print(stderr, FMT_STRING("{}\n\n"), this->to_string()); }
};

}  // namespace modle
