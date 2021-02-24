#pragma once

#include <absl/container/flat_hash_map.h>
#include <absl/strings/str_cat.h>
#include <absl/strings/str_split.h>
#include <fmt/format.h>

#include <cstdint>
#include <filesystem>
#include <string_view>
#include <thread>

namespace modle {
using namespace std::literals::string_view_literals;

struct config {
  // clang-format off
  // IO
  std::filesystem::path path_to_chr_sizes;
  std::filesystem::path path_to_chr_subranges;
  std::filesystem::path path_to_output_file;
  std::filesystem::path path_to_log_file;
  std::filesystem::path path_to_output_file_w_noise{};
  std::filesystem::path path_to_extr_barriers_bed;
  bool force{false};
  bool write_contacts_for_ko_chroms{false};
  bool write_raw_contacts{true};
  bool write_contacts_with_noise{false};

  // General settings
  uint32_t bin_size{1'000};  // NOLINT(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  uint32_t nthreads{std::thread::hardware_concurrency()};
  uint32_t diagonal_width{3'000'000 /* 3 Mbp */};  // NOLINT(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  uint32_t simulation_iterations{15'000'000};  // NOLINT(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  double target_contact_density{0.0};
  uint32_t number_of_lefs;
  uint32_t average_lef_processivity{100'000};  // NOLINT(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  double lef_unloader_strength{1.0};
  bool skip_burnin{false};
  uint64_t seed{0};

  // Misc probabilities
  double probability_of_lef_rebind{1.0};
  double probability_of_extrusion_unit_bypass{0.25};  // NOLINT(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
  double probability_of_extrusion_barrier_block{0};
  uint32_t contact_sampling_interval{1000};  // NOLINT(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)

  // Burn-in
  uint32_t min_n_of_burnin_rounds{0};
  uint32_t min_n_of_loops_per_lef{1};

  // Misc
  uint32_t number_of_randomly_gen_extr_barriers{0};
  bool exclude_chr_wo_extr_barriers{true};
  bool randomize_contact_sampling_interval{false};
  bool skip_output{false};
  double random_noise_mean{0};
  double random_noise_std{7'500 /* 7.5 Kbp */}; // NOLINT(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)

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
     {"Path to chr. sizes", fmt::format(FMT_STRING("{}"), std::filesystem::weakly_canonical(this->path_to_chr_sizes))},
     {"Path to extr. barriers", fmt::format(FMT_STRING("{}"), std::filesystem::weakly_canonical(this->path_to_extr_barriers_bed))},
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
     {"Randomize contact sampling interval", this->randomize_contact_sampling_interval ? "Yes" : "No"},
     {"Contact sampling interval", fmt::format(FMT_STRING("{}"), this->contact_sampling_interval)}}},
    {"LEF settings"sv,
    std::vector<cli_tokens>{
     {"# of LEFs", fmt::format(FMT_STRING("{}"), this->number_of_lefs)},
     {"Avg. LEF processivity", fmt::format(FMT_STRING("{}"), this->average_lef_processivity)},
     {"Prob. of LEF re-bind", fmt::format(FMT_STRING("{:.4G}"), this->probability_of_lef_rebind)},
     {"Prob. of LEF bypass", fmt::format(FMT_STRING("{:.4G}"), this->probability_of_extrusion_unit_bypass)},
     {"LEF unloader strength", fmt::format(FMT_STRING("{:.4G}"), this->lef_unloader_strength)}}},
    {"Extr. barrier settings"sv,
    std::vector<cli_tokens>{
     {"# of randomly gen extr. barriers", fmt::format(FMT_STRING("{}"), this->number_of_randomly_gen_extr_barriers)},
     {"Prob. of block", fmt::format(FMT_STRING("{}"), this->probability_of_extrusion_barrier_block)}}},
    {"Burn-in phase"sv,
    std::vector<cli_tokens>{
     {"Skip burn-in", this->skip_burnin ? "Yes" : "No"},
     {"Min. # of burn-in rounds", this->min_n_of_burnin_rounds != 0 ? fmt::format(FMT_STRING("{}"), this->min_n_of_burnin_rounds) : "Disabled"},
     {"Min. # of loops per LEF", fmt::format(FMT_STRING("{}"), this->min_n_of_loops_per_lef)}}}
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
