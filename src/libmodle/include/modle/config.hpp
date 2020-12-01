#pragma once

#include <absl/container/flat_hash_map.h>
#include <absl/strings/str_cat.h>
#include <absl/strings/str_split.h>
#include <fmt/printf.h>

#include <cstdint>
#include <filesystem>
#include <string_view>

namespace modle {
using namespace std::literals::string_view_literals;
struct config {
  std::string_view path_to_chr_sizes;
  std::string_view path_to_extr_barriers_bed;
  std::string_view output_dir;
  bool force{false};
  uint32_t bin_size{1'000};
  uint32_t diagonal_width{3'000'000};
  uint32_t simulation_iterations{15'000'000};
  double target_contact_density{0.0};
  uint32_t number_of_randomly_gen_extr_barriers{0};
  uint32_t number_of_lefs;
  uint32_t average_lef_processivity{100'000};
  double probability_of_extrusion_barrier_block{0};
  double probability_of_lef_rebind{1.0};
  double probability_of_extrusion_unit_bypass{0.25};
  uint64_t seed{0};
  bool skip_burnin{false};
  uint32_t min_n_of_burnin_rounds{0};
  uint32_t min_n_of_loops_per_lef{1};
  double lef_unloader_strength{1};
  uint32_t contact_sampling_interval{1000};
  bool randomize_contact_sampling_interval{false};
  bool skip_output{false};

  int argc;
  char** argv;

  [[nodiscard]] std::string to_string() const {
    struct cli_tokens {
      std::string name;
      std::string value;
    };

    // clang-format off
    // TODO: Remove FMT_STRING macro from trivial format strings
    absl::flat_hash_map<std::string_view, std::vector<cli_tokens>> tokens{
    {"Input/Output"sv,
    std::vector<cli_tokens>{
     {"Path to chr. sizes", fmt::format(FMT_STRING("{}"), std::filesystem::weakly_canonical(this->path_to_chr_sizes))},
     {"Path to extr. barriers", fmt::format(FMT_STRING("{}"), std::filesystem::weakly_canonical(this->path_to_extr_barriers_bed))},
     {"Output directory", fmt::format(FMT_STRING("{}"), std::filesystem::weakly_canonical(this->output_dir))},
     {"Skip output", this->skip_output ? "Yes" : "No"},
     {"Force overwrite", this->force ? "Yes" : "No"}}},
    {"General settings"sv,
    std::vector<cli_tokens>{
     {"Bin size (bp)", fmt::format(FMT_STRING("{}"), this->bin_size)},
     {"Diagonal width (bp)", fmt::format(FMT_STRING("{}"), this->diagonal_width)},
     {"# of iterations", fmt::format(FMT_STRING("{}"), this->target_contact_density != 0 ? "Disabled"sv : fmt::format("{}", this->simulation_iterations))},
     {"Target contact density", fmt::format(FMT_STRING("{}"), this->target_contact_density != 0 ? fmt::format("{}", this->simulation_iterations) : "Disabled"sv)},
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
    for (auto& [title, options] : tokens) {
      for (auto& opt : options) {
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
      for (auto& opt : tokens.at(title)) {
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

  void print() const { fmt::fprintf(stderr, FMT_STRING("%s\n\n"), this->to_string()); }
};

}  // namespace modle