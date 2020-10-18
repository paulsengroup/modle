#pragma once

#include <cstdint>
#include <execution>
#include <filesystem>
#include <string_view>

#include "absl/strings/str_format.h"
#include "absl/strings/str_split.h"

namespace modle {
struct config {
  std::string_view path_to_bed;
  std::string_view output_dir;
  bool force{false};

  uint32_t bin_size{1'000};
  uint32_t simulation_iterations{5'000'000};
  uint32_t number_of_barriers;
  uint32_t number_of_lefs;
  uint32_t average_lef_processivity{100'000};
  double probability_of_barrier_block{0.95};
  double probability_of_lef_rebind{1.0};
  double probability_of_extrusion_unit_bypass{0.25};
  bool make_heatmaps{true};
  uint64_t seed{0};
  bool skip_burnin{false};
  uint32_t min_n_of_burnin_rounds{0};
  uint32_t min_n_of_loops_per_lef{1};
  double lef_unloader_strength{1};
  uint32_t contact_sampling_interval{1000};
  bool randomize_contact_sampling{false};

  void print() const {
    const std::string padding_placeholder = "{modle-padding}";
    const std::string buf = absl::StrFormat(
        // clang-format off
        "{modle-padding} CONFIG SUMMARY {modle-padding}\n\n"
        "########## Input/Output {modle-padding}\n"
        "##    Input BED                        #  '%s'\n"
        "##    Output directory                 #  '%s'\n"
        "##    Overwrite existing output files  #  %s\n"
        "########## General settings {modle-padding}\n"
        "##    Bin size (bp)                    #  %lu\n"
        "##    # of iterations                  #  %lu\n"
        "##    Avg. LEF processivity (bp)       #  %lu\n"
        "##    LEF unloader strength            #  %.4f\n"
        "##    # of randomly generated barriers #  %lu\n"
        "##    # of randomly generated LEFs     #  %lu\n"
        "##    Skip burn-in                     #  %s\n"
        "##    Contact sampling interval        #  %lu\n"
        "########## Probabilities {modle-padding}\n"
        "##    Prob. of barrier block           #  %.4f\n"
        "##    Prob. of LEF rebind              #  %.4f\n"
        "##    Prob. of LEF bypass              #  %.4f\n"
        "##    Randomize contact sampling       #  %s\n"
        "########## Burn-in {modle-padding}\n"
        "##    Min. # of burn-in rounds         #  %lu\n"
        "##    Min. # of loops per LEF          #  %lu\n"
        "########## Various {modle-padding}\n"
        "##    Generate heatmaps                #  %s\n"
        "##    Seed                             #  %lu\n",
        // clang-format on
        std::filesystem::weakly_canonical(this->path_to_bed),
        std::filesystem::weakly_canonical(this->output_dir), this->force ? "Yes" : "No",
        this->bin_size, this->simulation_iterations, this->average_lef_processivity,
        this->lef_unloader_strength, this->number_of_barriers, this->number_of_lefs,
        this->skip_burnin ? "Yes" : "No", this->contact_sampling_interval,
        this->probability_of_barrier_block, this->probability_of_lef_rebind,
        this->probability_of_extrusion_unit_bypass, this->randomize_contact_sampling ? "Yes" : "No",
        this->min_n_of_burnin_rounds, this->min_n_of_loops_per_lef,
        this->make_heatmaps ? "Yes" : "No", this->seed);

    const auto& toks = absl::StrSplit(buf, '\n');
    uint32_t max_col_width =  // Find longest line (excluding {modle-padding})
        std::max_element(toks.begin(), toks.end(),
                         [&](std::string_view s1, std::string_view s2) {
                           if (s1.find_first_of(padding_placeholder) != std::string::npos)
                             return s1.size() - padding_placeholder.size() < s2.size();
                           if (s2.find_first_of(padding_placeholder) != std::string::npos)
                             return s1.size() < s2.size() - padding_placeholder.size();
                           return s1.size() < s2.size();
                         })
            ->size() +
        3;

    bool first_line = true;
    for (const auto& tok : toks) {
      if (first_line) {
        // Deal with the first line special case.
        // Example:
        // ###### Title ######
        // ###             ###
        std::string title(tok.begin() + tok.find(padding_placeholder) + padding_placeholder.size(),
                          tok.begin() + tok.rfind(padding_placeholder));
        double paddding_length = static_cast<double>(max_col_width - title.size()) / 2 + 1;
        absl::FPrintF(stderr, "%s%s%s\n", std::string(std::floor(paddding_length), '#'), title,
                      std::string(std::ceil(paddding_length), '#'));
        absl::FPrintF(stderr, "### %*s\n", max_col_width - 2, "###");
        first_line = false;
        continue;
      }

      if (tok == "" || tok == "\n") continue;  // Skip empty lines
      if (tok.find(padding_placeholder) ==
          std::string::npos) {  // Display the option and its value with the proper padding.
                                // Example:
                                // ##    Option 1  #        10 ##
        absl::FPrintF(stderr, "%s %*s\n", tok, max_col_width - tok.size() + 1, "##");
      } else {
        // Print Option group title with the appropriate padding.
        // Example:
        // ########## Group 1 ####################
        std::string title(tok.begin(), tok.begin() + tok.find(padding_placeholder));
        std::string rpad(max_col_width - title.size() + 2, '#');
        absl::FPrintF(stderr, "%s%s\n", title, rpad);
      }
    }
    std::string title = "   END OF CONFIG SUMMARY   ";
    double padding_length = static_cast<double>(max_col_width - title.size()) / 2 + 1;
    absl::FPrintF(stderr, "%s%s%s\n\n\n", std::string(std::floor(padding_length), '#'), title,
                  std::string(std::ceil(padding_length), '#'));
  }
};

}  // namespace modle