#include <absl/container/btree_set.h>  // for btree_set

#include <CLI/CLI.hpp>
#include <cstdint>  // for uint_fast8_t
#include <string>   // for string

#include "modle_tools/config.hpp"  // for config

namespace modle::tools {

class Cli {
 public:
  enum subcommand : uint_fast8_t {
    help,
    eval,
    filter_barriers,
    find_barrier_clusters,
    noisify,
    stats,
  };
  Cli(int argc, char** argv);
  [[nodiscard]] bool is_ok() const;
  [[nodiscard]] subcommand get_subcommand() const;
  [[nodiscard]] config parse_arguments();
  [[nodiscard]] int exit(const CLI::ParseError& e) const;

 private:
  int _argc;
  char** _argv;
  std::string _exec_name;
  int _exit_code{1};
  config _config{};
  CLI::App _cli{};
  subcommand _subcommand{subcommand::help};

  absl::btree_set<std::string> _filtering_criteria{"intersection", "pairwise-intersection"};

  void make_eval_subcommand();
  void make_filter_barriers_subcommand();
  void make_find_barrier_clusters_subcommand();
  void make_noisify_subcommand();
  void make_stats_subcommand();
  void make_cli();

  [[nodiscard]] std::string validate_eval_subcommand();
  [[nodiscard]] std::string validate_filter_barriers_subcommand() const;
  [[nodiscard]] std::string validate_find_barrier_clusters_subcommand() const;
  [[nodiscard]] std::string validate_noisify_subcommand() const;
  [[nodiscard]] std::string validate_stats_subcommand() const;
  [[nodiscard]] bool validate();
};

}  // namespace modle::tools
