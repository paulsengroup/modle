#include "convert.hpp"

void convert(modle::utils::config& c) {
  std::string argv;
  modle::tools::convert::check_convert_preconditions(c, argv);
  std::filesystem::create_directories(c.out_dir);
  if (c.convert_to_hic) modle::tools::convert::convert_to_hic(c, argv);
  if (c.convert_to_tsv) modle::tools::convert::convert_to_tsv(c);
}

int main(int argc, char** argv) {
  modle::utils::Cli cli(argc, argv);
  auto c = cli.parse_arguments();
  if (!cli.is_ok()) return c.exit_code;

  std::filesystem::create_directories(c.tmp_dir);
  try {
    switch (cli.get_subcommand()) {
      case modle::utils::Cli::subcommand::convert:
        convert(c);
        break;
      default:
        throw std::runtime_error(
            "Default branch in switch statement in main.cpp should be unreachable! If you see "
            "this message, please open an issue on GitHub");
    }
  } catch (const std::runtime_error& err) {
    absl::FPrintF(stderr, "FAILURE: %s.\n", err.what());
    if (!c.keep_tmp_files && std::filesystem::is_empty(c.tmp_dir))
      std::filesystem::remove_all(c.tmp_dir);
    return 1;
  }
  if (!c.keep_tmp_files && std::filesystem::is_empty(c.tmp_dir))
    std::filesystem::remove_all(c.tmp_dir);
  return 0;
}
