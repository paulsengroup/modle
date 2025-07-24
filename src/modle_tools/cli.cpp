// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "./cli.hpp"

#include <absl/strings/match.h>
#include <absl/strings/str_cat.h>
#include <absl/strings/str_format.h>
#include <absl/strings/str_split.h>
#include <absl/strings/strip.h>
#include <fmt/format.h>
#include <fmt/std.h>
#include <toml++/toml.h>

#include <CLI/CLI.hpp>
#include <algorithm>
#include <array>
#include <cassert>
#include <exception>
#include <filesystem>
#include <hictk/cooler.hpp>
#include <ostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <tuple>
#include <utility>
#include <variant>
#include <vector>

#include "modle/bed/bed.hpp"
#include "modle/bigwig/bigwig.hpp"
#include "modle/common/cli_utils.hpp"
#include "modle/common/common.hpp"
#include "modle/common/fmt_helpers.hpp"
#include "modle/common/utils.hpp"
#include "modle/config/version.hpp"
#include "modle_tools/modle_tools_config.hpp"

namespace modle::tools {

class CoolerFileValidator : public CLI::Validator {
 public:
  CoolerFileValidator() : Validator("Cooler") {
    func_ = [](std::string& uri) -> std::string {
      if (!hictk::cooler::utils::is_cooler(uri)) {
        if (hictk::cooler::utils::is_multires_file(uri)) {
          return "URI points to a .mcool file: " + uri;
        }
        return "Not a valid Cooler: " + uri;
      }
      return "";
    };
  }
};

class BigWigFileValidator : public CLI::Validator {
 public:
  BigWigFileValidator() : Validator("BigWig") {
    func_ = [](std::string& path) -> std::string {
      try {
        io::bigwig::Reader bw{std::filesystem::path(path)};
        if (bw.chromosomes().empty()) {
          throw std::runtime_error("");
        }
      } catch (const std::exception& e) {
        return "Not a valid bigWig file: " + path;
      }
      return "";
    };
  }
};

inline const auto IsValidCoolerFile = CoolerFileValidator();
inline const auto IsValidBigWigFile = BigWigFileValidator();

namespace config = modle::config;

// These stream operators are needed to properly serialize enums when writing the config file to
// TOML or JSON
std::ostream& operator<<(std::ostream& os, const eval_config::Metric& metric) {
  os << eval_config::metric_map.at(metric);
  return os;
}

std::ostream& operator<<(std::ostream& os, const transform_config::Transformation& method) {
  os << transform_config::transformation_map.at(method);
  return os;
}

Cli::Cli(int argc, char** argv) : _argc(argc), _argv(argv), _exec_name(*argv) { make_cli(); }

void Cli::make_annotate_barriers_subcommand() {
  auto& sc =
      *_cli.add_subcommand("annotate-barriers",
                           "Helper tool to generate an extrusion barrier annotation for MoDLE.")
           ->fallthrough()
           ->preparse_callback([this]([[maybe_unused]] usize i) {
             assert(_config.index() == 0);
             _config = annotate_barriers_config{};
           });

  _config = annotate_barriers_config{};
  auto& c = std::get<annotate_barriers_config>(_config);

  auto& io = *sc.add_option_group("IO", "");
  auto& gen = *sc.add_option_group("Generic", "");

  // clang-format off
  io.add_option(
     "chip-bigwig",
     c.path_to_bigwig,
     "Path to a BigWig file.\n"
     "This is usually a ChIP-seq BigWig file with the fold-change over control for CTCF/RAD21.\n"
     "Values from this file are used to compute extrusion barrier occupancies.")
     ->check(IsValidBigWigFile)
     ->required();

  io.add_option(
     "binding-sites",
     c.path_to_bed,
     "Path to a BED6+ file with the list of candidate extrusion barriers.\n"
     "This is usually the ChIP-seq broad/narrowPeak for CTCF/RAD21.")
     ->check(CLI::ExistingFile)
     ->required();

  gen.add_option(
     "--scaling-factor",
     c.scaling_factor,
     "Scaling factor to use in the logistic transformation.\n"
     "This value is used to divide peak summits read from the bigWig file provided as input.\n"
     "Refer to section \"Converting CTCF ChIP-Seq to Orientation-Specific BE Permeability\" from\n"
     "\"https://www.cell.com/cell-reports/pdf/S2211-1247(16)30530-7.pdf\" for more details.")
     ->check(CLI::PositiveNumber)
     ->capture_default_str();

  gen.add_option(
     "--occupancy-lower-bound,--lb",
     c.occupancy_lb,
     "Lower bound for the barrier occupancy.\n"
     "Bins with an occupancy smaller than --occupancy-lower-bound will be discarded or clamped,\n"
     "depending on whether --clamp-occupancy was specified.")
     ->check(CLI::Bound(0.0, 1.0))
     ->capture_default_str();

  gen.add_option(
     "--occupancy-upper-bound,--ub",
     c.occupancy_ub,
     "Upper bound for the barrier occupancy.\n"
     "Bins with an occupancy greater than --occupancy-upper-bound will be discarded or clamped,\n"
     "depending on whether --clamp-occupancy was specified.")
     ->check(CLI::Bound(0.0, 1.0))
     ->capture_default_str();

  gen.add_flag(
     "--clamp-occupancy,!--no-clamp-occupancy",
     c.clamp_occupancy,
     "Clamp occupancy values falling outside of the interval specified through --lb and --ub.")
     ->capture_default_str();
  // clang-format on

  _config = std::monostate{};
}

void Cli::make_eval_subcommand() {
  auto& sc =
      *_cli.add_subcommand("evaluate", "Helper tool to compare contact matrices in cooler format.")
           ->fallthrough()
           ->preparse_callback([this]([[maybe_unused]] usize i) {
             assert(_config.index() == 0);
             _config = eval_config{};
           });
  sc.alias("eval");

  _config = eval_config{};
  auto& c = std::get<eval_config>(_config);

  auto& io = *sc.add_option_group("IO", "");
  auto& gen = *sc.add_option_group("Generic", "");

  // clang-format off
  io.add_option(
     "-i,--input-cooler",
     c.input_cooler_uri,
     "Path or URI to a Cooler file to use as input.")
     ->check(IsValidCoolerFile)
     ->required();

  io.add_option(
     "-o,--output-prefix",
     c.output_prefix,
     "Output prefix.\n"
     "Can be an absolute or relative path including the file name but without the extension.")
     ->required();

  io.add_option(
     "-r,--reference-cooler",
     c.reference_cooler_uri,
     "Path or URI to a Cooler file to use as reference.")
     ->check(IsValidCoolerFile)
     ->required();

  io.add_option(
     "--regions-of-interest",
     c.path_to_regions_of_interest_bed,
     "Path to a BED file with the list of regions of interest to be processed.")
     ->check(CLI::ExistingFile);

  io.add_option(
     "-c,--chrom-sizes",
     c.path_to_chrom_sizes,
     "Path to a .chrom.sizes file.\n"
     "This option is required when the Coolers passed with options --input-cooler and --reference-cooler"
     "do not have the same number of chromosomes.")
     ->check(CLI::ExistingFile);

  io.add_option(
     "--weight-file",
     c.path_to_weights,
     "Path to a TSV file containing the weights to apply to pixels before computing the\n"
     "correlation. The TSV should have a header with at least the following columns:\n"
     " - chrom\n"
     " - diag\n"
     "Storing chromosome name and distance from the diagonal in bins to which the weight applies\n"
     "respectively. This file is often generated using cooltools expected-cis.\n"
     "The name of the column storing the weights can be specified through --weight-column-name.")
     ->check(CLI::ExistingFile);

  io.add_flag(
     "-f,--force",
     c.force,
     "Overwrite existing files (if any).")
     ->capture_default_str();

  gen.add_option(
      "--weight-column-name",
      c.weight_column_name,
      "Label for the weight column (see --weight-file description for more details).")
      ->check(CLI::Number)
      ->capture_default_str();

  gen.add_flag(
     "--normalize,!--no-normalize",
     c.normalize,
     "Normalize matrices to range 0-1 before comparing them.")
     ->capture_default_str();

  gen.add_option(
      "-m,--metric",
      c.metric,
      fmt::format("Comparison metric.\n"
                  "Supported metrics:\n"
                  " - {}",
                  fmt::join(eval_config::metric_map.keys(), "\n - ")))
      ->transform(CLI::CheckedTransformer(eval_config::metric_map, CLI::ignore_case))
      ->capture_default_str();

  gen.add_option(
     "-w,--diagonal-width",
     c.diagonal_width,
     "Width of the subdiagonal window to consider for comparison.")
     ->check(CLI::PositiveNumber)
     ->transform(utils::cli::TrimTrailingZerosFromDecimalDigit | utils::cli::AsGenomicDistance)
     ->required();

  gen.add_flag(
      "--reciprocal-weights",
      c.reciprocal_weights,
      "Use the weights reciprocal instead of the actual weights.")
      ->capture_default_str();

  gen.add_flag(
      "--exclude-zero-pixels,!--include-zero-pixels",
      c.exclude_zero_pxls,
      "When comparing rows or columns of pixels, ignore pairs of pixels where one or\n"
      "both pixels are zero.")
      ->capture_default_str();

  gen.add_option(
     "-t,--threads",
     c.nthreads,
     "Number of worker threads.")
     ->check(CLI::PositiveNumber)
     ->transform(CLI::Bound(1U, std::thread::hardware_concurrency()));

  // clang-format on
  gen.get_option("--reciprocal-weights")->needs(io.get_option("--weight-file"));
  gen.get_option("--weight-column-name")->needs(io.get_option("--weight-file"));
  _config = std::monostate{};
}

void Cli::make_transform_subcommand() {
  auto& sc = *_cli.add_subcommand("transform",
                                  "Transform contact matrices in cooler format using one of the "
                                  "supported transformation methods.")
                  ->fallthrough()
                  ->preparse_callback([this]([[maybe_unused]] usize i) {
                    assert(_config.index() == 0);
                    _config = transform_config{};
                  });

  _config = transform_config{};
  auto& c = std::get<transform_config>(_config);

  auto& io = *sc.add_option_group("Input/Output", "");
  auto& trans = *sc.add_option_group("Transformations", "");
  auto& mat = *sc.add_option_group("Contact matrix", "");
  auto& gen = *sc.add_option_group("Generic", "");

  // clang-format off
  io.add_option(
      "-i,--input-cooler",
      c.input_cooler_uri,
      "Path or URI to a Cooler file to use as input.")
      ->check(IsValidCoolerFile)
      ->required();

  io.add_option(
      "-o,--output-cooler",
      c.output_cooler_uri,
      "Path or URI to to the output Cooler.")
      ->required();

  io.add_flag(
      "-f,--force",
      c.force,
      "Overwrite existing files (if any).")
      ->capture_default_str();

  trans.add_option(
      "-m,--method",
      c.method,
      fmt::format("Transformation method to apply to the input contact matrix.\n"
                  "Supported methods:\n"
                  " - {}",
      fmt::join(transform_config::transformation_map.keys(), "\n - ")))
      ->transform(CLI::CheckedTransformer(transform_config::transformation_map, CLI::ignore_case))
      ->required();

  trans.add_option(
      "--normalization-range",
      c.normalization_range,
      "Lower and upper bound for the normalization range as a comma-separate pair of values.")
      ->check(CLI::NonNegativeNumber)
      ->delimiter(',')
      ->capture_default_str();

  trans.add_option(
      "--saturation-range",
      c.saturation_range,
      "Lower and upper bound for the saturation range as a comma-separate pair of values.")
      ->delimiter(',')
      ->capture_default_str();

  trans.add_option(
      "--gaussian-blur-sigma",
      c.gaussian_blur_sigma,
      "Gaussian kernel sigma.\n"
      "Controls the level of blurriness applied to contact matrices when using the\n"
      "difference of Gaussians as transformation method.")
      ->check(CLI::PositiveNumber)
      ->capture_default_str();

  trans.add_option(
      "--gaussian-blur-multiplier",
      c.gaussian_blur_sigma_multiplier,
      "Gaussian blur multiplier.\n"
      "Multiplier used to compute the sigma of the Gaussian kernel used to generate the more\n"
      "blurry contact matrix.")
      ->check(CLI::Bound(1.0, std::numeric_limits<double>::infinity()))
      ->capture_default_str();

  trans.add_option(
      "--binary-discretization-value",
      c.discretization_val,
      "Threshold for the step function used to discretize pixels.\n"
      "Pixels above the threshold are mapped to 1, while all the others are mapped to 0.")
      ->check(utils::cli::IsFinite);

  trans.add_option(
      "--discretization-ranges-tsv",
      c.path_to_discretization_ranges_tsv,
      "Path to a TSV with the discretization ranges.\n"
      "The TSV should have three columns.\n"
      "The first two columns contain the first and last values (open-close interval)\n"
      "of a discretization interval, while the third column contains the discretization value.\n"
      "Values not falling in any discretization interval are left unmodified.\n"
      "To make an interval start/end at infinity specify -inf or inf.")
      ->check(CLI::ExistingFile);

  trans.add_flag(
      "--float-contacts,!--integer-contacts",
      c.floating_point,
      "Use floating point or integer numbers to represent contacts.")
      ->capture_default_str();

  mat.add_option(
      "-w,--diagonal-width",
      c.diagonal_width,
      "Width of the subdiagonal window to transform.\n"
      "Pixels outside the window are set to 0.")
      ->check(CLI::PositiveNumber)
      ->transform(utils::cli::TrimTrailingZerosFromDecimalDigit | utils::cli::AsGenomicDistance)
      ->required();

  gen.add_option(
      "-t,--threads",
      c.nthreads,
      "Number of worker threads.")
      ->check(CLI::PositiveNumber)
      ->transform(CLI::Bound(1U, std::thread::hardware_concurrency()));

  // clang-format on

  trans.get_option("--binary-discretization-value")
      ->excludes(trans.get_option("--discretization-ranges-tsv"));

  _config = std::monostate{};
}

void Cli::make_cli() {
  _cli.name(_exec_name);
  _cli.description(
      "MoDLE's helper tool.\n"
      "This tool can be used to perform common pre and post-process operations\n"
      "on MoDLE's input/output files.");
  _cli.set_version_flag("-V,--version",
                        "MoDLE-tools-" + std::string{modle::tools::config::version::str()});
  _cli.require_subcommand(1);
  _cli.formatter(std::make_shared<utils::cli::Formatter>());

  make_annotate_barriers_subcommand();
  make_eval_subcommand();
  make_transform_subcommand();
}

void Cli::validate_annotate_barriers_subcomand() const {
  assert(_cli.get_subcommand("annotate-barriers")->parsed());
  std::vector<std::string> errors;
  const auto& c = std::get<annotate_barriers_config>(_config);

  if (c.occupancy_lb >= c.occupancy_ub) {
    errors.emplace_back(
        fmt::format("occupancy lower bound should be smaller than the upper bound.\n"
                    "   - Lower bound: {}\n"
                    "   - Upper bound: {}\n",
                    c.occupancy_lb, c.occupancy_ub));
  }

  if (!errors.empty()) {
    throw std::runtime_error(
        fmt::format("the following error(s) where encountered while validating CLI "
                    "arguments and input file(s):\n - {}",
                    fmt::join(errors, "\n - ")));
  }
}

void Cli::validate_eval_subcommand() const {
  assert(_cli.get_subcommand("eval")->parsed());
  std::vector<std::string> errors;
  const auto& c = std::get<eval_config>(_config);

  try {
    const hictk::cooler::File f1(c.input_cooler_uri.string());
    const hictk::cooler::File f2(c.reference_cooler_uri.string());
    if (f1.resolution() != f2.resolution()) {
      errors.emplace_back(fmt::format(
          "Coolers at URIs {} and {} have different resolutions ({} and {} respectively)",
          c.input_cooler_uri, c.reference_cooler_uri, f1.resolution(), f2.resolution()));
    }

    const auto& io_group = *_cli.get_subcommand("eval")->get_option_group("IO");
    if (constexpr auto* name = "--chrom-sizes";
        f1.chromosomes() != f2.chromosomes() && !io_group.get_option(name)->empty()) {
      errors.emplace_back(fmt::format(
          "{} is required when input and reference Coolers do not have the same chromosomes",
          name));
    }

  } catch (const std::exception& e) {
    errors.emplace_back(e.what());
  }

  if (c.metric == eval_config::Metric::custom) {
    const auto& io_group = *_cli.get_subcommand("eval")->get_option_group("IO");
    const auto& gen_group = *_cli.get_subcommand("eval")->get_option_group("Generic");

    constexpr std::array<std::string_view, 4> gen_option_names{
        "--weight-column-name", "--reciprocal-weights", "--exclude-zero-pixels", "--normalize"};

    if (constexpr auto* name = "--weight-file"; !io_group.get_option(name)->empty()) {
      errors.emplace_back(fmt::format("{} is not allowed when --metric=custom.", name));
    }

    for (const auto& name : gen_option_names) {
      if (!gen_group.get_option(std::string{name})->empty()) {
        errors.emplace_back(fmt::format("{} is not allowed when --metric=custom.", name));
      }
    }
  }

  if (!c.force) {
    constexpr std::array<std::string_view, 2> stripe_directions{"horizontal", "vertical"};
    constexpr std::array<std::string_view, 2> extensions{"bw", "tsv.gz"};

    std::vector<std::string> name_collisions;
    const auto metric_name = [&]() {
      const auto* it = std::find_if(eval_config::metric_map.begin(), eval_config::metric_map.end(),
                                    [&](const auto p) { return p.second == c.metric; });
      assert(it != eval_config::metric_map.end());
      return it->first;
    }();

    DISABLE_WARNING_PUSH
    DISABLE_WARNING_RANGE_LOOP_ANALYSIS
    for (const auto direction : stripe_directions) {
      for (const auto extension : extensions) {
        const std::filesystem::path fname =
            fmt::format("{}_{}_{}.{}", c.output_prefix.string(), metric_name, direction, extension);
        auto collision =
            utils::detect_path_collision(fname, c.force, std::filesystem::file_type::regular);
        if (!collision.empty()) {
          name_collisions.emplace_back(std::move(collision));
        }
      }
    }
    DISABLE_WARNING_PUSH
    if (!name_collisions.empty()) {
      errors.emplace_back(
          fmt::format("detected {} file name collision(s): refusing to proceed. Pass --force to "
                      "overwrite existing file(s).\n   Colliding file name(s):\n    - {}",
                      name_collisions.size(), fmt::join(name_collisions, ".\n    - ")));
    }
  }

  if (!errors.empty()) {
    throw std::runtime_error(
        fmt::format("the following error(s) where encountered while validating CLI "
                    "arguments and input file(s):\n - {}",
                    fmt::join(errors, "\n - ")));
  }
}

void Cli::validate_transform_subcommand() const {
  assert(_cli.get_subcommand("transform")->parsed());
  std::vector<std::string> errors;
  const auto& c = std::get<transform_config>(_config);

  if (auto collision = utils::detect_path_collision(c.output_cooler_uri, c.force,
                                                    std::filesystem::file_type::regular);
      !collision.empty()) {
    errors.push_back(collision);
  }

  using t = transform_config::Transformation;
  if (c.method == t::normalize) {
    if (const auto [lb, ub] = c.normalization_range; lb >= ub) {
      errors.push_back(
          fmt::format("the upper bound for the normalization range specified through "
                      "--normalization range is expected to be strictly larger than the "
                      "lower bound.\n   - Lower bound: {}\n   - Upper bound: {}\n",
                      lb, ub));
    }
  }

  if (const auto [lb, ub] = c.normalization_range; !std::isinf(lb) || !std::isinf(ub)) {
    if (lb >= ub) {
      errors.push_back(
          fmt::format("the upper bound for the saturation range specified through "
                      "--saturation range is expected to be strictly larger than the "
                      "lower bound.\n   - Lower bound: {}\n   - Upper bound: {}\n",
                      lb, ub));
    }
  }

  if (!errors.empty()) {
    throw std::runtime_error(
        fmt::format("The following error(s) where encountered while validating CLI "
                    "arguments and input file(s):\n - {}",
                    fmt::join(errors, "\n - ")));
  }
}

void Cli::validate() const {
  if (_cli.get_subcommand("annotate-barriers")->parsed()) {
    validate_annotate_barriers_subcomand();
  } else if (_cli.get_subcommand("eval")->parsed()) {
    validate_eval_subcommand();
  } else if (_cli.get_subcommand("transform")->parsed()) {
    validate_transform_subcommand();
  } else {
    MODLE_UNREACHABLE_CODE;
  }
}

Cli::subcommand Cli::get_subcommand() const noexcept { return _subcommand; }
std::string_view Cli::get_printable_subcommand() const noexcept {
  return Cli::subcommand_to_str(get_subcommand());
}
modle::tools::modle_tools_config Cli::parse_arguments() {
  _cli.name(_exec_name);
  _cli.parse(_argc, _argv);

  try {
    if (_cli.get_subcommand("annotate-barriers")->parsed()) {
      _subcommand = subcommand::annotate_barriers;
    } else if (_cli.get_subcommand("evaluate")->parsed()) {
      _subcommand = subcommand::eval;
    } else if (_cli.get_subcommand("transform")->parsed()) {
      _subcommand = subcommand::transform;
    } else {
      _subcommand = subcommand::help;
    }
  } catch (const CLI::ParseError& e) {
    //  This takes care of formatting and printing error messages (if any)
    _exit_code = _cli.exit(e);
    return _config;
  } catch (const std::exception& e) {
    _exit_code = 1;
    throw std::runtime_error(fmt::format(
        "An unexpected error has occurred while parsing CLI arguments: {}. If you see this "
        "message, please file an issue on GitHub",
        e.what()));

  } catch (...) {
    _exit_code = 1;
    throw std::runtime_error(
        "An unknown error occurred while parsing CLI arguments! If you see this message, please "
        "file an issue on GitHub");
  }
  validate();

  std::visit(
      [&, this](auto& config) {
        if constexpr (!std::is_same_v<absl::remove_cvref_t<decltype(config)>, std::monostate>) {
          config.args = absl::MakeSpan(_argv, static_cast<usize>(_argc));
          config.args_json = to_json();
        }
      },
      _config);

  _exit_code = 0;
  return _config;
}

int Cli::exit(const CLI::ParseError& e) const { return _cli.exit(e); }

std::string Cli::to_json() const {
  const auto prefix = std::string{get_printable_subcommand()} + ".";
  std::string buff;
  for (const auto& line : absl::StrSplit(_cli.config_to_str(true, false), '\n')) {
    if (line.empty() || line.front() == '[') {
      continue;
    }

    if (buff.empty()) {
      buff = fmt::format("[{}]\n", get_printable_subcommand());
    }

    if (line.find(prefix) != std::string::npos) {
      assert(line.find('=') != std::string::npos);
      absl::StrAppend(&buff, line.substr(prefix.size()), "\n");
    }
  }

  assert(!buff.empty());
  try {
    auto tt = toml::parse(buff);
    std::stringstream ss;
    ss << toml::json_formatter{tt};
    return ss.str();
  } catch (const std::exception& e) {
    throw std::runtime_error(fmt::format(
        "The following error occurred while converting MoDLE's config from TOML to JSON: {}",
        e.what()));
  }
}

std::string_view Cli::subcommand_to_str(subcommand s) noexcept {
  switch (s) {
    case annotate_barriers:
      return "annotate-barriers";
    case eval:
      return "evaluate";
    case transform:
      return "transform";
    default:
      assert(s == help);
      return "--help";
  }
}

}  // namespace modle::tools
