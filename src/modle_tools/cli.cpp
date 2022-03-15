// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "./cli.hpp"  // for Cli

#include <absl/strings/match.h>       // for EndsWith, EndsWithIgnoreCase
#include <absl/strings/str_cat.h>     // for StrAppend, StrCat
#include <absl/strings/str_format.h>  // for StrAppendFormat
#include <absl/strings/str_split.h>   // for SplitIterator, Splitter, StrSplit, operator!=
#include <absl/strings/strip.h>       // for StripPrefix, StripSuffix
#include <absl/types/variant.h>       // for get, monostate
#include <fmt/format.h>               // for format, join, FMT_STRING, make_format_args
#include <fmt/ostream.h>              // for formatbuf<>::int_type
#include <toml++/toml.h>              // for array::operator[], operator<<, parse, print_to...

#include <CLI/CLI.hpp>
#include <algorithm>                         // for max
#include <array>                             // for array
#include <boost/filesystem/file_status.hpp>  // for directory_file, file_status, regular_file
#include <boost/filesystem/operations.hpp>   // for exists, is_directory, is_regular_file, status
#include <boost/filesystem/path.hpp>         // for path, operator<<, path::iterator, operator==
#include <cassert>                           // for assert
#include <exception>                         // for exception
#include <ostream>                           // for streamsize, stringstream, basic_ostream
#include <stdexcept>                         // for runtime_error, invalid_argument, out_of_range
#include <string>                            // for string, allocator, basic_string
#include <string_view>                       // for basic_string_view, string_view, basic_stri...
#include <thread>                            // for hardware_concurrency
#include <tuple>                             // for ignore
#include <utility>                           // for move
#include <vector>                            // for vector

#include "modle/bed/bed.hpp"  // for Parser, bed_dialects, str_to_bed_dialect_m...
#include "modle/common/cli_utils.hpp"
#include "modle/common/common.hpp"  // for usize
#include "modle/common/utils.hpp"   // for str_float_to_str_int, ConstMap::begin, Con...
#include "modle/config/version.hpp"
#include "modle/cooler/cooler.hpp"             // for Cooler, Cooler::READ_ONLY
#include "modle_tools/modle_tools_config.hpp"  // for eval_config, find_barrier_clusters_config

namespace modle::tools {

namespace config = modle::config;

Cli::Cli(int argc, char** argv) : _argc(argc), _argv(argv), _exec_name(*argv) { this->make_cli(); }

void Cli::make_eval_subcommand() {
  auto& sc = *this->_cli
                  .add_subcommand("evaluate",
                                  "Compare MoDLE's output with other contact matrices using "
                                  "various correlation tests.")
                  ->fallthrough()
                  ->preparse_callback([this]([[maybe_unused]] usize i) {
                    assert(this->_config.index() == 0);
                    this->_config = eval_config{};
                  });
  sc.alias("eval");

  this->_config = eval_config{};
  auto& c = absl::get<eval_config>(this->_config);

  auto& io = *sc.add_option_group("Input/Output", "");
  auto& gen = *sc.add_option_group("Generic", "");

  // clang-format off
  io.add_option(
     "-i,--input-matrix",
     c.path_to_input_matrix,
     "Path to a contact matrix in Cooler format.")
     ->check(CLI::ExistingFile)
     ->required();

  io.add_option(
     "-o,--output-prefix",
     c.output_prefix,
     "Output prefix (including directories) to use for output.")
     ->required();

  io.add_option(
     "-r,--reference-matrix",
     c.path_to_reference_matrix,
     "Path to contact matrix in Cooler format to use as reference.\n")
     ->required();

  io.add_option(
     "--chromosome-subrange-file",
     c.path_to_chrom_subranges,
     "Path to BED file with subranges of the chromosomes to be processed.")
     ->check(CLI::ExistingFile);

  io.add_option(
     "--weight-file",
     c.path_to_weights,
     "Path to a TSV file containing the weights to apply to pixels before computing the correlation values.\n"
     "Other TSV files should have an header with column names. Columns \"chrom\" and \"diag\" are required"
     " and should store chromosome names and distance from the diagonal in bins respectively.\n"
     "The name of the column storing the weights can be specified through --weight-column-name.")
     ->check(CLI::ExistingFile);

  io.add_flag(
     "-f,--force",
     c.force,
     "Overwrite existing file(s).")
     ->capture_default_str();

  gen.add_option(
      "--weight-column-name",
      c.weight_column_name,
      "Name of the column from the file specified though --weight-file containing the weight values.")
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
      "Transformation method to apply to the input contact matrix.")
      ->transform(CLI::CheckedTransformer(eval_config::metric_map, CLI::ignore_case))
      ->capture_default_str();

  gen.add_option("-b,--bin-size",
     c.bin_size,
     "Bin size in base pairs.\n"
     "Only used when the contact matrix passed to --reference-matrix is in .mcool format.")
     ->check(CLI::PositiveNumber)
     ->transform(utils::str_float_to_str_int);

  gen.add_option(
     "-w,--diagonal-width",
     c.diagonal_width,
     "Diagonal width to use when computing correlation coefficients.")
     ->check(CLI::PositiveNumber)
     ->transform(utils::str_float_to_str_int)
     ->required();

  gen.add_flag(
      "--reciprocal-weights",
      c.reciprocal_weights,
      "Weights are computed by taking the reciprocal of numbers read from the file specified through --weight-file.")
      ->capture_default_str();

  gen.add_flag(
      "--exclude-zero-pixels,!--include-zero-pixels",
      c.exclude_zero_pxls,
      "When computing the correlation between a pair of rows/columns of pixels, exclude bins corresponding to pixels that are 0 in either of the rows or columns.")
      ->capture_default_str();

  gen.add_option(
     "-t,--threads",
     c.nthreads,
     "CPU threads to allocate.")
     ->check(CLI::PositiveNumber)
     ->transform(CLI::Bound(1U, std::thread::hardware_concurrency()));

  // clang-format on
  gen.get_option("--reciprocal-weights")->needs(io.get_option("--weight-file"));
  gen.get_option("--weight-column-name")->needs(io.get_option("--weight-file"));
  this->_config = absl::monostate{};
}

void Cli::make_find_barrier_clusters_subcommand() {
  auto& sc = *this->_cli
                  .add_subcommand("find-barrier-clusters",
                                  "Detect clusters of extrusion barriers given a BED file.")
                  ->fallthrough()
                  ->preparse_callback([this]([[maybe_unused]] usize i) {
                    assert(this->_config.index() == 0);
                    this->_config = find_barrier_clusters_config{};
                  });

  sc.alias("fbcl");

  this->_config = find_barrier_clusters_config{};
  auto& c = absl::get<find_barrier_clusters_config>(this->_config);

  auto& io = *sc.add_option_group("Input/Output", "");
  auto& cluster = *sc.add_option_group("Cluster properties", "");

  // clang-format off
  io.add_option(
      "-i,--input",
      c.path_to_input_barriers,
      "Path to a BED file with the collection of extrusion barriers to be processed.")
      ->check(CLI::ExistingFile)
      ->required();

  io.add_option(
      "-o,--output-name",
      c.path_to_output,
      "Path to output file. When not provided, barrier clusters will be written to stdout.");

  io.add_option(
      "-b,--breaking-points",
      c.path_to_breaking_points,
      "Path to a BED file listing a the coordinates of cluster breaking points. "
      "These could for instance be a list of genes or enhancers.")
      ->check(CLI::ExistingFile);

  io.add_flag(
      "-f,--force",
      c.force,
      "Overwrite existing file(s).")
      ->capture_default_str();

  io.add_flag(
      "-q,--quiet",
      c.quiet,
      "Don't print warnings to stderr.")
      ->capture_default_str();

  cluster.add_option(
      "-w,--extension-window",
      c.extension_window,
      "Size of the extension window in bp.\n"
      "Clusters of extrusion barriers are identified by finding the first barrier not belonging to any cluster, "
      "then extending the cluster by --extension-window bp, until extending the cluster does not increase the "
      "number of barriers in the cluster.")
      ->check(CLI::PositiveNumber);

  cluster.add_option(
      "--min-cluster-span",
      c.min_cluster_span,
      "The minimum span in bp of a cluster of barriers.\n"
      "--min-cluster-span=0 can be used to allow clusters of any size less than --max-cluster-span.")
      ->check(CLI::NonNegativeNumber)
      ->capture_default_str();

  cluster.add_option(
      "--max-cluster-span",
      c.max_cluster_span,
      "The maximum span in bp of a cluster of barriers.\n"
      "--max-cluster-span=0 can be used to allow clusters to grow indefinitely.")
      ->check(CLI::NonNegativeNumber)
      ->capture_default_str();

  cluster.add_option(
      "--min-cluster-size",
      c.min_cluster_size,
      "The minimum number of barriers that can belong to a given cluster.\n")
      ->check(CLI::NonNegativeNumber)
      ->capture_default_str();

  cluster.add_option(
      "--max-cluster-size",
      c.max_cluster_size,
      "The maximum number of barriers that can belong to a given cluster.\n"
      "--max-cluster-size=0 can be used to allow clusters with an unlimited number of barriers.")
      ->check(CLI::NonNegativeNumber)
      ->capture_default_str();
  // clang-format on
  this->_config = absl::monostate{};
}

void Cli::make_noisify_subcommand() {
  auto& sc =
      *this->_cli
           .add_subcommand("noisify", "Add noise to MoDLE's contact matrix in Cooler format.")
           ->fallthrough()
           ->preparse_callback([this]([[maybe_unused]] usize i) {
             assert(this->_config.index() == 0);
             this->_config = noisify_config{};
           });

  this->_config = noisify_config{};
  auto& c = absl::get<noisify_config>(this->_config);

  auto& io = *sc.add_option_group("Input/Output", "");
  auto& gen = *sc.add_option_group("Generic", "");
  auto& noise = *sc.add_option_group("Noise properties", "");

  // clang-format off
  io.add_option(
     "-i,--input",
     c.path_to_input_matrix,
     "Path to a contact matrix in Cooler format.")
     ->check(CLI::ExistingFile)
     ->required();

  io.add_option(
     "-o,--output-name",
     c.path_to_output_matrix,
     "Output file name (possibly including directories) to use for output.")
     ->required();

  io.add_flag(
     "-f,--force",
     c.force,
     "Overwrite existing file(s).")
     ->capture_default_str();

  gen.add_option(
     "-w,--diagonal-width",
     c.diagonal_width,
     "Diagonal width of the input contact matrix.")
     ->check(CLI::PositiveNumber)
     ->transform(utils::str_float_to_str_int)
     ->required();

  gen.add_option(
      "--bin-size",
      c.bin_size,
      "Bin size of the matrix to noisify. Required and only used when the input matrix is in .mcool format.")
      ->check(CLI::PositiveNumber);

  noise.add_option(
     "--mu,--location",
     c.genextreme_mu,
     "Location parameter (mu) of the generalized extreme value used to.add noise to the contact matrix.")
     ->check(CLI::NonNegativeNumber)
     ->capture_default_str();

  noise.add_option(
     "--sigma,--scale",
     c.genextreme_sigma,
     "Scale parameter (sigma) of the generalized extreme value used to.add noise to the contact matrix.")
     ->check(CLI::NonNegativeNumber)
     ->capture_default_str();

  noise.add_option(
     "--xi,--shape",
     c.genextreme_xi,
     "Shape parameter (xi) of the generalized extreme value used to.add noise to the contact matrix.")
     ->check(CLI::NonNegativeNumber)
     ->capture_default_str();

  noise.add_option(
     "--seed",
     c.seed,
     "Seed used to initialize the PRNG.")
     ->check(CLI::NonNegativeNumber)
     ->capture_default_str();
  // clang-format on
  this->_config = absl::monostate{};
}

void Cli::make_transform_subcommand() {
  auto& sc = *this->_cli
                  .add_subcommand("transform",
                                  "Transform contacts in a matrix in Cooler format for a given "
                                  "resolution and bin size.")
                  ->fallthrough()
                  ->preparse_callback([this]([[maybe_unused]] usize i) {
                    assert(this->_config.index() == 0);
                    this->_config = transform_config{};
                  });

  this->_config = transform_config{};
  auto& c = absl::get<transform_config>(this->_config);

  auto& io = *sc.add_option_group("Input/Output", "");
  auto& trans = *sc.add_option_group("Transformations", "");
  auto& mat = *sc.add_option_group("Contact matrix", "");
  auto& gen = *sc.add_option_group("Generic", "");

  // clang-format off
  io.add_option(
      "-i,--input-matrix",
      c.path_to_input_matrix,
      "Path to a contact matrix in Cooler format.")
      ->check(CLI::ExistingFile)
      ->required();

  io.add_option(
      "-o,--output-matrix",
      c.path_to_output_matrix,
      "Path to output matrix.")
      ->required();

  io.add_flag(
      "-f,--force",
      c.force,
      "Overwrite existing file(s).")
      ->capture_default_str();

  trans.add_option(
      "-m,--method",
      c.method,
      "Transformation method to apply to the input contact matrix.")
      ->transform(CLI::CheckedTransformer(transform_config::transformation_map, CLI::ignore_case))
      ->required();

  trans.add_option(
      "--normalization-range",
      c.normalization_range,
      "Pair of comma separated numbers representing the range of normalized values.")
      ->check(CLI::NonNegativeNumber)
      ->delimiter(',')
      ->capture_default_str();

  trans.add_option(
      "--saturation-range",
      c.saturation_range,
      "Pair of comma separated numbers representing the saturation values. Value saturation is applied prior to normalization.")
      ->delimiter(',')
      ->capture_default_str();

  trans.add_option(
      "--gaussian-blur-sigma",
      c.gaussian_blur_sigma,
      "Sigma parameter of the gaussian kernel used during transformation.")
      ->check(CLI::PositiveNumber)
      ->capture_default_str();

  trans.add_option(
      "--gaussian-blur-multiplier",
      c.gaussian_blur_sigma_multiplier,
      "Multiplier to use when transforming using the difference of gaussians. This multiplier is used to compute the sigma value of the more blurry contact matrix based on the sigma value specifie through --gaussian-blur-sigma.")
      ->check(CLI::PositiveNumber)
      ->capture_default_str();

  trans.add_option(
      "--binary-discretization-value",
      c.discretization_val,
      "Partition value used to discretize the transformed matrix to binary values. Values strictly larger than the threshold will be discretized to one while all other values will be discretized to zero.");

  trans.add_option(
      "--discretization-ranges-tsv",
      c.path_to_discretization_ranges_tsv,
      "Path to a TSV with the discratization ranges.\n"
      "The TSV should have three columns.\n"
      "The first two columns contain the first and last values (open-close inteval) of a discretization interval, while the third column contains the discretization value.\n"
      "Values not falling in any discretization interval are left unmodified.\n"
      "To make an interval start/end at infinity specify -inf or inf.")
      ->check(CLI::ExistingFile);

  trans.add_flag(
      "--float-contacts,!--no-float-contacts",
      c.floating_point,
      "Treat contacts as floating point numbers.")
      ->capture_default_str();

  mat.add_option(
      "-b,--bin-size",
      c.bin_size,
      "Bin size in base pairs.\n"
      "Only used when the input contact matrix is in .mcool format.")
      ->check(CLI::PositiveNumber)
      ->transform(utils::str_float_to_str_int);

  mat.add_option(
      "-w,--diagonal-width",
      c.diagonal_width,
      "Diagonal width to use during transformation.")
      ->check(CLI::PositiveNumber)
      ->transform(utils::str_float_to_str_int)
      ->required();

  gen.add_option(
      "-t,--threads",
      c.nthreads,
      "CPU threads to allocate.")
      ->check(CLI::PositiveNumber)
      ->transform(CLI::Bound(1U, std::thread::hardware_concurrency()));

  // clang-format on

  trans.get_option("--binary-discretization-value")
      ->excludes(trans.get_option("--discretization-ranges-tsv"));

  this->_config = absl::monostate{};
}

void Cli::make_cli() {
  this->_cli.name(this->_exec_name);
  this->_cli.description(
      "MoDLE's helper tool.\n"
      "This tool can be used to perform common pre and post-process operations on MoDLE's "
      "input and output files.");
  this->_cli.set_version_flag("-V,--version",
                              "MoDLE-tools-" + std::string{modle::tools::config::version::str()});
  this->_cli.require_subcommand(1);

  this->make_eval_subcommand();
  this->make_find_barrier_clusters_subcommand();
  this->make_noisify_subcommand();
  this->make_transform_subcommand();
}

void Cli::validate_eval_subcommand() const {
  assert(this->_cli.get_subcommand("eval")->parsed());
  std::vector<std::string> errors;
  const auto& c = absl::get<eval_config>(this->_config);

  try {
    std::ignore =
        cooler::Cooler(c.path_to_input_matrix, cooler::Cooler<>::IO_MODE::READ_ONLY, c.bin_size);
    std::ignore = cooler::Cooler(c.path_to_reference_matrix, cooler::Cooler<>::IO_MODE::READ_ONLY,
                                 c.bin_size);
  } catch (const std::exception& e) {
    if (absl::StartsWith(e.what(), "Cooler::open_file(): bin_size cannot be 0")) {
      errors.emplace_back(
          fmt::format(FMT_STRING("File {} appears to be in .mcool format. --bin-size is required"
                                 "when one or both contact matrices are in .mcool format"),
                      c.path_to_reference_matrix));
    } else {
      errors.emplace_back(fmt::format(FMT_STRING("{}"), e.what()));
    }
  }

  if (c.metric == eval_config::Metric::custom) {
    const auto& io_group = *this->_cli.get_subcommand("eval")->get_option_group("Input/Output");
    const auto& gen_group = *this->_cli.get_subcommand("eval")->get_option_group("Generic");

    constexpr std::array<std::string_view, 4> gen_option_names{
        "--weight-column-name", "--reciprocal-weights", "--exclude-zero-pixels", "--normalize"};

    if (const auto* name = "--weight-file"; !io_group.get_option(name)->empty()) {
      errors.emplace_back(fmt::format(FMT_STRING("{} is not allowed when --metric=custom."), name));
    }

    for (const auto& name : gen_option_names) {
      if (!gen_group.get_option(std::string{name})->empty()) {
        errors.emplace_back(
            fmt::format(FMT_STRING("{} is not allowed when --metric=custom."), name));
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
        const boost::filesystem::path fname = fmt::format(
            FMT_STRING("{}_{}_{}.{}"), c.output_prefix.string(), metric_name, direction, extension);
        auto collision =
            utils::detect_path_collision(fname, c.force, boost::filesystem::regular_file);
        if (!collision.empty()) {
          name_collisions.emplace_back(std::move(collision));
        }
      }
    }
    DISABLE_WARNING_PUSH
    if (!name_collisions.empty()) {
      errors.emplace_back(fmt::format(
          FMT_STRING("Detected {} file name collision(s): refusing to proceed. Pass --force to "
                     "overwrite existing file(s).\n   Colliding file name(s):\n    - {}"),
          name_collisions.size(), fmt::join(name_collisions, ".\n    - ")));
    }
  }

  if (!errors.empty()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("The following error(s) where encountered while validating CLI "
                               "arguments and input file(s):\n - {}"),
                    fmt::join(errors, "\n - ")));
  }
}

void Cli::validate_find_barrier_clusters_subcommand() const {
  assert(this->_cli.get_subcommand("find-barrier-clusters")->parsed());
  std::vector<std::string> errors;
  const auto& c = absl::get<find_barrier_clusters_config>(this->_config);

  assert(boost::filesystem::exists(c.path_to_input_barriers));

  if (auto collision =
          utils::detect_path_collision(c.path_to_output, c.force, boost::filesystem::regular_file);
      !collision.empty()) {
    errors.push_back(collision);
  }

  if (c.min_cluster_span >= c.max_cluster_span && c.max_cluster_span != 0) {
    errors.emplace_back(fmt::format(
        FMT_STRING("--min-cluster-span should be strictly less than --max-cluster-span.\n"
                   "Found:\n"
                   "  --min-cluster-span={}\n"
                   "  --max-cluster-span={}"),
        c.min_cluster_span, c.max_cluster_span));
  }

  if (c.min_cluster_size >= c.max_cluster_size && c.max_cluster_size != 0) {
    errors.emplace_back(fmt::format(
        FMT_STRING("--min-cluster-size should be strictly less than --max-cluster-size.\n"
                   "Found:\n"
                   "  --min-cluster-size={}\n"
                   "  --max-cluster-size={}"),
        c.min_cluster_size, c.max_cluster_size));
  }

  if (!errors.empty()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("The following error(s) where encountered while validating CLI "
                               "arguments and input file(s):\n - {}"),
                    fmt::join(errors, "\n - ")));
  }
}

void Cli::validate_noisify_subcommand() const {
  std::string errors;
  assert(this->_cli.get_subcommand("noisify")->parsed());
  const auto& c = absl::get<noisify_config>(this->_config);

  assert(boost::filesystem::exists(c.path_to_input_matrix));
  try {
    cooler::Cooler f(c.path_to_input_matrix, cooler::Cooler<>::IO_MODE::READ_ONLY, c.bin_size);
  } catch (const std::runtime_error& e) {
    if (absl::EndsWith(e.what(),
                       "A bin size other than 0 is required when calling "
                       "Cooler::validate_multires_cool_flavor()")) {
      absl::StrAppendFormat(
          &errors,
          "File \"%s\" appears to be a multi-resolution Cooler. --bin-size is a "
          "mandatory argument and must be different than 0 when processing .mcool files.\n",
          c.path_to_input_matrix.string());
    } else {
      absl::StrAppendFormat(&errors,
                            "Validation of file \"%s\" failed with the following error: %s.\n",
                            c.path_to_input_matrix.string(), e.what());
    }
  }

  if (!c.force && boost::filesystem::exists(c.path_to_output_matrix)) {
    const auto path_type = boost::filesystem::status(c.path_to_output_matrix).type();
    if (path_type == boost::filesystem::directory_file) {
      absl::StrAppendFormat(&errors,
                            "File \"%s\" already exists and is actually a directory. Please remove "
                            "the directory and try again.",
                            c.path_to_output_matrix.string());
    } else {
      absl::StrAppendFormat(&errors, "File \"%s\" already exists. Pass --force to overwrite",
                            c.path_to_output_matrix.string());
    }
  }
  if (!errors.empty()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("The following error(s) where encountered while validating CLI "
                               "arguments and input file(s):\n{}"),
                    errors));
  }
}

void Cli::validate_transform_subcommand() const {
  assert(this->_cli.get_subcommand("transform")->parsed());
  std::vector<std::string> errors;
  const auto& c = absl::get<transform_config>(this->_config);

  if (auto collision = utils::detect_path_collision(c.path_to_output_matrix, c.force,
                                                    boost::filesystem::regular_file);
      !collision.empty()) {
    errors.push_back(collision);
  }

  try {
    std::ignore =
        cooler::Cooler(c.path_to_input_matrix, cooler::Cooler<>::IO_MODE::READ_ONLY, c.bin_size);
  } catch (const std::exception& e) {
    if (absl::StartsWith(e.what(), "Cooler::open_file(): bin_size cannot be 0")) {
      errors.emplace_back(
          fmt::format(FMT_STRING("File {} appears to be in .mcool format. --bin-size is required"
                                 "when --input-matrix is in .mcool format"),
                      c.path_to_input_matrix));
    } else {
      errors.emplace_back(fmt::format(FMT_STRING("{}"), e.what()));
    }
  }

  using t = transform_config::Transformation;
  if (c.method == t::normalize) {
    if (const auto [lb, ub] = c.normalization_range; lb >= ub) {
      errors.push_back(
          fmt::format(FMT_STRING("The upper bound for the normalization range specified through "
                                 "--normalization range is expected to be strictly larger than the "
                                 "lower bound.\n   - Lower bound: {}\n   - Upper bound: {}\n"),
                      lb, ub));
    }
  }

  if (const auto [lb, ub] = c.normalization_range; !std::isinf(lb) || !std::isinf(ub)) {
    if (lb >= ub) {
      errors.push_back(
          fmt::format(FMT_STRING("The upper bound for the saturation range specified through "
                                 "--saturation range is expected to be strictly larger than the "
                                 "lower bound.\n   - Lower bound: {}\n   - Upper bound: {}\n"),
                      lb, ub));
    }
  }

  if (!errors.empty()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("The following error(s) where encountered while validating CLI "
                               "arguments and input file(s):\n - {}"),
                    fmt::join(errors, "\n - ")));
  }
}

void Cli::validate() const {
  if (this->_cli.get_subcommand("eval")->parsed()) {
    this->validate_eval_subcommand();
  } else if (this->_cli.get_subcommand("find-barrier-clusters")->parsed()) {
    this->validate_find_barrier_clusters_subcommand();
  } else if (this->_cli.get_subcommand("noisify")->parsed()) {
    this->validate_noisify_subcommand();
  } else if (this->_cli.get_subcommand("transform")->parsed()) {
    this->validate_transform_subcommand();
  } else {
    MODLE_UNREACHABLE_CODE;
  }
}

bool Cli::is_ok() const noexcept {
  return (this->_exit_code != 0) && this->_subcommand != subcommand::help;
}
Cli::subcommand Cli::get_subcommand() const noexcept { return this->_subcommand; }
std::string_view Cli::get_printable_subcommand() const noexcept {
  return Cli::subcommand_to_str(this->get_subcommand());
}
modle::tools::modle_tools_config Cli::parse_arguments() {
  this->_cli.name(this->_exec_name);
  this->_cli.parse(this->_argc, this->_argv);

  try {
    if (this->_cli.get_subcommand("evaluate")->parsed()) {
      this->_subcommand = subcommand::eval;
    } else if (this->_cli.get_subcommand("find-barrier-clusters")->parsed()) {
      this->_subcommand = subcommand::fbcl;
    } else if (this->_cli.get_subcommand("noisify")->parsed()) {
      this->_subcommand = subcommand::noisify;
    } else if (this->_cli.get_subcommand("transform")->parsed()) {
      this->_subcommand = subcommand::transform;
    } else {
      this->_subcommand = subcommand::help;
    }
  } catch (const CLI::ParseError& e) {
    //  This takes care of formatting and printing error messages (if any)
    this->_exit_code = this->_cli.exit(e);
    return this->_config;
  } catch (const std::exception& e) {
    this->_exit_code = 1;
    throw std::runtime_error(fmt::format(
        FMT_STRING(
            "An unexpected error has occurred while parsing CLI arguments: {}. If you see this "
            "message, please file an issue on GitHub"),
        e.what()));

  } catch (...) {
    this->_exit_code = 1;
    throw std::runtime_error(
        "An unknown error occurred while parsing CLI arguments! If you see this message, please "
        "file an issue on GitHub");
  }
  this->validate();
  this->_exit_code = 0;
  return this->_config;
}

int Cli::exit(const CLI::ParseError& e) const { return this->_cli.exit(e); }

std::string Cli::to_json() const {
  std::string buff;
  for (const auto& line : absl::StrSplit(this->_cli.config_to_str(true, false), '\n')) {
    if (line.empty()) {
      continue;
    }
    if (line.front() == '[') {  // Begin of the section for the active subcommand
      absl::StrAppend(&buff, line, "\n");
      continue;
    }
    // Given two subcommands named comm1 and comm2, assuming comm1 was parsed while comm2 was not,
    // the TOML produced by CLI11 will have values for comm2 formatted as comm2.myarg1=1,
    // comm2.myarg2="a" etc.
    // All we are doing here is to look for an argument name containing '.'.
    // In this way we can filter out entry corresponding to arguments for inactive subcommands
    const auto arg = line.substr(0, line.find('='));
    assert(!arg.empty());
    if (arg.find('.') == decltype(arg)::npos) {
      absl::StrAppend(&buff, line, "\n");
    }
  }

  try {
    auto tt = toml::parse(buff);
    std::stringstream ss;
    ss << toml::json_formatter{tt};
    return ss.str();
  } catch (const std::exception& e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING(
            "The following error occurred while converting MoDLE's config from TOML to JSON: {}"),
        e.what()));
  }
}

std::string_view Cli::subcommand_to_str(subcommand s) noexcept {
  switch (s) {
    case eval:
      return "evaluate";
    case fbcl:
      return "find-barrier-clusters";
    case noisify:
      return "noisify";
    case transform:
      return "transform";
    default:
      assert(s == help);
      return "--help";
  }
}

}  // namespace modle::tools
