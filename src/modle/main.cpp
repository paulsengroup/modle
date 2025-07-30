// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#include <fmt/ranges.h>
#include <fmt/std.h>
#include <spdlog/common.h>
#include <spdlog/logger.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

#include <CLI/CLI.hpp>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <exception>
#include <filesystem>
#include <iosfwd>
#include <memory>
#include <new>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "./cli.hpp"
#include "modle/common/chrono.hpp"
#include "modle/common/common.hpp"
#include "modle/common/fmt_helpers.hpp"
#include "modle/common/simulation_config.hpp"
#include "modle/common/string_utils.hpp"
#include "modle/config/version.hpp"
#include "modle/simulation.hpp"

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
static std::atomic<bool> logger_ready{false};

void setup_logger_console(const std::uint8_t verbosity) {
  auto stderr_sink = std::make_shared<spdlog::sinks::stderr_color_sink_mt>();
  //                        [2021-08-12 17:49:34.581] [info]: my log msg
  stderr_sink->set_pattern("[%Y-%m-%d %T.%e] %^[%l]%$: %v");
  stderr_sink->set_level(spdlog::level::level_enum{verbosity});

  auto main_logger = std::make_shared<spdlog::logger>("main_logger", stderr_sink);
  main_logger->set_level(spdlog::level::debug);

  spdlog::set_default_logger(main_logger);
  SPDLOG_INFO("running MoDLE v{}", modle::config::version::str());

  logger_ready = true;
}

void setup_logger_file(const std::filesystem::path& path_to_log_file) {
  SPDLOG_INFO("complete log will be written to file {}", path_to_log_file);

  auto file_sink =
      std::make_shared<spdlog::sinks::basic_file_sink_mt>(path_to_log_file.string(), true);
  //                      [2021-08-12 17:49:34.581] [139797797574208] [info]: my log msg
  file_sink->set_pattern("[%Y-%m-%d %T.%e] [%t] %^[%l]%$: %v");

  spdlog::logger("tmp_logger", file_sink).info("running MoDLE v{}", modle::config::version::str());

  spdlog::default_logger()->sinks().emplace_back(std::move(file_sink));

  logger_ready = true;
}

static std::string concat_args(std::span<char*> args) {
  assert(!args.empty());
  std::string s{args.front()};
  for (std::size_t i = 1; i < args.size(); ++i) {
    s.append(" ");
    s.append(args[i]);
  }

  return s;
}

void write_param_summary_to_log(const modle::Config& c) {
  SPDLOG_INFO("command: {}", concat_args(c.args));
  SPDLOG_INFO("simulation will use up to {} out of {} available CPU cores.", c.nthreads,
              std::thread::hardware_concurrency());
  if (c.stopping_criterion == modle::Config::StoppingCriterion::contact_density) {
    SPDLOG_INFO("using --target-contact-density={:.2f} as stopping criterion.",
                c.target_contact_density);
  } else {
    SPDLOG_INFO("using --target-number-of-epochs={} as stropping criterion.",
                c.target_simulation_epochs);
  }
  SPDLOG_INFO("contact sampling strategy: {}.",
              modle::Cli::contact_sampling_strategy_map.at(c.contact_sampling_strategy));
  SPDLOG_INFO("contact matrix resolution: {}bp", c.bin_size);
}

std::tuple<int, modle::Cli::subcommand, modle::Config> parse_cli_and_setup_logger(int argc,
                                                                                  char** argv) {
  std::unique_ptr<modle::Cli> cli{nullptr};
  try {
    cli = std::make_unique<modle::Cli>(argc, argv);
    auto config = cli->parse_arguments();
    const auto subcmd = cli->get_subcommand();
    setup_logger_console(config.verbosity);

    if (const auto collisions = cli->detect_path_collisions(config); !collisions.empty()) {
      throw std::filesystem::filesystem_error(fmt::format("{}", fmt::join(collisions, "\n")),
                                              std::make_error_code(std::errc::file_exists));
    }

    if (!config.skip_output) {
      assert(!config.path_to_log_file.empty());
      if (const auto& output_dir = config.path_to_output_prefix.parent_path();
          !output_dir.empty()) {
        std::filesystem::create_directories(output_dir.string());
      }
      setup_logger_file(config.path_to_log_file);

      if (!cli->config_file_parsed()) {
        SPDLOG_INFO("writing simulation parameters to config file {}", config.path_to_config_file);
        cli->write_config_file();
      }
    }

    cli->log_warnings();

    return std::make_tuple(0, subcmd, config);
  } catch (const CLI::ParseError& e) {
    assert(cli);
    //  This takes care of formatting and printing error messages (if any)
    return std::make_tuple(cli->exit(e), modle::Cli::subcommand::help, modle::Config{});

  } catch (const std::filesystem::filesystem_error& e) {
    SPDLOG_ERROR("FAILURE! {}", modle::str_strip_suffix(e.what(), ": File exists"));
    return std::make_tuple(1, modle::Cli::subcommand::help, modle::Config{});
  } catch (const spdlog::spdlog_ex& e) {
    fmt::println(stderr,
                 "FAILURE! An error occurred while setting up the main application logger: {}.",
                 e.what());
    return std::make_tuple(1, modle::Cli::subcommand::help, modle::Config{});
  }
}

template <typename... Args>
void try_log_fatal_error(fmt::format_string<Args...> fmt, Args&&... args) {
  if (logger_ready) {
    assert(spdlog::default_logger());
    SPDLOG_ERROR(fmt, std::forward<Args>(args)...);
    spdlog::shutdown();
  } else {
    fmt::print(stderr, fmt, std::forward<Args>(args)...);
  }
}

int main(int argc, char** argv) noexcept {
  try {
    const auto [ec, subcmd, config] = parse_cli_and_setup_logger(argc, argv);
    if (ec != 0 || subcmd == modle::Cli::subcommand::help) {
      return ec;
    }

    assert(spdlog::default_logger());
    write_param_summary_to_log(config);
    const auto t0 = std::chrono::steady_clock::now();
    modle::Simulation sim(config);

    assert(subcmd == modle::Cli::subcommand::simulate);
    sim.run_simulate();

    SPDLOG_INFO("simulation terminated without errors in {}!\n\nBye.",
                modle::format_duration(std::chrono::steady_clock::now() - t0));
  } catch (const std::bad_alloc& e) {
    try_log_fatal_error("FAILURE! Unable to allocate enough memory: {}.", e.what());
    return 1;
  } catch (const std::exception& e) {
    try_log_fatal_error("FAILURE! An error occurred during simulation: {}.", e.what());
    return 1;
  } catch (...) {
    try_log_fatal_error(
        "FAILURE! An error occurred during simulation: Caught an unhandled exception! "
        "If you see this message, please file an issue on GitHub.");
    return 1;
  }
  spdlog::shutdown();
  return 0;
}
