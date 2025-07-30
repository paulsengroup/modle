// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#include <fmt/std.h>
#include <spdlog/common.h>
#include <spdlog/logger.h>
#include <spdlog/sinks/sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

#include <CLI/CLI.hpp>
#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstring>
#include <exception>
#include <memory>
#include <new>
#include <stdexcept>
#include <string_view>
#include <variant>
#include <vector>

#include "./cli.hpp"
#include "modle_tools/tools.hpp"

int main(int argc, char** argv) {
  using namespace modle::tools;
  std::unique_ptr<Cli> cli{nullptr};
  spdlog::set_default_logger(std::make_shared<spdlog::logger>("main_logger"));
  {
    auto stderr_sink = spdlog::default_logger()->sinks().emplace_back(
        std::make_shared<spdlog::sinks::stderr_color_sink_mt>());
    //                        [2021-08-12 17:49:34.581] [info]: my log msg
    stderr_sink->set_pattern("[%Y-%m-%d %T.%e] %^[%l]%$: %v");
  }
  try {
    cli = std::make_unique<Cli>(argc, argv);
    const auto config = cli->parse_arguments();

    {
      using subcmd = Cli::subcommand;
      switch (cli->get_subcommand()) {
        case subcmd::annotate_barriers:
          annotate_barriers_subcmd(std::get<annotate_barriers_config>(config));
          return 0;
        case subcmd::eval:
          eval_subcmd(std::get<eval_config>(config));
          return 0;
        case subcmd::transform:
          transform_subcmd(std::get<transform_config>(config));
          return 0;
        default:
          throw std::runtime_error(
              "Default branch in switch statement in modle_tools::main() should be unreachable! If "
              "you see this message, please file an issue on GitHub");
      }
    }
  } catch (const CLI::ParseError& e) {
    assert(cli);
    return cli->exit(e);  //  This takes care of formatting and printing error messages (if any)
  } catch (const std::bad_alloc& err) {
    SPDLOG_ERROR("FAILURE! Unable to allocate enough memory: {}", err.what());
    return 1;
  } catch (const spdlog::spdlog_ex& e) {
    fmt::println(stderr,
                 "FAILURE! An error occurred while setting up the main application logger: {}.",
                 e.what());
    return 1;
  } catch (const std::exception& e) {
    assert(cli);
    SPDLOG_ERROR("FAILURE! modle_tools {} encountered the following error: {}.",
                 cli->get_printable_subcommand(), e.what());
    return 1;
  } catch (...) {
    SPDLOG_ERROR(
        "FAILURE! modle_tools {} encountered the following error: Caught an "
        "unhandled exception! "
        "If you see this message, please file an issue on GitHub.",
        cli->get_printable_subcommand());
    throw;
    return 1;
  }
  return 0;
}
