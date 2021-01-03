#pragma once

#include <absl/strings/str_cat.h>
#include <absl/time/clock.h>
#include <fmt/compile.h>
#include <fmt/format.h>
#include <fmt/ostream.h>

#include <atomic>
#include <boost/process/child.hpp>
#include <boost/process/io.hpp>
#include <boost/process/pipe.hpp>
#include <random>
#include <string_view>
#include <thread>

#include "modle/contacts.hpp"

namespace modle::dcic4dn {

std::string short_record::to_string() const {
  return absl::StrCat(this->readID, " ", this->chr1, " ", this->position1, " ", this->chr2, " ",
                      this->position2, " ", this->strand1, " ", this->strand2);
}

void converter::from_modle(const std::vector<std::string>& path_to_cmatrices,
                           std::string_view path_to_output_file, uint64_t noise_range,
                           bool force_write, double noise_stddev, uint8_t compression_lvl,
                           std::size_t nthreads) {
  // Initialize random number generators
  std::unique_ptr<std::normal_distribution<double>> noise_generator = nullptr;
  std::unique_ptr<std::mt19937_64> rand_eng = nullptr;
  if (noise_range != 0) {
    auto init_noise_generator = [&]() {
      constexpr double stdev_99_9_ci = 6.0;  // 99.9% CI
      return std::normal_distribution<double>(
          0, noise_stddev == 0 ? static_cast<double>(noise_range) / stdev_99_9_ci : noise_stddev);
    };
    noise_generator = std::make_unique<std::normal_distribution<double>>(init_noise_generator());
    std::random_device rand_dev;
    std::seed_seq seeder{rand_dev()};
    rand_eng = std::make_unique<std::mt19937_64>(seeder);
  }

  // Useful lambdas
  auto process_pos = [&](const double pos, const double chr_start, const double chr_end) {
    if (noise_generator != nullptr) {
      assert(rand_eng);  // NOLINT
      return static_cast<uint64_t>(
          std::clamp(std::round(pos + (*noise_generator)(*rand_eng)), chr_start, chr_end));
    }
    return static_cast<uint64_t>(pos);
  };

  auto argv = boost::process::search_path("pigz").string();
  if (argv.empty()) {
    argv = boost::process::search_path("gzip").string();
    if (argv.empty()) {
      throw std::runtime_error("Unable to find pigz or gzip in system PATH");
    }
  }

  assert(compression_lvl >= 1 && compression_lvl <= 9);  // NOLINT
  if (absl::EndsWith(argv, "pigz")) {
    absl::StrAppendFormat(&argv, " -p %d", nthreads);
  }
  absl::StrAppendFormat(&argv, " -%d -c%s", compression_lvl, force_write ? " --force" : "");

  boost::process::opstream compressor_stdin;
  boost::process::ipstream compressor_stderr;
  std::string stderr_msg;
  fmt::print(stderr, "Writing contacts to temporary file '{}'...\n{}\n", path_to_output_file, argv);
  boost::process::child compressor(
      argv,
      boost::process::std_in<compressor_stdin, boost::process::std_out> path_to_output_file.data(),
      boost::process::std_err > compressor_stderr);

  std::thread t([&]() {
    std::string buff;
    while (compressor.running() && std::getline(compressor_stderr, buff)) {
      absl::StrAppend(&stderr_msg, buff, "\n");
    }
  });
  fmt::print(compressor_stdin,
             "## pairs format v1.0\n#columns: readID chr1 pos1 chr2 pos2 strand1 strand2\n");
  std::string line;
  std::vector<std::string_view> toks;
  uint64_t n_contacts = 0;
  for (const auto& path_to_cmatrix : path_to_cmatrices) {
    fmt::print(stderr, "Processing {}...", path_to_cmatrix);
    n_contacts = 0;
    const auto t0 = absl::Now();
    std::ifstream fp_in(path_to_cmatrix, std::ios_base::binary);
    if (!fp_in) {
      throw fmt::system_error(errno, "Unable to open file '{}' for reading", path_to_cmatrix);
    }
    boost::iostreams::filtering_istream in;
    try {
      in.push(boost::iostreams::bzip2_decompressor());
      in.push(fp_in);
      if (!in) {
        throw std::runtime_error("An error occurred while initializing decompression stream");
      }
      auto header = modle::ContactMatrix<uint32_t>::parse_header(path_to_cmatrix, in, false);

      uint64_t n{};
      uint64_t i{};
      for (i = 0UL; std::getline(in, line); ++i) {
        if (!in || !fp_in) {
          throw fmt::system_error(errno, "IO error while reading file '{}'", path_to_cmatrix);
        }
        if (toks = absl::StrSplit(line, '\t'); toks.size() != header.ncols) {
          throw std::runtime_error(fmt::format(
              "File '{}' appears to be malformed: expected {} fields at row {}, but got {}",
              path_to_cmatrix, header.ncols, i, toks.size()));
        }
        const auto chr_start = static_cast<double>(header.start);
        const auto chr_end = static_cast<double>(header.end);
        const auto bin_size = static_cast<double>(header.bin_size);
        for (auto j = 0UL; j < toks.size(); ++j) {
          modle::utils::parse_numeric_or_throw(toks[j], n);
          const double pos1 = std::clamp(
              std::round(2.0 * (chr_start + (static_cast<double>(j - i) * bin_size)) + bin_size) /
                  2.0,
              chr_start, chr_end);
          const double pos2 = std::clamp(
              std::round(2.0 * (chr_start + (static_cast<double>(j) * bin_size)) + bin_size) / 2.0,
              chr_start, chr_end);
          for (auto k = n; k > 0UL; --k) {
            // readID chr1 position1 chr2 position2 strand1 strand2
            fmt::print(compressor_stdin,
                       fmt::format(FMT_COMPILE(". {} {} {} {} . .\n"), header.chr_name,
                                   process_pos(pos1, chr_start, chr_end), header.chr_name,
                                   process_pos(pos2, chr_start, chr_end)));
          }
          n_contacts += n;
        }
      }

      if (i != header.nrows) {
        throw std::runtime_error(fmt::format("Expected {} rows, got {}", header.nrows, i));
      }
      if (!in.eof()) {
        throw std::runtime_error(
            fmt::format("IO error while reading file '{}': {}", path_to_cmatrix,
                        in.fail() ? std::strerror(errno) : "Reading stopped before reaching EOF"));
      }

      fp_in.close();
      if (!fp_in) {
        throw fmt::system_error(errno, "IO error while reading file '{}'", path_to_cmatrix);
      }
    } catch (const boost::iostreams::bzip2_error& err) {
      throw std::runtime_error(fmt::format("An error occurred while decompressing file '{}': {}",
                                           path_to_cmatrix, err.what()));
    }

    fmt::print(stderr, " DONE! Written {} contacts in {}\n", n_contacts,
               absl::FormatDuration(absl::Now() - t0));
  }

  compressor_stdin.close();
  compressor_stdin.pipe().close();  // Required for the child to return
  compressor.wait();
  t.join();
  if (auto ec = compressor.exit_code(); ec != 0) {
    throw std::runtime_error(
        fmt::format("Juicer tools terminated with exit code {}:\n{}", ec, stderr_msg));
  }
}

}  // namespace modle::dcic4dn