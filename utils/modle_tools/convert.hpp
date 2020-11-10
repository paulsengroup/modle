#pragma once

#include <algorithm>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/process.hpp>
#include <filesystem>
#include <random>

#include "absl/time/clock.h"
#include "cli.hpp"
#include "modle/contacts.hpp"
#include "utils.hpp"

namespace modle::tools::convert {

std::normal_distribution<float> init_noise_generator(uint64_t range) {
  return std::normal_distribution<float>(0, range / 6.0);  // 99.9 CI
}

void check_convert_preconditions(modle::utils::config& c, std::string& argv) {
  if (c.convert_to_hic) {
    if (c.path_to_juicer.empty()) {
      for (std::string file : {"juicer_tools", "juicer_tools.sh", "juicer_tools.exe"}) {
        if (auto p = boost::process::search_path("juicer_tools").string(); !p.empty()) {
          c.path_to_juicer = std::move(p);
          break;
        }
      }
    }
    if (c.path_to_juicer.empty()) {
      throw std::runtime_error(
          "--path-to-juicer-tools was not specified and we were unable to find Juicer "
          "tools in your path");
    }
    argv = c.path_to_juicer;
    if (absl::EndsWith(argv, ".jar")) {
      // TODO: Check that java >= 1.7
      auto java = boost::process::search_path("java").string();
      if (java.empty())
        throw std::runtime_error(
            "--path-to-juicer-tools points to a jar file, but we were unable to find java in "
            "your "
            "path");
      argv = absl::StrCat(java, " -Xms512m -Xmx2048m -jar ", c.path_to_juicer);
    }
  }
}

std::string prepare_files_for_juicer_tools(modle::utils::config& c) {
  assert(c.convert_to_hic);
  auto tmp_file_name = modle::tools::utils::generate_random_path_name(c.tmp_dir);
  absl::StrAppend(&tmp_file_name, ".txt.gz");
  absl::FPrintF(stderr, "Writing contacts to temporary file '%s'...\n", tmp_file_name);
  std::ifstream fp_in(c.path_to_input_matrix, std::ios_base::binary);
  if (!fp_in)
    throw std::runtime_error(absl::StrFormat("Unable to open file '%s' for reading: %s",
                                             c.path_to_input_matrix, std::strerror(errno)));
  std::ofstream fp_out(tmp_file_name, std::ios_base::binary);
  if (!fp_out)
    throw std::runtime_error(absl::StrFormat("Unable to create temporary file '%s': %s",
                                             tmp_file_name, std::strerror(errno)));
  boost::iostreams::filtering_istream in;
  boost::iostreams::filtering_ostream out;
  std::string line;
  std::vector<std::string_view> toks;
  try {
    in.push(boost::iostreams::bzip2_decompressor());
    in.push(fp_in);
    out.push(boost::iostreams::gzip_compressor());
    out.push(fp_out);
    if (!in) throw std::runtime_error("An error occurred while initializing decompression stream");
    if (!out) throw std::runtime_error("An error occurred while initializing compression stream");

    auto header = modle::ContactMatrix<uint32_t>::parse_header(c.path_to_input_matrix, in, false);
    std::unique_ptr<std::normal_distribution<float>> noise_generator = nullptr;
    std::unique_ptr<std::mt19937_64> rand_eng = nullptr;
    if (c.add_noise) {
      noise_generator =
          std::make_unique<std::normal_distribution<float>>(init_noise_generator(c.noise_range));
      rand_eng = std::make_unique<std::mt19937_64>(c.seed);
    }
    uint64_t n, i, n_contacts = 0;
    std::string record;
    auto t0 = absl::Now();
    for (i = 0UL; std::getline(in, line); ++i) {
      if (!in || !fp_in) {
        if (in.eof() && fp_in.eof()) break;  // This should never happen
        throw std::runtime_error(absl::StrFormat("IO error while reading from file '%s': %s",
                                                 c.path_to_input_matrix, std::strerror(errno)));
      }
      if (toks = absl::StrSplit(line, '\t'); toks.size() != header.ncols) {
        throw std::runtime_error(absl::StrFormat(
            "File '%s' appears to be malformed: expected %lu in row %lu, but got %lu",
            c.path_to_input_matrix, header.ncols, i, toks.size()));
      }
      for (auto j = 0UL; j < toks.size(); ++j) {
        modle::utils::parse_numeric_or_throw(toks[j], n);
        uint64_t pos1 = std::clamp(static_cast<uint64_t>(std::round(
                                       ((2 * (j - i) * header.bin_size) + header.bin_size) / 2.0)),
                                   0UL, header.end);
        uint64_t pos2 = std::clamp(
            static_cast<uint64_t>(std::round(((2 * j * header.bin_size) + header.bin_size) / 2.0)),
            0UL, header.end);
        record = absl::StrFormat("1 %s %lu 0 0 %s %lu 1\n", header.chr_name, pos1, header.chr_name,
                                 pos2);
        for (auto k = n; k > 0UL; --k) {
          if (c.add_noise) {
            record = absl::StrFormat(
                "1 %s %lu 0 0 %s %lu 1\n", header.chr_name,
                std::clamp<uint64_t>(std::round((*noise_generator)(*rand_eng) + pos1), 0.0,
                                     static_cast<float>(header.end)),
                header.chr_name,
                std::clamp<uint64_t>(std::round((*noise_generator)(*rand_eng) + pos2), 0.0,
                                     static_cast<float>(header.end)));
          }
          out.write(record.data(), record.size());
          if (!out || !fp_out) {
            throw std::runtime_error(
                absl::StrFormat("IO error while writing to temporary file '%s': %s", tmp_file_name,
                                std::strerror(errno)));
          }
        }
        n_contacts += n;
      }
    }

    if (i != header.nrows) {
      throw std::runtime_error(absl::StrFormat("Expected %lu rows, got %lu", header.nrows, i));
    }
    if (!in.eof()) {
      throw std::runtime_error(absl::StrFormat(
          "IO error while reading from file '%s': %s", c.path_to_input_matrix,
          in.fail() ? std::strerror(errno) : "Reading stopped before reaching EOF"));
    }
    fp_in.close();
    if (!fp_in) {
      throw std::runtime_error(absl::StrFormat("IO error while reading from file '%s': %s",
                                               c.path_to_input_matrix, std::strerror(errno)));
    }
    if (!out || !fp_out) {
      throw std::runtime_error(absl::StrFormat("IO error while writing to temporary file '%s': %s",
                                               c.path_to_input_matrix, std::strerror(errno)));
    }
    out.reset();  // Very important!! Needed to actually call gzip_compressor::close(). The stream
    // won't be usable after this
    fp_out.close();
    if (!fp_out || !out)
      throw std::runtime_error(absl::StrFormat("IO error while closing file '%s': %s",
                                               tmp_file_name, std::strerror(errno)));
    absl::FPrintF(stderr, "DONE! Written %lu contacts in %s!\n", n_contacts,
                  absl::FormatDuration(absl::Now() - t0));
  } catch (const boost::iostreams::bzip2_error& err) {
    throw std::runtime_error(absl::StrFormat("An error occurred while decompressing file '%s': %s",
                                             c.path_to_input_matrix, err.what()));
  } catch (const boost::iostreams::gzip_error& err) {
    throw std::runtime_error(absl::StrFormat(
        "An error occurred while compressing temporary file '%s': %s", tmp_file_name, err.what()));
  }

  return tmp_file_name;
}

void convert_to_hic(modle::utils::config& c, std::string& argv) {
  assert(c.convert_to_hic);
  try {
    auto tmp_file_name = prepare_files_for_juicer_tools(c);
    absl::StrAppendFormat(&argv, " pre %s %s/%s.hic %s", tmp_file_name, c.out_dir, c.base_name,
                          c.chr_sizes);
    const auto t0 = absl::Now();
    absl::FPrintF(stderr, "Running: %s\n", argv);
    boost::process::ipstream juicer_tools_stderr;
    std::string buff, stderr_msg, line;
    boost::process::child juicer_tools(argv, boost::process::std_in.close(),
                                       boost::process::std_out > boost::process::null,
                                       boost::process::std_err > juicer_tools_stderr);

    while (juicer_tools.running() && std::getline(juicer_tools_stderr, line)) {
      absl::StrAppendFormat(&stderr_msg, "%s\n", line);
    }
    juicer_tools.wait();
    if (!c.keep_tmp_files) std::filesystem::remove(tmp_file_name);
    if (auto ec = juicer_tools.exit_code(); ec != 0) {
      throw std::runtime_error(
          absl::StrFormat("Juicer tools terminated with exit code %lu: %s", ec, stderr_msg));
    }
    absl::FPrintF(stderr, "DONE! Conversion took %s!\n", absl::FormatDuration(absl::Now() - t0));
  } catch (const std::runtime_error& err) {
    throw std::runtime_error(
        absl::StrFormat("An error occurred while converting file '%s' to .hic format: %s",
                        c.path_to_input_matrix, err.what()));
  }
}

void convert_to_tsv(modle::utils::config& c) {
  auto out_file = absl::StrCat(c.out_dir, "/", c.base_name, ".tsv.bz2");
  absl::FPrintF(stderr, "Writing contacts as a symmetric matrix in TSV format to file '%s'...\n",
                out_file);

  const auto t0 = absl::Now();
  auto header = ContactMatrix<uint32_t>::parse_header(c.path_to_input_matrix);
  auto noise_generator =
      c.add_noise
          ? std::make_unique<std::normal_distribution<float>>(init_noise_generator(c.noise_range))
          : nullptr;

  auto m = ContactMatrix<uint32_t>(c.path_to_input_matrix, noise_generator.get(), c.seed);
  auto [bytes_in, bytes_out] = m.write_full_matrix_to_tsv(out_file);
  absl::FPrintF(stderr, "DONE! Compression took %s (compression rate %.2f)\n",
                absl::FormatDuration(absl::Now() - t0), static_cast<double>(bytes_in) / bytes_out);
}
}  // namespace modle::tools::convert