#include "modle_tools/convert.hpp"

#include <absl/strings/str_cat.h>     // for StrAppend, StrCat
#include <absl/strings/str_format.h>  // for StrAppendFormat
#include <absl/strings/str_split.h>
#include <absl/time/clock.h>  // for Now
#include <absl/time/time.h>   // for FormatDuration, operator-, Time
#include <fmt/compile.h>      // for FMT_COMPILE
#include <fmt/format.h>       // for format, system_error

#include <algorithm>                         // for clamp, max
#include <boost/iostreams/filter/bzip2.hpp>  // for basic_bzip2_decompressor, bzip2_error, bzip2_decompressor
#include <boost/iostreams/filter/gzip.hpp>  // for basic_gzip_compressor, gzip_error, gzip_compressor
#include <boost/iostreams/filtering_stream.hpp>  // for filtering_ostream, filtering_istream
#include <boost/process/child.hpp>
#include <boost/process/io.hpp>    // for std_out_, null, std_err, std_in, std_in_, std_out
#include <boost/process/pipe.hpp>  // for ipstream
#include <cassert>
#include <cerrno>
#include <cmath>       // for round
#include <cstdint>     // for uint*_t
#include <cstdio>      // for stderr
#include <cstring>     // for strerror
#include <filesystem>  // for remove
#include <istream>     // for basic_ios, ofstream, basic_istream, ifstream, ios_base
#include <memory>      // for unique_ptr, make_unique, make_shared
#include <random>      // for normal_distribution, mt19937_64, generate_canonical
#include <stdexcept>   // for runtime_error
#include <string_view>
#include <vector>

#include "modle/contacts.hpp"     // for ContactMatrix<>::Header, ContactMatrix
#include "modle_tools/cli.hpp"    // for config
#include "modle_tools/utils.hpp"  // for generate_random_path_name

namespace modle::tools {
std::normal_distribution<double> init_noise_generator(uint64_t range) {
  constexpr double stdev_99_9_ci = 6.0;  // 99.9% CI
  return std::normal_distribution<double>(0, static_cast<double>(range) / stdev_99_9_ci);
}

std::string prepare_files_for_juicer_tools(const modle::tools::config& c) {
  assert(c.convert_to_hic);
  auto tmp_file_name = modle::tools::utils::generate_random_path_name(c.tmp_dir);
  absl::StrAppend(&tmp_file_name, ".txt.gz");
  fmt::print(stderr, "Writing contacts to temporary file '{}'...\n", tmp_file_name);
  std::ifstream fp_in(c.path_to_input_matrix, std::ios_base::binary);
  if (!fp_in) {
    throw fmt::system_error(errno, "Unable to open file '{}' for reading", c.path_to_input_matrix);
  }
  std::ofstream fp_out(tmp_file_name, std::ios_base::binary);
  if (!fp_out) {
    throw fmt::system_error(errno, "Unable to create temporary file '{}'", tmp_file_name);
  }
  boost::iostreams::filtering_istream in;
  boost::iostreams::filtering_ostream out;
  std::string line;
  std::vector<std::string_view> toks;
  try {
    in.push(boost::iostreams::bzip2_decompressor());
    in.push(fp_in);
    out.push(boost::iostreams::gzip_compressor());
    out.push(fp_out);
    if (!in) {
      throw std::runtime_error("An error occurred while initializing decompression stream");
    }
    if (!out) {
      throw std::runtime_error("An error occurred while initializing compression stream");
    }

    auto header = modle::ContactMatrix<uint32_t>::parse_header(c.path_to_input_matrix, in, false);
    std::unique_ptr<std::normal_distribution<double>> noise_generator = nullptr;
    std::unique_ptr<std::mt19937_64> rand_eng = nullptr;
    if (c.add_noise) {
      noise_generator =
          std::make_unique<std::normal_distribution<double>>(init_noise_generator(c.noise_range));
      rand_eng = std::make_unique<std::mt19937_64>(c.seed);
    }
    uint64_t n{};
    uint64_t i{};
    uint64_t n_contacts = 0;
    std::string record;
    const auto t0 = absl::Now();
    for (i = 0UL; std::getline(in, line); ++i) {
      if (!in || !fp_in) {
        throw fmt::system_error(errno, "IO error while reading file '{}'", c.path_to_input_matrix);
      }
      if (toks = absl::StrSplit(line, '\t'); toks.size() != header.ncols) {
        throw std::runtime_error(fmt::format(
            "File '{}' appears to be malformed: expected {} fields at row {}, but got {}",
            c.path_to_input_matrix, header.ncols, i, toks.size()));
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
        record = fmt::format(FMT_COMPILE("1 {} {:.0f} 0 0 {} {:.0f} 1\n"), header.chr_name, pos1,
                             header.chr_name, pos2);
        for (auto k = n; k > 0UL; --k) {
          if (c.add_noise) {
            record = fmt::format(
                FMT_COMPILE("1 {} {:.0f} 0 0 {} {:.0f} 1\n"), header.chr_name,
                std::clamp(std::round(pos1 + (*noise_generator)(*rand_eng)), chr_start, chr_end),
                header.chr_name,
                std::clamp(std::round(pos2 + (*noise_generator)(*rand_eng)), chr_start, chr_end));
          }
          fmt::print(out, record);
          if (!out || !fp_out) {
            throw fmt::system_error(errno, "IO error while writing to temporary file '{}'",
                                    tmp_file_name);
          }
        }
        n_contacts += n;
      }
    }

    if (i != header.nrows) {
      throw std::runtime_error(fmt::format("Expected {} rows, got {}", header.nrows, i));
    }
    if (!in.eof()) {
      throw std::runtime_error(
          fmt::format("IO error while reading file '{}': {}", c.path_to_input_matrix,
                      in.fail() ? std::strerror(errno) : "Reading stopped before reaching EOF"));
    }
    fp_in.close();
    if (!fp_in) {
      throw fmt::system_error(errno, "IO error while reading file '{}'", c.path_to_input_matrix);
    }
    if (!out || !fp_out) {
      throw fmt::system_error(errno, "IO error while writing to temporary file '{}'",
                              c.path_to_input_matrix);
    }
    out.reset();  // Very important!! Needed to actually call gzip_compressor::close(). The stream
    // won't be usable after this
    fp_out.close();
    if (!fp_out || !out) {
      throw fmt::system_error(errno, "IO error while closing file '{}'", tmp_file_name);
    }
    fmt::fprintf(stderr, "DONE! Written %lu contacts in %s!\n", n_contacts,
                 absl::FormatDuration(absl::Now() - t0));
  } catch (const boost::iostreams::bzip2_error& err) {
    throw std::runtime_error(fmt::format("An error occurred while decompressing file '{}': {}",
                                         c.path_to_input_matrix, err.what()));
  } catch (const boost::iostreams::gzip_error& err) {
    throw std::runtime_error(fmt::format(
        "An error occurred while compressing temporary file '{}': {}", tmp_file_name, err.what()));
  }

  return tmp_file_name;
}

void convert_to_hic(const modle::tools::config& c, std::string& argv) {
  assert(c.convert_to_hic);
  try {
    auto tmp_file_name = prepare_files_for_juicer_tools(c);
    absl::StrAppendFormat(&argv, " pre -t %s %s %s/%s.hic %s", c.tmp_dir, tmp_file_name, c.out_dir,
                          c.base_name, c.chr_sizes);
    const auto t0 = absl::Now();
    fmt::fprintf(stderr, "Running: %s\n", argv);
    boost::process::ipstream juicer_tools_stderr;
    std::string buff;
    std::string stderr_msg;
    std::string line;
    boost::process::child juicer_tools(argv, boost::process::std_in.close(),
                                       boost::process::std_out > boost::process::null,
                                       boost::process::std_err > juicer_tools_stderr);

    while (juicer_tools.running() && std::getline(juicer_tools_stderr, line)) {
      absl::StrAppend(&stderr_msg, line, "\n");
    }
    if (!juicer_tools_stderr && !juicer_tools_stderr.eof()) {
      juicer_tools.terminate();
      throw std::runtime_error("An error occurred while reading from Juicer Tools stderr");
    }
    juicer_tools.wait();
    if (!c.keep_tmp_files) {
      std::filesystem::remove(tmp_file_name);
    }
    if (auto ec = juicer_tools.exit_code(); ec != 0) {
      throw std::runtime_error(
          fmt::format("Juicer tools terminated with exit code {}:\n{}", ec, stderr_msg));
    }
    fmt::fprintf(stderr, "DONE! Conversion took %s!\n", absl::FormatDuration(absl::Now() - t0));
  } catch (const std::runtime_error& err) {
    throw std::runtime_error(
        fmt::format("An error occurred while converting file '{}' to .hic format: {}",
                    c.path_to_input_matrix, err.what()));
  }
}

void convert_to_tsv(const modle::tools::config& c) {
  auto out_file = absl::StrCat(c.out_dir, "/", c.base_name, "_square_cmatrix.tsv.bz2");
  assert(c.force || !std::filesystem::exists(out_file));
  fmt::fprintf(stderr, "Writing contacts as a symmetric matrix in TSV format to file '%s'...\n",
               out_file);

  const auto t0 = absl::Now();
  auto header = ContactMatrix<uint32_t>::parse_header(c.path_to_input_matrix);
  auto noise_generator =
      c.add_noise
          ? std::make_unique<std::normal_distribution<double>>(init_noise_generator(c.noise_range))
          : nullptr;

  auto m = ContactMatrix<uint32_t>(c.path_to_input_matrix, noise_generator.get(), c.seed);
  auto [bytes_in, bytes_out] = m.write_full_matrix_to_tsv(out_file);
  fmt::fprintf(stderr, "DONE! Compression took %s (compression rate %.2f)\n",
               absl::FormatDuration(absl::Now() - t0), static_cast<double>(bytes_in) / bytes_out);
}

}  // namespace modle::tools