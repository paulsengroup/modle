#include "modle/juicer_contacts.hpp"

#include <absl/strings/match.h>  // for StartsWith
#include <absl/strings/str_cat.h>
#include <absl/strings/str_format.h>  // for StrAppendFormat
#include <absl/strings/str_split.h>
#include <absl/time/clock.h>  // for Now
#include <absl/time/time.h>   // for FormatDuration, ToDoubleSeconds, Time
#include <fmt/format.h>       // for format, print, system_error

#include <boost/asio/io_context.hpp>
#include <boost/process/child.hpp>
#include <boost/process/io.hpp>    // for std_err, std_in, std_out
#include <boost/process/pipe.hpp>  // for ipstream
#include <cassert>
#include <cerrno>
#include <cmath>       // for lround
#include <cstdio>      // for stderr
#include <filesystem>  // for exists, is_directory
#include <future>
#include <new>
#include <stdexcept>  // for runtime_error, logic_error
#include <string>     // for string, getline
#include <string_view>
#include <thread>
#include <type_traits>  // for declval, is_arithmetic, is_integral
#include <utility>      // for forward, move
#include <vector>

#include "modle/contacts.hpp"  // for ContactMatrix<>::Header, ContactMatrix

namespace modle::juicer_contacts {

template <typename I, typename R>
Contact::Contact(I b1, I b2, R contacts_num) {
  static_assert(std::is_integral<I>());
  static_assert(std::is_arithmetic<R>());
  if (b1 < 0) {
    throw std::runtime_error(
        fmt::format("Exception thrown in Contact constructor: b1 should be a positive integral "
                    "number, is '{}'",
                    b1));
  }
  if (b2 < 0) {
    throw std::runtime_error(
        fmt::format("Exception thrown in Contact constructor: b2 should be a positive integral "
                    "number, is '{}'",
                    b2));
  }
  if (contacts_num < 0) {
    throw std::runtime_error(fmt::format(
        "Exception thrown in Contact constructor: b2 should be a positive number, is '{}'",
        contacts_num));
  }
  this->bin1 = static_cast<uint64_t>(b1);
  this->bin2 = static_cast<uint64_t>(b2);
  this->contacts = static_cast<double>(contacts_num);
}

ContactsParser::ContactsParser(std::string contact_file) : _path(std::move(contact_file)) {
  if (!std::filesystem::exists(this->_path)) {
    throw std::runtime_error(fmt::format("File '{}' does not exist", this->_path));
  }
  if (std::filesystem::is_directory(this->_path)) {
    throw std::runtime_error(fmt::format("'{}' is a directory. Expected a file", this->_path));
  }
}

bool ContactsParser::parse_next(std::string& buff, Contact& c, std::string_view sep) {
  assert(this->_f);
  if (std::getline(this->_f, buff); buff.empty()) {
    if (this->_f.eof()) {
      return false;
    }
    throw std::logic_error(
        "If you see this error, please report it to the developers on GitHub: "
        "ContactsParser::parse_next() "
        "reached a branch that should not be reachable!");
  }
  if (!this->_f && !this->_f.eof()) {
    throw fmt::system_error(errno, "IO error while reading file {}", this->_path);
  }
  std::vector<std::string_view> tmp = absl::StrSplit(buff, sep);
  if (tmp.size() != 3) {
    throw std::runtime_error(
        fmt::format("Malformed record in file '{}': expected 3 fields, got {} in line '{}'",
                    this->_path, tmp.size(), buff));
  }
  utils::parse_numeric_or_throw(tmp[0], c.bin1);
  utils::parse_numeric_or_throw(tmp[1], c.bin2);
  utils::parse_real_or_throw(tmp[2], c.contacts);
  return false;
}

ContactsParser::MatrixProperties ContactsParser::get_matrix_properties(std::string_view sep) {
  Contact c{};
  ContactsParser::MatrixProperties prop{};
  std::string buff;
  while (this->parse_next(buff, c, sep)) {
    prop.min_bin1 = std::min(prop.min_bin1, c.bin1);
    prop.min_bin2 = std::min(prop.min_bin2, c.bin2);
    prop.max_bin1 = std::max(prop.max_bin1, c.bin1);
    prop.max_bin2 = std::max(prop.max_bin2, c.bin2);
    if (prop.bin_size == 0) {
      prop.bin_size = c.bin1 > c.bin2 ? c.bin1 - c.bin2 : c.bin2 - c.bin1;
      continue;
    }
    if (uint64_t new_bin_size = c.bin1 > c.bin2 ? c.bin1 - c.bin2 : c.bin2 - c.bin1;
        new_bin_size != prop.bin_size) {
      throw std::runtime_error(
          fmt::format("Expected bin of size {}, got {}", prop.bin_size, new_bin_size));
    }
  }
  if (!this->_f && !this->_f.eof()) {
    throw fmt::system_error(errno, "IO error while reading file {}", this->_path);
  }
  this->_f.clear();
  this->_f.seekg(0);
  return prop;
}

ContactMatrix<uint32_t> ContactsParser::parse_into_contact_matrix(uint64_t width,
                                                                  std::string_view sep) {
  std::string buff;
  Contact c{};
  const auto matrix_prop = this->get_matrix_properties(sep);
  ContactMatrix<uint32_t> m((matrix_prop.max_bin1 - matrix_prop.min_bin1) / matrix_prop.bin_size,
                            width / matrix_prop.bin_size);
  this->_f = std::ifstream(this->_path);
  if (!this->_f) {
    throw fmt::system_error(errno, "Unable to open file '{}' for reading", this->_path);
  }
  try {
    while (this->parse_next(buff, c, sep)) {
      uint64_t row = (c.bin1 - matrix_prop.min_bin1) / matrix_prop.bin_size;
      uint64_t col = (c.bin2 - matrix_prop.min_bin2) / matrix_prop.bin_size;
      assert(c.contacts >= 0);
      m.set(row, col, static_cast<uint32_t>(std::lround(c.contacts)));
    }
    if (!this->_f.eof()) {
      throw std::runtime_error("An IO error occurred");
    }
  } catch (const std::exception& err) {
    throw std::runtime_error(
        fmt::format("An error occurred while reading file '{}': {}", this->_path, err.what()));
  }
  return m;
}

ContactMatrix<uint32_t> run_juicer_dump_and_parse_contacts(
    std::string_view path_to_input_matrix, std::string_view path_to_reference_matrix,
    std::string_view chr_name_hic, uint64_t chr_offset_hic, std::string_view path_to_juicer_tools,
    uint64_t juicer_tools_mem) {
  auto t0 = absl::Now();
  const auto header = modle::ContactMatrix<uint32_t>::parse_header(path_to_input_matrix);

  std::string argv = modle::utils::init_juicer_tools_argv(path_to_juicer_tools, juicer_tools_mem);
  const auto chr = chr_name_hic.empty() ? header.chr_name : chr_name_hic;
  const auto start = chr_offset_hic != UINT64_MAX ? chr_offset_hic + header.start : header.start;
  const auto end = chr_offset_hic != UINT64_MAX ? chr_offset_hic + header.end : header.end;
  const auto coords = absl::StrCat(chr, ":", start, ":", end - 1);
  absl::StrAppendFormat(&argv, " dump observed NONE %s %s %s BP %lu --", path_to_reference_matrix,
                        coords, coords, header.bin_size);
  ContactMatrix<uint32_t> cmatrix(header.nrows, header.ncols);
  try {
    boost::asio::io_context ios;
    std::future<std::string> juicer_tools_stderr;
    fmt::print(stderr, "Running {}\n", argv);
    std::string buff;
    std::string stderr_msg;
    boost::process::ipstream juicer_tools_stdout;
    boost::process::child juicer_tools(argv, boost::process::std_in.close(),
                                       boost::process::std_out > juicer_tools_stdout,
                                       boost::process::std_err > juicer_tools_stderr, ios);
    uint64_t records_parsed = 0;
    auto juicer_stdout_parser = std::thread([&] {
      auto t1 = absl::Now();
      std::vector<std::string_view> toks;

      uint64_t i{};
      uint64_t j{};
      double contacts{};
      while (juicer_tools.running() && std::getline(juicer_tools_stdout, buff)) {
        if (absl::StartsWith(buff, "INFO") || absl::StartsWith(buff, "WARN") ||
            absl::StartsWith(buff, "ERROR")) {
          continue;
        }
        if (toks = absl::StrSplit(buff, '\t'); toks.size() != 3) {
          throw std::runtime_error(fmt::format(
              "Malformed file: expected 3 fields, got {}: line that triggered the error: '{}'",
              toks.size(), buff));
        }
        modle::utils::parse_numeric_or_throw(toks[0], i);
        modle::utils::parse_numeric_or_throw(toks[1], j);
        modle::utils::parse_real_or_throw(toks[2], contacts);
        i = (i - start) / header.bin_size;
        j = (j - start) / header.bin_size;
        cmatrix.set(i, j, std::lround(contacts));
        // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
        if (auto step = 1'000'000U; ++records_parsed % step == 0U) {
          fmt::print(stderr, "Parsed {} records ({:.2f} records/s)...\n", records_parsed,
                     static_cast<double>(step) / absl::ToDoubleSeconds(absl::Now() - t1));
        }
      }

      if (!juicer_tools_stdout && !juicer_tools_stdout.eof()) {
        throw fmt::system_error(errno, "An error occurred while reading Juicer Tools stdout");
      }
    });

    ios.run();
    juicer_tools.wait();
    juicer_stdout_parser.join();
    if (auto ec = juicer_tools.exit_code(); ec != 0) {
      throw std::runtime_error(fmt::format("Juicer Tools terminated with exit code {}: {}", ec,
                                           juicer_tools_stderr.get()));
    }
    fmt::print(stderr, "Parsed {} records in {}!\n", records_parsed,
               absl::FormatDuration(absl::Now() - t0));
  } catch (const std::runtime_error& err) {
    throw std::runtime_error(
        fmt::format("An error occurred while running juicer_tools dump: {}", err.what()));
  }

  return cmatrix;
}

}  // namespace modle::juicer_contacts