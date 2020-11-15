#include "modle_tools/eval.hpp"

#include <boost/asio/io_context.hpp>
#include <boost/process.hpp>
#include <future>
#include <thread>

#include "absl/strings/str_format.h"
#include "absl/time/clock.h"
#include "modle/contacts.hpp"
#include "modle_tools/cli.hpp"
#include "modle_tools/utils.hpp"

namespace modle::tools::eval {

ContactMatrix<uint32_t> parse_hic_matrix(const modle::tools::config &c) {
  auto t0 = absl::Now();
  const auto header = modle::ContactMatrix<uint32_t>::parse_header(c.path_to_input_matrix);

  std::string argv = modle::tools::utils::init_juicer_tools_argv(c);
  const auto chr = c.chr_name_hic.empty() ? header.chr_name : c.chr_name_hic;
  const auto start = header.start + ((c.chr_offset_hic != UINT64_MAX) * c.chr_offset_hic);
  const auto end = header.end + ((c.chr_offset_hic != UINT64_MAX) * c.chr_offset_hic) - 1;
  const auto coords = absl::StrCat(chr, ":", start, ":", end);
  absl::StrAppendFormat(&argv, " dump observed NONE %s %s %s BP %lu --",
                        c.path_to_reference_cmatrix, coords, coords, header.bin_size);
  ContactMatrix<uint32_t> cmatrix(header.nrows, header.ncols);
  try {
    boost::asio::io_context ios;
    std::future<std::string> juicer_tools_stderr;
    absl::FPrintF(stderr, "Running %s\n", argv);
    std::string buff, stderr_msg;
    boost::process::ipstream juicer_tools_stdout;
    boost::process::child juicer_tools(argv, boost::process::std_in.close(),
                                       boost::process::std_out > juicer_tools_stdout,
                                       boost::process::std_err > juicer_tools_stderr, ios);
    uint64_t records_parsed = 0;
    auto juicer_stdout_parser = std::thread([&] {
      std::vector<std::string_view> toks;
      std::string buff;

      uint64_t i, j;
      double contacts;
      while (juicer_tools.running() && std::getline(juicer_tools_stdout, buff)) {
        if (absl::StartsWith(buff, "INFO") || absl::StartsWith(buff, "WARN") ||
            absl::StartsWith(buff, "ERROR")) {
          continue;
        }
        if (toks = absl::StrSplit(buff, '\t'); toks.size() != 3) {
          throw std::runtime_error(absl::StrFormat(
              "Malformed file: expected 3 fields, got %lu: line that triggered the error: '%s'",
              toks.size(), buff));
        }
        modle::utils::parse_numeric_or_throw(toks[0], i);
        modle::utils::parse_numeric_or_throw(toks[1], j);
        modle::utils::parse_real_or_throw(toks[2], contacts);
        i = (i - start) / header.bin_size;
        j = (j - start) / header.bin_size;
        cmatrix.set(i, j, std::round<uint64_t>(contacts));
        if (++records_parsed % 1'000'000 == 0) {
          absl::FPrintF(stderr, "Parsed %lu records...\n", records_parsed);
        }
      }

      if (!juicer_tools_stdout && !juicer_tools_stdout.eof()) {
        throw std::runtime_error(absl::StrCat(
            "An error occurred while reading from Juicer Tools stdout: ", std::strerror(errno)));
      }
    });
    ios.run();
    juicer_tools.wait();
    juicer_stdout_parser.join();
    if (auto ec = juicer_tools.exit_code(); ec != 0) {
      throw std::runtime_error(absl::StrFormat("Juicer Tools terminated with exit code %lu: %s", ec,
                                               juicer_tools_stderr.get()));
    }
    absl::FPrintF(stderr, "Parsed %lu records in %s!\n", records_parsed,
                  absl::FormatDuration(absl::Now() - t0));
  } catch (const std::runtime_error &err) {
    throw std::runtime_error(
        absl::StrFormat("An error occurred while running juicer_tools dump: %s", err.what()));
  }

  return cmatrix;
}

}  // namespace modle::utils::eval
