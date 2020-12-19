#include "modle_tools/convert.hpp"

#include <absl/strings/str_cat.h>  // for StrAppend, StrCat
#include <absl/time/clock.h>       // for Now
#include <absl/time/time.h>        // for FormatDuration, operator-, Time
#include <fmt/format.h>            // for format, system_error

#include <boost/asio/thread_pool.hpp>
#include <boost/process/child.hpp>
#include <boost/process/io.hpp>    // for std_out_, null, std_err, std_in, std_in_, std_out
#include <boost/process/pipe.hpp>  // for ipstream
#include <cassert>
#include <cstdint>     // for uint*_t
#include <cstdio>      // for stderr
#include <filesystem>  // for remove
#include <memory>      // for unique_ptr, make_unique, make_shared
#include <random>      // for normal_distribution, mt19937_64, generate_canonical
#include <stdexcept>   // for runtime_error
#include <string_view>

#include "modle/4dn_dcic.hpp"
#include "modle/contacts.hpp"  // for ContactMatrix<>::Header, ContactMatrix
#include "modle_tools/config.hpp"  // for config
#include "modle_tools/utils.hpp"   // for generate_random_path_name

namespace modle::tools {

std::normal_distribution<double> init_noise_generator(uint64_t noise_range, double noise_stddev) {
  constexpr double stdev_99_9_ci = 6.0;  // 99.9% CI
  return std::normal_distribution<double>(
      0, noise_stddev == 0 ? static_cast<double>(noise_range) / stdev_99_9_ci : noise_stddev);
}

void convert_to_hic(const modle::tools::config& c, std::string_view template_argv) {
  assert(c.convert_to_hic);  // NOLINT
  try {
    const auto tmp_file_name = modle::tools::utils::generate_random_path(
        c.tmp_dir, ".gz");  // Adding .gz at the end is important, otherwise Juicer Tools will crash
    modle::dcic4dn::converter::from_modle(c.path_to_input_matrices, tmp_file_name,
                                          c.add_noise ? c.noise_range : 0U, c.force, c.noise_stdev,
                                          1, c.nthreads);

    const auto juicer_tools_ver = modle::tools::utils::detect_juicer_tools_version(template_argv);

    std::string argv = fmt::format(
        FMT_STRING("{} pre -t {}{}{} {}.hic {}"), template_argv, c.tmp_dir,
        modle::tools::utils::juicer_tools_version_is_greater_or_equal(juicer_tools_ver, {1, 22, 01})
            ? fmt::format(" -j {} ", c.nthreads)
            : "",
        tmp_file_name, c.output_base_name, c.chr_sizes);

    const auto t0 = absl::Now();
    fmt::print(stderr, "Running Juicer Tools pre. This may take a while...\n{}\n", argv);
    boost::process::ipstream juicer_tools_stderr;
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
    fmt::print(stderr, "DONE! Conversion of {} files to HIC format took {}!\n",
               c.path_to_input_matrices.size(), absl::FormatDuration(absl::Now() - t0));
  } catch (const std::runtime_error& err) {
    throw std::runtime_error(
        fmt::format("An error occurred while converting files to .hic format: {}", err.what()));
  }
}

void convert_to_tsv(boost::asio::thread_pool& tpool, const modle::tools::config& c) {
  for (auto i = 0UL; i < c.path_to_input_matrices.size(); ++i) {
    boost::asio::post(tpool, [&, i]() {
      const auto& path_to_cmatrix = c.path_to_input_matrices[i];
      auto out_file =
          absl::StrCat(c.output_base_name, "/",
                       std::string_view(path_to_cmatrix.data(), path_to_cmatrix.rfind(".tsv.bz2")),
                       "_square_cmatrix.tsv.bz2");
      assert(c.force || !std::filesystem::exists(out_file));
      fmt::print(stderr, "Writing contacts as a symmetric matrix in TSV format to file '{}'...\n",
                 out_file);

      const auto t0 = absl::Now();
      auto header = ContactMatrix<uint32_t>::parse_header(path_to_cmatrix);
      auto noise_generator = c.add_noise ? std::make_unique<std::normal_distribution<double>>(
                                               init_noise_generator(c.noise_range, c.noise_stdev))
                                         : nullptr;

      auto m = ContactMatrix<uint32_t>(path_to_cmatrix, noise_generator.get(), c.seed);
      auto [bytes_in, bytes_out] = m.write_full_matrix_to_tsv(out_file);
      fmt::print(stderr, "DONE! Compression took {} (compression rate {:.2f})\n",
                 absl::FormatDuration(absl::Now() - t0), static_cast<double>(bytes_in) / bytes_out);
    });
  }
}

}  // namespace modle::tools