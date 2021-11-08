// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <absl/time/clock.h>                      // for Now
#include <absl/time/time.h>                       // for FormatDuration, operator-, Time
#include <fmt/format.h>                           // for format, make_format_args, vfor...
#include <readerwriterqueue/readerwriterqueue.h>  // for BlockingReaderWriterQueue
#include <spdlog/spdlog.h>                        // for info

#include <algorithm>    // for clamp, min, max, max_element
#include <cassert>      // for assert
#include <chrono>       // for milliseconds
#include <cmath>        // for round
#include <cstdio>       // for usize
#include <exception>    // for exception
#include <functional>   // for ref, hash, reference_wrapper
#include <limits>       // for numeric_limits
#include <stdexcept>    // for runtime_error
#include <string>       // for basic_string
#include <string_view>  // for hash, string_view
#include <thread>       // for sleep_for, thread
#include <type_traits>  // for add_const<>::type
#include <utility>      // for tuple_element<>::type
#include <vector>       // for vector

#include "modle/common/genextreme_value_distribution.hpp"  // for genextreme_value_distribution
#include "modle/common/random.hpp"                         // for PRNG
#include "modle/contacts.hpp"                              // for ContactMatrix
#include "modle/cooler/cooler.hpp"                         // for Cooler::Pixel, Cooler, Cooler:...
#include "modle_tools/config.hpp"                          // for noisify_config
#include "modle_tools/tools.hpp"                           // for noisify_subcmd

namespace modle::tools {

void noisify_contacts(const noisify_config& c) {
  using pixel_queue_t = moodycamel::BlockingReaderWriterQueue<modle::cooler::Cooler<>::Pixel>;

  const auto PIXEL_BATCH_SIZE =
      modle::cooler::Cooler<>::DEFAULT_HDF5_BUFFER_SIZE / sizeof(modle::cooler::Cooler<>::Pixel);
  pixel_queue_t pixel_queue(PIXEL_BATCH_SIZE);
  modle::ContactMatrix<> cmatrix{};

  const auto END_OF_PIXEL_QUEUE = modle::cooler::Cooler<>::Pixel{
      (std::numeric_limits<contacts_t>::max)(), (std::numeric_limits<contacts_t>::max)(),
      (std::numeric_limits<contacts_t>::max)()};

  auto input_cool = modle::cooler::Cooler(c.path_to_input_matrix,
                                          cooler::Cooler<>::IO_MODE::READ_ONLY, c.bin_size);

  const auto bin_size = input_cool.get_bin_size();
  const auto max_chrom_name_size = [&]() {
    const auto& chrom_names = input_cool.get_chrom_names();
    const auto it =
        std::max_element(chrom_names.begin(), chrom_names.end(),
                         [](const auto& c1, const auto& c2) { return c1.size() < c2.size(); });
    return it->size();
  }();

  modle::cooler::Cooler output_cool(c.path_to_output_matrix, cooler::Cooler<>::IO_MODE::WRITE_ONLY,
                                    input_cool.get_bin_size(), max_chrom_name_size);

  for (const auto& [chrom_name, chrom_size] : input_cool.get_chroms()) {
    const auto t0 = absl::Now();
    spdlog::info(FMT_STRING("Processing contacts for {}..."), chrom_name);
    const auto ncols = (chrom_size / bin_size) + static_cast<usize>(chrom_size % bin_size != 0);
    const auto nrows = std::min(ncols, (c.diagonal_width / bin_size) +
                                           static_cast<usize>(c.diagonal_width % bin_size != 0));

    std::thread t([&, chrom_name = chrom_name]() {
      input_cool.stream_contacts_for_chrom(std::ref(pixel_queue), chrom_name, c.diagonal_width,
                                           bin_size);
      while (!pixel_queue.try_enqueue(END_OF_PIXEL_QUEUE)) {
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
      }
    });

    cmatrix.unsafe_resize(nrows, ncols);
    cmatrix.unsafe_reset();
    modle::cooler::Cooler<>::Pixel pixel{};  // NOLINT
    const auto seed =
        c.seed + std::hash<std::string_view>{}(chrom_name) + std::hash<usize>{}(ncols);

    auto rang_eng = random::PRNG(seed);
    auto genextreme = modle::genextreme_value_distribution<double>{
        c.genextreme_mu, c.genextreme_sigma, c.genextreme_xi};
    while (true) {
      pixel_queue.wait_dequeue(pixel);
      if (pixel == END_OF_PIXEL_QUEUE) {
        break;
      }
      assert(pixel.row <= pixel.col);  // NOLINT
      for (usize i = 0; i < pixel.count; ++i) {
        const auto pos1 = static_cast<double>(pixel.row * bin_size) + genextreme(rang_eng);
        const auto pos2 = static_cast<double>(pixel.col * bin_size) - genextreme(rang_eng);
        const auto bin1 = std::clamp(
            static_cast<usize>(std::round(pos1 / static_cast<double>(bin_size))), 0UL, ncols - 1);
        const auto bin2 = std::clamp(
            static_cast<usize>(std::round(pos2 / static_cast<double>(bin_size))), 0UL, ncols - 1);
        cmatrix.increment(std::min(bin1, bin2), std::max(bin1, bin2));
      }
    }
    t.join();

    output_cool.write_or_append_cmatrix_to_file(cmatrix, chrom_name, 0UL, chrom_size, chrom_size);
    spdlog::info(FMT_STRING("DONE procesing contacts for {} in {}."), chrom_name,
                 absl::FormatDuration(absl::Now() - t0));
  }
}

void noisify_subcmd(const modle::tools::noisify_config& c) { modle::tools::noisify_contacts(c); }

}  // namespace modle::tools
