#pragma once

// IWYU pragma: private, include "modle_tools/noisify_contacts.hpp"

#include <absl/time/clock.h>  // for Now
#include <absl/time/time.h>   // for FormatDuration
#include <fmt/format.h>       // for print
#include <readerwriterqueue/readerwriterqueue.h>

#include <boost/asio/post.hpp>         // IWYU pragma: keep for post
#include <boost/asio/thread_pool.hpp>  // for thread_pool
#include <random>                      // for gamma_distribution, uniform_distribution
#include <thread>                      // for thread, this_thread::sleep_for

#include "include/modle_tools/config.hpp"
#include "modle/common.hpp"
#include "modle/contacts.hpp"
#include "modle/cooler.hpp"

namespace modle::tools {

void noisify_contacts(const config& c) {
  using pixel_queue_t = moodycamel::BlockingReaderWriterQueue<modle::cooler::Cooler::Pixel>;

  constexpr auto PIXEL_BATCH_SIZE =
      modle::cooler::Cooler::DEFAULT_HDF5_BUFFER_SIZE / sizeof(modle::cooler::Cooler::Pixel);
  pixel_queue_t pixel_queue(PIXEL_BATCH_SIZE);
  modle::ContactMatrix<uint32_t> cmatrix{};

  constexpr auto END_OF_PIXEL_QUEUE = modle::cooler::Cooler::Pixel{
      std::numeric_limits<std::size_t>::max(), std::numeric_limits<std::size_t>::max(),
      std::numeric_limits<std::size_t>::max()};

  auto input_cool =
      modle::cooler::Cooler(c.path_to_input_matrix, cooler::Cooler::READ_ONLY, c.bin_size);

  const auto bin_size = input_cool.get_bin_size();
  const auto max_chrom_name_size = [&]() {
    const auto& chrom_names = input_cool.get_chrom_names();
    const auto it =
        std::max_element(chrom_names.begin(), chrom_names.end(),
                         [](const auto& c1, const auto& c2) { return c1.size() < c2.size(); });
    return it->size();
  }();

  modle::cooler::Cooler output_cool(c.path_to_output_matrix, cooler::Cooler::WRITE_ONLY,
                                    input_cool.get_bin_size(), max_chrom_name_size);

  for (const auto& [chrom_name, chrom_size] : input_cool.get_chroms()) {
    const auto t0 = absl::Now();
    fmt::print(stderr, FMT_STRING("Processing contacts for {}..."), chrom_name);
    const auto ncols =
        (chrom_size / bin_size) + static_cast<std::size_t>(chrom_size % bin_size != 0);
    const auto nrows =
        std::min(ncols, (c.diagonal_width / bin_size) +
                            static_cast<std::size_t>(c.diagonal_width % bin_size != 0));

    std::thread t([&, chrom_name = chrom_name]() {
      input_cool.stream_contacts_for_chrom(std::ref(pixel_queue), chrom_name, c.diagonal_width,
                                           bin_size);
      while (!pixel_queue.try_enqueue(END_OF_PIXEL_QUEUE)) {
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
      }
    });

    cmatrix.reset();
    cmatrix.resize(nrows, ncols);
    modle::cooler::Cooler::Pixel pixel;  // NOLINT
    const auto seed =
        c.seed + std::hash<std::string_view>{}(chrom_name) + std::hash<std::size_t>{}(ncols);

    auto rang_eng = PRNG(seed);
    auto gammad = std::gamma_distribution<double>{c.gamma_k, c.gamma_theta};
    auto uniformd = std::uniform_real_distribution<double>{static_cast<double>(bin_size) / -2.0,
                                                           static_cast<double>(bin_size) / 2.0};
    while (true) {
      pixel_queue.wait_dequeue(pixel);
      if (pixel == END_OF_PIXEL_QUEUE) {
        break;
      }
      assert(pixel.row <= pixel.col);  // NOLINT
      for (auto i = 0UL; i < pixel.count; ++i) {
        const auto pos1 =
            static_cast<double>(pixel.row * bin_size) + uniformd(rang_eng) + gammad(rang_eng);
        const auto pos2 =
            static_cast<double>(pixel.col * bin_size) + uniformd(rang_eng) - gammad(rang_eng);
        const auto bin1 =
            std::clamp(static_cast<std::size_t>(std::round(pos1 / static_cast<double>(bin_size))),
                       pixel.row, pixel.col);
        const auto bin2 =
            std::clamp(static_cast<std::size_t>(std::round(pos2 / static_cast<double>(bin_size))),
                       pixel.row, pixel.col);
        cmatrix.increment(bin1, bin2);
      }
    }
    t.join();

    output_cool.write_or_append_cmatrix_to_file(cmatrix, chrom_name, 0UL, chrom_size, chrom_size,
                                                true);
    fmt::print(stderr, FMT_STRING(" DONE in {}.\n"), absl::FormatDuration(absl::Now() - t0));
  }
}
}  // namespace modle::tools
