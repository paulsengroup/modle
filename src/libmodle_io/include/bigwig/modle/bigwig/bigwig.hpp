// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// IWYU pragma: no_include "modle/src/libio/bigwig_impl.hpp"

#include <absl/types/span.h>  // for Span
#include <fmt/std.h>
#include <parallel_hashmap/btree.h>

#include <filesystem>   // for path
#include <string_view>  // for string_view
#include <type_traits>  // for enable_if

#include "modle/common/common.hpp"  // for u32, usize, u64

// clang-format off
#include "modle/common/suppress_compiler_warnings.hpp"  // for DISABLE_WARNING_PADDED, DISABLE_W...

// clang-format on
DISABLE_WARNING_PUSH
DISABLE_WARNING_PADDED
#include "libbigwig/bigWig.h"  // for bigWigFile_t
DISABLE_WARNING_POP

namespace modle::io::bigwig {

class Reader {
 public:
  // TODO replace with flat_hash after issues with ASAN and abseil have been addressed
  // using Chromosomes = phmap::flat_hash_map<std::string, bp_t>;
  using Chromosomes = phmap::btree_map<std::string, bp_t>;
  using Intervals =
      std::unique_ptr<bwOverlappingIntervals_t, decltype(&bwDestroyOverlappingIntervals)>;
  using StatsType = bwStatsType;

 private:
  std::filesystem::path _fname{};
  bigWigFile_t* _fp{nullptr};
  Chromosomes _chroms{};

 public:
  Reader() = default;
  Reader(const Reader&) = delete;
  Reader(Reader&& other) noexcept;
  explicit Reader(std::filesystem::path name);
  ~Reader();

  Reader& operator=(const Reader&) = delete;
  Reader& operator=(Reader&& other) noexcept;

  [[nodiscard]] constexpr explicit operator bool() const noexcept;
  [[nodiscard]] constexpr const std::filesystem::path& path() const noexcept;
  [[nodiscard]] constexpr auto chromosomes() const noexcept -> const Chromosomes&;

  [[nodiscard]] std::vector<float> read_values(const std::string& chrom, bp_t start, bp_t end);
  [[nodiscard]] auto get_intervals(const std::string& chrom, bp_t start, bp_t end) -> Intervals;
  template <StatsType stat>
  [[nodiscard]] double stats(const std::string& chrom, bp_t start, bp_t end);
  template <StatsType stat>
  [[nodiscard]] std::vector<double> stats(const std::string& chrom, bp_t start, bp_t end,
                                          bp_t window_size);
  template <StatsType stat>
  void stats(const std::string& chrom, bp_t start, bp_t end, bp_t window_size,
             std::vector<double>& buff);

 private:
  void validate_query(const std::string& chrom, bp_t start, bp_t end);
};

class Writer {
  static constexpr uint_fast8_t DEFAULT_ZOOM_LEVELS{10};

  std::filesystem::path _fname{};
  bigWigFile_t* _fp{nullptr};
  uint_fast8_t _zoom_levels{DEFAULT_ZOOM_LEVELS};
  bool _initialized{false};

 public:
  Writer() = default;
  Writer(const Writer&) = delete;
  Writer(Writer&& other) noexcept;
  explicit Writer(std::filesystem::path name, uint_fast8_t zoom_levels = DEFAULT_ZOOM_LEVELS);
  ~Writer();

  Writer& operator=(const Writer&) = delete;
  Writer& operator=(Writer&& other) noexcept;

  [[nodiscard]] constexpr explicit operator bool() const noexcept;
  [[nodiscard]] constexpr const std::filesystem::path& path() const noexcept;

  template <class Chromosomes>
  void write_chromosomes(Chromosomes& chroms);
  void write_chromosomes(const std::vector<std::string>& chrom_names,
                         const std::vector<u32>& chrom_sizes);
  void write_chromosomes(const char* const* chrom_names, const u32* chrom_sizes, usize num_chroms);

  template <class N, class = std::enable_if<std::is_arithmetic_v<N>>>
  void write_range(std::string_view chrom_name, absl::Span<N> values, u64 span, u64 step,
                   u64 offset = 0);
};

}  // namespace modle::io::bigwig

#include "../../../../bigwig_impl.hpp"  // IWYU pragma: export
