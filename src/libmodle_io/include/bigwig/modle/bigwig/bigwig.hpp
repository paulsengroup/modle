// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// IWYU pragma: no_include "modle/src/libio/bigwig_impl.hpp"

#include <absl/types/span.h>  // for Span

#include <filesystem>   // for path
#include <mutex>        // for mutex
#include <string_view>  // for string_view
#include <type_traits>  // for enable_if

#include "modle/common/common.hpp"  // for u32, usize, u64

// clang-format off
#include "modle/common/suppress_compiler_warnings.hpp"  // for DISABLE_WARNING_PADDED, DISABLE_W...

// clang-format on
DISABLE_WARNING_PUSH
DISABLE_WARNING_PADDED
#include "libBigWig/bigWig.h"  // for bigWigFile_t
DISABLE_WARNING_POP

namespace modle::io::bigwig {

class Writer {
  static inline u32 _global_bigwig_files_opened{0};
  // This mutex protects read and write access to _global_bigwig_files_opened
  static inline std::mutex _global_state_mutex;

  static constexpr uint_fast8_t DEFAULT_ZOOM_LEVELS{10};
  static constexpr usize DEFAULT_BUFFER_SIZE{1U << 17U};  // 128 KiB

  std::filesystem::path _fname{};
  bigWigFile_t* _fp{nullptr};
  uint_fast8_t _zoom_levels{DEFAULT_ZOOM_LEVELS};
  usize _buff_size{DEFAULT_BUFFER_SIZE};
  bool _initialized{false};
  u32 _offset{0};

 public:
  Writer() = default;
  Writer(const Writer&) = delete;
  Writer(Writer&& other) noexcept;
  explicit Writer(std::filesystem::path name, uint_fast8_t zoom_levels = DEFAULT_ZOOM_LEVELS,
                  usize buff_size = DEFAULT_BUFFER_SIZE);
  ~Writer();

  Writer& operator=(const Writer&) = delete;
  Writer& operator=(Writer&& other) noexcept;

  template <class Chromosomes>
  inline void write_chromosomes(const Chromosomes& chroms);
  template <class Str>
  inline void write_chromosomes(absl::Span<Str> chrom_names, absl::Span<const u32> chrom_sizes);
  void write_chromosomes(const char* const* chrom_names, const u32* chrom_sizes, usize num_chroms);

  template <class N, class = std::enable_if<std::is_arithmetic_v<N>>>
  inline void write_range(std::string_view chrom_name, absl::Span<N> values, u64 span, u64 step);

  [[nodiscard]] const std::filesystem::path& path() const noexcept;
};

}  // namespace modle::io::bigwig

#include "../../../../bigwig_impl.hpp"  // IWYU pragma: export
