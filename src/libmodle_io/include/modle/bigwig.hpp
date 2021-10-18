// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// IWYU pragma: no_include "modle/src/libio/bigwig_impl.hpp"

#include <absl/container/flat_hash_map.h>  // for flat_hash_map

#include <iterator>     // for pair
#include <memory>       // for unique_ptr
#include <string>       // for string
#include <string_view>  // for string_view
#include <vector>       // for vector

#include "modle/common/common.hpp"  // for u64, i32, u32, uint...

// clang-format off
#include "modle/common/suppress_compiler_warnings.hpp"  // for DISABLE_WARNING_PADDED, DISABLE_W...

// clang-format on
DISABLE_WARNING_PUSH
DISABLE_WARNING_PADDED
#include "libBigWig/bigWig.h"  // for bigWigFile_t

DISABLE_WARNING_POP

namespace modle::bigwig {

enum : std::uint_fast8_t { DEFAULT_ZOOM_LEVELS = 10 };
inline constexpr usize DEFAULT_BUFF_SIZE{1U << 17U};  // 128 KiB

void close_bigwig_file(bigWigFile_t* fp);
using file = std::unique_ptr<bigWigFile_t, decltype(&bigwig::close_bigwig_file)>;

template <typename N1, typename N2>
inline u64 write_range(absl::flat_hash_map<std::pair<std::string, N1>, std::vector<N2>>& data,
                       u64 offset, u64 span, u64 step, std::string_view output_file);

template <typename N1, typename N2>
inline void write_range(absl::flat_hash_map<std::pair<std::string, N1>, std::vector<N2>>& data,
                        u64 offset, u64 span, u64 step, bigwig::file& bigwig_fp);

void write_range(const std::string& chrom_name, const std::vector<double>& vals, u64 offset,
                 u64 span, u64 step, bigwig::file& bigwig_fp);

template <typename I>  // This overload is not as efficient as the others. This shouldn't be an
                       // issue, unless the caller is initializing the file with hundreds of
                       // thousands of chromosomes
[[nodiscard]] bigwig::file init_bigwig_file(std::string_view output_path,
                                            std::vector<std::pair<std::string, I>>& chromosomes,
                                            i32 zoom_levels = DEFAULT_ZOOM_LEVELS,
                                            usize buff_size = DEFAULT_BUFF_SIZE);

[[nodiscard]] bigwig::file init_bigwig_file(std::string_view output_path,
                                            std::vector<char*>& chrom_names,
                                            std::vector<u32>& chrom_sizes,
                                            i32 zoom_levels = DEFAULT_ZOOM_LEVELS,
                                            usize buff_size = DEFAULT_BUFF_SIZE);

[[nodiscard]] bigwig::file init_bigwig_file(std::string_view output_path, std::string& chrom_name,
                                            u64 chrom_size, i32 zoom_levels = DEFAULT_ZOOM_LEVELS,
                                            usize buff_size = DEFAULT_BUFF_SIZE);

void close_bigwig_file(bigwig::file fp);
}  // namespace modle::bigwig

#include "../../bigwig_impl.hpp"  // IWYU pragma: export
