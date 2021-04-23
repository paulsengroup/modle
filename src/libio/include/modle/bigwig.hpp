#pragma once

#include <absl/container/flat_hash_map.h>  // for flat_hash_map
#include <bigWig.h>                        // for bigWigFile_t

#include <cstddef>      // IWYU pragma: keep for size_t
#include <cstdint>      // for uint64_t, int32_t, uint32_t, uint_fast8_t
#include <iterator>     // for pair
#include <memory>       // for unique_ptr
#include <string>       // for string
#include <string_view>  // for string_view
#include <vector>       // for vector

namespace modle::bigwig {

enum : uint_fast8_t { DEFAULT_ZOOM_LEVELS = 10 };
inline constexpr size_t DEFAULT_BUFF_SIZE{1U << 17U};  // 128 KiB

inline void close_bigwig_file(bigWigFile_t* fp);
using file = std::unique_ptr<bigWigFile_t, decltype(&bigwig::close_bigwig_file)>;

template <typename N1, typename N2>
inline uint64_t write_range(absl::flat_hash_map<std::pair<std::string, N1>, std::vector<N2>>& data,
                            uint64_t offset, uint64_t span, uint64_t step,
                            std::string_view output_file);

template <typename N1, typename N2>
inline void write_range(absl::flat_hash_map<std::pair<std::string, N1>, std::vector<N2>>& data,
                        uint64_t offset, uint64_t span, uint64_t step, bigwig::file& bigwig_fp);

inline void write_range(const std::string& chrom_name, const std::vector<double>& vals,
                        uint64_t offset, uint64_t span, uint64_t step, bigwig::file& bigwig_fp);

template <typename I>  // This overload is not as efficient as the others. This shouldn't be an
                       // issue, unless the caller is initializing the file with hundreds of
                       // thousands of chromosomes
[[nodiscard]] inline bigwig::file init_bigwig_file(
    std::string_view output_path, std::vector<std::pair<std::string, I>>& chromosomes,
    int32_t zoom_levels = DEFAULT_ZOOM_LEVELS, size_t buff_size = DEFAULT_BUFF_SIZE);

[[nodiscard]] inline bigwig::file init_bigwig_file(std::string_view output_path,
                                                   std::vector<char*>& chrom_names,
                                                   std::vector<uint32_t>& chrom_sizes,
                                                   int32_t zoom_levels = DEFAULT_ZOOM_LEVELS,
                                                   size_t buff_size = DEFAULT_BUFF_SIZE);

[[nodiscard]] inline bigwig::file init_bigwig_file(std::string_view output_path,
                                                   std::string& chrom_name, uint64_t chrom_size,
                                                   int32_t zoom_levels = DEFAULT_ZOOM_LEVELS,
                                                   size_t buff_size = DEFAULT_BUFF_SIZE);

inline void close_bigwig_file(bigwig::file fp);
}  // namespace modle::bigwig

#include "../../bigwig_impl.hpp"  // IWYU pragma: keep
