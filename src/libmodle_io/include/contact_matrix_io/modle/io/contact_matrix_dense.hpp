// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once
#include <filesystem>
#include <hictk/cooler.hpp>
#include <string_view>
#include <vector>

#include "modle/common/common.hpp"
#include "modle/contact_matrix_dense.hpp"

namespace modle::io {

template <class N, usize chunk_size = 64U << 10U>
void append_contact_matrix_to_cooler(hictk::cooler::File& f, std::string_view chrom_name,
                                     const ContactMatrixDense<N>& matrix, bp_t offset = 0);

template <class N, usize chunk_size = 64U << 10U>
void append_contact_matrix_to_cooler(hictk::cooler::File& f, usize chrom_id,
                                     const ContactMatrixDense<N>& matrix, bp_t offset = 0);

namespace internal {
template <class N, class PixelT, usize chunk_size>
void append_contact_matrix_to_cooler(hictk::cooler::File& f, usize chrom_id,
                                     const ContactMatrixDense<N>& matrix, std::vector<PixelT>& buff,
                                     bp_t offset = 0);
}  // namespace internal

template <class N, class ChromIt>
[[nodiscard]] hictk::cooler::File init_cooler_file(const std::filesystem::path& path,
                                                   bool force_overwrite, ChromIt first_chrom,
                                                   ChromIt last_chrom, usize bin_size,
                                                   std::string_view assembly,
                                                   std::string_view generated_by,
                                                   std::string_view metadata_str);
template <class N, class ChromIt>
[[nodiscard]] hictk::cooler::File init_cooler_file(const std::filesystem::path& path,
                                                   bool force_overwrite, ChromIt first_chrom,
                                                   ChromIt last_chrom,
                                                   hictk::cooler::StandardAttributes attrs);

template <class N, class ChromNameIt, class ChromSizeIt>
[[nodiscard]] hictk::cooler::File init_cooler_file(const std::filesystem::path& path,
                                                   bool force_overwrite, ChromNameIt first_name,
                                                   ChromNameIt last_name, ChromSizeIt first_size,
                                                   usize bin_size, std::string_view assembly,
                                                   std::string_view generated_by,
                                                   std::string_view metadata_str);

template <class N, class ChromNameIt, class ChromSizeIt>
[[nodiscard]] hictk::cooler::File init_cooler_file(const std::filesystem::path& path,
                                                   bool force_overwrite, ChromNameIt first_name,
                                                   ChromNameIt last_name, ChromSizeIt first_size,
                                                   hictk::cooler::StandardAttributes attrs);

template <class N>
[[nodiscard]] hictk::cooler::File init_cooler_file(const std::filesystem::path& path,
                                                   bool force_overwrite, hictk::Reference chroms,
                                                   hictk::cooler::StandardAttributes attrs);

template <class N>
[[nodiscard]] ContactMatrixDense<N> read_contact_matrix_from_cooler(
    const hictk::cooler::File& f, std::string_view chrom_name, bp_t start_pos = 0,
    bp_t end_pos = std::numeric_limits<bp_t>::max(), bp_t diagonal_width = 0);
template <class N>
[[nodiscard]] ContactMatrixDense<N> read_contact_matrix_from_cooler(
    const hictk::cooler::File& f, usize chrom_id, bp_t start_pos = 0,
    bp_t end_pos = std::numeric_limits<bp_t>::max(), bp_t diagonal_width = 0);

template <class N>
[[nodiscard]] ContactMatrixDense<N> read_contact_matrix_from_cooler(
    const std::filesystem::path& path_to_cooler, std::string_view chrom_name, bp_t start_pos = 0,
    bp_t end_pos = std::numeric_limits<bp_t>::max(), bp_t diagonal_width = 0);
template <class N>
[[nodiscard]] ContactMatrixDense<N> read_contact_matrix_from_cooler(
    const std::filesystem::path& path_to_cooler, usize chrom_id, bp_t start_pos = 0,
    bp_t end_pos = std::numeric_limits<bp_t>::max(), bp_t diagonal_width = 0);

[[nodiscard]] bool query_returns_no_pixels(const hictk::cooler::File& f,
                                           std::string_view chrom_name, bp_t start_pos,
                                           bp_t end_pos);

}  // namespace modle::io

#include "../../../../contact_matrix_dense_io_impl.hpp"
