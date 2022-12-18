#pragma once
#include <coolerpp/coolerpp.hpp>
#include <string_view>
#include <vector>

#include "modle/common/common.hpp"
#include "modle/contact_matrix_dense.hpp"

namespace modle::io {
template <class N, usize chunk_size = 64U << 10U>
void append_contact_matrix_to_cooler(coolerpp::File& f, std::string_view chrom_name,
                                     const ContactMatrixDense<N>& matrix);

template <class N, usize chunk_size = 64U << 10U>
void append_contact_matrix_to_cooler(coolerpp::File& f, usize chrom_id,
                                     const ContactMatrixDense<N>& matrix);

namespace internal {
template <class N, class PixelT, usize chunk_size>
void append_contact_matrix_to_cooler(coolerpp::File& f, usize chrom_id,
                                     const ContactMatrixDense<N>& matrix,
                                     std::vector<PixelT>& buff);
}  // namespace internal

template <class N>
[[nodiscard]] ContactMatrixDense<N> read_contact_matrix_from_cooler(
    const coolerpp::File& f, std::string_view chrom_name, bp_t start_pos = 0,
    bp_t end_pos = std::numeric_limits<bp_t>::max(), bp_t diagonal_width = 0);
template <class N>
[[nodiscard]] ContactMatrixDense<N> read_contact_matrix_from_cooler(
    const coolerpp::File& f, usize chrom_id, bp_t start_pos = 0,
    bp_t end_pos = std::numeric_limits<bp_t>::max(), bp_t diagonal_width = 0);

template <class N>
[[nodiscard]] ContactMatrixDense<N> read_contact_matrix_from_cooler(
    const std::filesystem::path& path_to_cooler, std::string_view chrom_name, bp_t start_pos = 0,
    bp_t end_pos = std::numeric_limits<bp_t>::max(), bp_t diagonal_width = 0);
template <class N>
[[nodiscard]] ContactMatrixDense<N> read_contact_matrix_from_cooler(
    const std::filesystem::path& path_to_cooler, usize chrom_id, bp_t start_pos = 0,
    bp_t end_pos = std::numeric_limits<bp_t>::max(), bp_t diagonal_width = 0);

[[nodiscard]] bool query_returns_no_pixels(const coolerpp::File& f, std::string_view chrom_name,
                                           bp_t start_pos, bp_t end_pos);

}  // namespace modle::io

#include "../../../../contact_matrix_dense_io_impl.hpp"
