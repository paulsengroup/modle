// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <filesystem>
#include <libcuckoo/cuckoohash_map.hh>
#include <type_traits>
#include <vector>

#include "modle/common/common.hpp"

namespace modle {

template <class N = contacts_t>
class ContactMatrixSparse {
  static_assert(std::is_arithmetic_v<N>);

 public:
  using value_type = N;
  template <class M>
  friend class ContactMatrixSparse;

  struct ChunkSize {
    std::size_t size{0};
  };

 private:
  using BlockT = libcuckoo::cuckoohash_map<std::size_t, N>;
  using SumT = typename std::conditional<std::is_floating_point_v<N>, double, std::int64_t>::type;
  static constexpr ChunkSize default_max_chunk_size{1'000'000};
  std::uint64_t _nrows{0};
  std::uint64_t _ncols{0};
  std::size_t _cols_per_chunk{0};
  mutable std::vector<BlockT> _contact_blocks{};
  mutable std::atomic<SumT> _tot_contacts{0};
  mutable std::atomic<std::size_t> _nnz{0};
  mutable std::atomic<bool> _global_stats_outdated{false};
  std::atomic<std::size_t> _updates_missed{0};

 public:
  ContactMatrixSparse() = default;
  ContactMatrixSparse(ContactMatrixSparse<N>&& other) noexcept = default;
  inline ContactMatrixSparse(const ContactMatrixSparse<N>& other);
  inline ContactMatrixSparse(std::size_t nrows, std::size_t ncols,
                             ChunkSize max_chunk_size = default_max_chunk_size);
  inline ContactMatrixSparse(bp_t length, bp_t diagonal_width, bp_t bin_size,
                             ChunkSize max_chunk_size = default_max_chunk_size);
  ~ContactMatrixSparse() = default;

  // Operators
  inline ContactMatrixSparse<N>& operator=(const ContactMatrixSparse<N>& other);
  ContactMatrixSparse<N>& operator=(ContactMatrixSparse<N>&& other) noexcept = default;

  // Thread-safe count getters and setters
  [[nodiscard]] inline N get(std::size_t row, std::size_t col) const;
  inline void set(std::size_t row, std::size_t col, N n);
  inline void add(std::size_t row, std::size_t col, N n);
  inline void subtract(std::size_t row, std::size_t col, N n);
  inline void increment(std::size_t row, std::size_t col);
  inline void decrement(std::size_t row, std::size_t col);

  [[nodiscard]] constexpr std::size_t ncols() const;
  [[nodiscard]] constexpr std::size_t nrows() const;
  [[nodiscard]] constexpr std::size_t npixels() const;

  [[nodiscard]] constexpr std::size_t get_n_of_missed_updates() const noexcept;
  [[nodiscard]] inline double get_fraction_of_missed_updates() const;
  [[nodiscard]] inline auto get_tot_contacts() const -> SumT;
  [[nodiscard]] inline std::size_t get_nnz() const;
  [[nodiscard]] inline double get_avg_contact_density() const;

  [[nodiscard]] inline auto unsafe_get_tot_contacts() const -> SumT;
  [[nodiscard]] inline double unsafe_get_fraction_of_missed_updates() const;
  [[nodiscard]] inline std::size_t unsafe_get_nnz() const;
  [[nodiscard]] inline double unsafe_get_avg_contact_density() const;

  [[nodiscard]] inline std::size_t num_blocks() const noexcept;

 private:
  [[nodiscard]] static constexpr std::size_t compute_cols_per_chunk(std::size_t nrows,
                                                                    ChunkSize max_chunk_size);
  [[nodiscard]] static constexpr std::size_t compute_num_chunks(std::size_t ncols,
                                                                std::size_t cols_per_chunk);
  [[nodiscard]] inline std::size_t compute_block_idx(std::size_t col) const noexcept;
  [[nodiscard]] inline auto get_block(std::size_t col) noexcept -> BlockT&;
  [[nodiscard]] inline auto get_block(std::size_t col) const noexcept -> const BlockT&;
  [[nodiscard]] inline auto lock_tables() const -> std::vector<typename BlockT::locked_table>;
  inline void update_global_stats(const std::vector<typename BlockT::locked_table>& tables) const;

  void bound_check_coords(std::size_t row, std::size_t col) const;
};

}  // namespace modle

#include "../../contact_matrix_sparse_impl.hpp"  // IWYU pragma: export
// IWYU pragma: "../../contact_matrix_sparse_impl.hpp"
