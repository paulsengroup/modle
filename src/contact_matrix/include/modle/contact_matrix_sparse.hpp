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

template <class N>
class ContactMatrixSerde;

template <class N = contacts_t>
class ContactMatrixSparse {
  static_assert(std::is_arithmetic_v<N>);

  friend class ContactMatrixSerde<N>;

 public:
  using value_type = N;
  template <class M>
  friend class ContactMatrixSparse;

  struct ChunkSize {
    usize size{0};
  };

 private:
  using BlockT = libcuckoo::cuckoohash_map<usize, N>;
  using SumT = typename std::conditional<std::is_floating_point_v<N>, double, i64>::type;
  static constexpr ChunkSize default_max_chunk_size{1'000'000};
  u64 _nrows{0};
  u64 _ncols{0};
  usize _cols_per_chunk{0};
  mutable std::vector<BlockT> _contact_blocks{};
  mutable std::atomic<SumT> _tot_contacts{0};
  mutable std::atomic<usize> _nnz{0};
  mutable std::atomic<bool> _global_stats_outdated{false};
  std::atomic<usize> _updates_missed{0};

 public:
  ContactMatrixSparse() = default;
#if defined(__clang__) && __clang_major__ < 9
  ContactMatrixSparse(ContactMatrixSparse<N>&& other) = default;
#else
  ContactMatrixSparse(ContactMatrixSparse<N>&& other) noexcept = default;
#endif
  inline ContactMatrixSparse(const ContactMatrixSparse<N>& other);
  inline ContactMatrixSparse(usize nrows, usize ncols,
                             ChunkSize max_chunk_size = default_max_chunk_size);
  inline ContactMatrixSparse(bp_t length, bp_t diagonal_width, bp_t bin_size,
                             ChunkSize max_chunk_size = default_max_chunk_size);
  ~ContactMatrixSparse() = default;

  // Operators
  inline ContactMatrixSparse<N>& operator=(const ContactMatrixSparse<N>& other);
#if defined(__clang__) && __clang_major__ < 9
  ContactMatrixSparse<N>& operator=(ContactMatrixSparse<N>&& other) = default;
#else
  ContactMatrixSparse<N>& operator=(ContactMatrixSparse<N>&& other) noexcept = default;
#endif

  // Thread-safe count getters and setters
  [[nodiscard]] inline N get(usize row, usize col) const;
  inline void set(usize row, usize col, N n);
  inline void add(usize row, usize col, N n);
  inline void subtract(usize row, usize col, N n);
  inline void increment(usize row, usize col);
  inline void decrement(usize row, usize col);

  [[nodiscard]] constexpr usize ncols() const;
  [[nodiscard]] constexpr usize nrows() const;
  [[nodiscard]] constexpr usize npixels() const;

  [[nodiscard]] constexpr usize get_n_of_missed_updates() const noexcept;
  [[nodiscard]] inline double get_fraction_of_missed_updates() const;
  [[nodiscard]] inline auto get_tot_contacts() const -> SumT;
  [[nodiscard]] inline usize get_nnz() const;
  [[nodiscard]] inline double get_avg_contact_density() const;

  [[nodiscard]] inline auto unsafe_get_tot_contacts() const -> SumT;
  [[nodiscard]] inline double unsafe_get_fraction_of_missed_updates() const;
  [[nodiscard]] inline usize unsafe_get_nnz() const;
  [[nodiscard]] inline double unsafe_get_avg_contact_density() const;

  [[nodiscard]] inline usize num_blocks() const noexcept;

 private:
  [[nodiscard]] static constexpr usize compute_cols_per_chunk(usize nrows,
                                                              ChunkSize max_chunk_size);
  [[nodiscard]] static constexpr usize compute_num_chunks(usize ncols, usize cols_per_chunk);
  [[nodiscard]] inline usize compute_block_idx(usize col) const noexcept;
  [[nodiscard]] inline auto get_block(usize col) noexcept -> BlockT&;
  [[nodiscard]] inline auto get_block(usize col) const noexcept -> const BlockT&;
  [[nodiscard]] inline auto lock_tables() const -> std::vector<typename BlockT::locked_table>;
  inline void update_global_stats(const std::vector<typename BlockT::locked_table>& tables) const;

  void bound_check_coords(usize row, usize col) const;
};

}  // namespace modle

#include "../../contact_matrix_sparse_impl.hpp"  // IWYU pragma: export
// IWYU pragma: "../../contact_matrix_sparse_impl.hpp"
