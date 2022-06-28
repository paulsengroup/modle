// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <boost/serialization/access.hpp>
#include <istream>
#include <ostream>
#include <type_traits>
#include <vector>

#include "modle/common/common.hpp"
#include "modle/contact_matrix_dense.hpp"
#include "modle/contact_matrix_sparse.hpp"
#include "modle/interval_tree.hpp"

namespace modle {

template <class N>
class ContactMatrixSerde {
  static_assert(std::is_arithmetic_v<N>);
  using SumT = typename std::conditional<std::is_floating_point_v<N>, double, i64>::type;

  struct Header;
  class ChunkBuffer;
  struct ChunkMetadata;

  Header _header{};
  ChunkBuffer _chunk{};
  std::vector<usize> _sorting_buff{};
  u32 _compression_lvl{1};

 public:
  ContactMatrixSerde() = default;

  ContactMatrixSerde(const ContactMatrixSerde& other) = delete;
  ContactMatrixSerde(ContactMatrixSerde&& other) noexcept = default;

  ~ContactMatrixSerde() noexcept = default;

  ContactMatrixSerde& operator=(const ContactMatrixSerde& other) = delete;
  ContactMatrixSerde& operator=(ContactMatrixSerde&& other) noexcept = default;

  template <class ContactMatrix>
  inline void read(const std::filesystem::path& path, ContactMatrix& matrix);
  template <class ContactMatrix>
  inline void read(std::istream& stream, ContactMatrix& matrix);

  template <class ContactMatrix>
  inline void write(const std::filesystem::path& out_path, const ContactMatrix& matrix,
                    std::string_view metadata);
  template <class ContactMatrix>
  inline void write(std::ostream& stream, const ContactMatrix& matrix, std::string_view metadata);

  template <class ContactMatrix>
  inline void write_destructive(const std::filesystem::path& out_path, ContactMatrix& matrix,
                                std::string_view metadata);
  template <class ContactMatrix>
  inline void write_destructive(std::ostream& stream, ContactMatrix& matrix,
                                std::string_view metadata);

 private:
  inline void write_empty_header(std::ostream& stream);
  inline void read_header(std::istream& stream);
  inline void write_header(std::ostream& stream);

  inline void read_internal(std::istream& stream, ContactMatrixDense<N>& matrix);
  inline void read_internal(std::istream& stream, ContactMatrixSparse<N>& matrix);

  inline void write_internal(std::ostream& stream, const ContactMatrixDense<N>& matrix,
                             std::string_view metadata);
  inline void write_internal(std::ostream& stream, const ContactMatrixSparse<N>& matrix,
                             std::string_view metadata);
  inline void write_internal_destructive(std::ostream& stream, ContactMatrixDense<N>& matrix,
                                         std::string_view metadata);
  inline void write_internal_destructive(std::ostream& stream, ContactMatrixSparse<N>& matrix,
                                         std::string_view metadata);

  inline void read_chunk(std::istream& stream, i64 offset = -1);
  inline void write_chunk(std::ostream& stream, i64 offset = -1);
  inline void reset() noexcept;

  struct ChunkSpan {
    usize start;
    usize end;
  };

  struct Header {
    friend class boost::serialization::access;
    std::string metadata{};
    u64 nrows{0};
    u64 ncols{0};
    SumT tot_contacts{0};
    usize nnz{0};
    usize updates_missed{0};
    modle::IITree<usize, ChunkMetadata> chunk_metadata{};

   private:
    usize _num_chunks{0};

   public:
    Header() = default;
    template <class ContactMatrix>
    inline Header(const ContactMatrix& m, std::string_view metadata_str);
    [[nodiscard]] inline usize num_chunks() const noexcept;

    inline void clear() noexcept;

   private:
    template <class BoostArchive>
    [[maybe_unused]] inline void save(BoostArchive& ar, unsigned int version) const;
    template <class BoostArchive>
    [[maybe_unused]] inline void load(BoostArchive& ar, unsigned int version);
    template <class BoostArchive>
    [[maybe_unused]] inline void serialize(BoostArchive& ar, unsigned int version);

    static constexpr usize compute_num_cols_per_chunk(
        usize nrows, usize max_chunk_size_bytes = 64ULL * (1024ULL << 10U)) noexcept;
  };

  struct ChunkMetadata {
    i64 offset{-1};
    usize nnz{};
    SumT sum{};

    constexpr bool operator<(const ChunkMetadata& other) noexcept;
    constexpr bool operator==(const ChunkMetadata& other) noexcept;

    friend class boost::serialization::access;

   private:
    template <class BoostArchive>
    inline void serialize(BoostArchive& ar, unsigned int version);
  };

  class ChunkBuffer {
    std::vector<usize> _idx_buff{};
    std::vector<N> _count_buff{};

   public:
    static constexpr usize compute_size_from_bytes(usize bytes) noexcept;
    static constexpr usize DEFAULT_CAPACITY{compute_size_from_bytes(64ULL * (1024ULL << 10U))};

    friend class ContactMatrixSerde<N>;
    friend class boost::serialization::access;
    inline ChunkBuffer();

    [[nodiscard]] inline usize size() const noexcept;
    [[nodiscard]] inline usize capacity() const noexcept;
    [[nodiscard]] inline bool empty() const noexcept;
    [[nodiscard]] inline bool full() const noexcept;
    inline void reserve(usize new_size);

    inline void sort(std::vector<usize>& sorting_buff);
    inline auto sum_contacts() const noexcept -> SumT;

    inline void clear() noexcept;
    inline void push_back(usize i, N count) noexcept;

    inline i64 serialize(std::ostream& stream, u32 compression_level = 1) const;
    static inline i64 deserialize(std::istream& stream, ChunkBuffer& buff);

    template <class BoostArchive>
    inline void serialize(BoostArchive& ar, unsigned int version);
  };
};
}  // namespace modle

#include "../../contact_matrix_serde_impl.hpp"
