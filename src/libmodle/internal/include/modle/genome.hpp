// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <absl/types/span.h>
#include <fmt/format.h>
#include <parallel_hashmap/btree.h>
#include <xxhash.h>

#include <filesystem>
#include <iterator>
#include <limits>
#include <memory>
#include <optional>
#include <shared_mutex>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

#include "modle/bed/bed.hpp"
#include "modle/common/common.hpp"
#include "modle/common/utils.hpp"
#include "modle/contact_matrix_dense.hpp"
#include "modle/extrusion_barriers.hpp"
#include "modle/interval_tree.hpp"

namespace modle {

struct ExtrusionBarrier;

namespace internal {
class ContactMatrixLazy {
 public:
  using ContactMatrix = ContactMatrixDense<contacts_t>;

 private:
  mutable std::optional<ContactMatrix> _matrix{};
  mutable std::once_flag _alloc_flag{};
  std::once_flag _dealloc_flag{};
  std::uint64_t _nrows{0};
  std::uint64_t _ncols{0};

 public:
  ContactMatrixLazy() = default;
  ContactMatrixLazy(bp_t length, bp_t diagonal_width, bp_t bin_size) noexcept;
  explicit ContactMatrixLazy(ContactMatrix matrix) noexcept;
  ContactMatrixLazy(const ContactMatrixLazy& other) = delete;
  ContactMatrixLazy(ContactMatrixLazy&& other) noexcept;

  ~ContactMatrixLazy() = default;

  ContactMatrixLazy& operator=(const ContactMatrixLazy& other) = delete;
  ContactMatrixLazy& operator=(ContactMatrixLazy&& other) noexcept;

  [[nodiscard]] auto operator()() const noexcept -> const ContactMatrix&;
  [[nodiscard]] auto operator()() noexcept -> ContactMatrix&;
  [[nodiscard]] constexpr explicit operator bool() const noexcept;

  [[nodiscard]] constexpr std::uint64_t nrows() const noexcept;
  [[nodiscard]] constexpr std::uint64_t ncols() const noexcept;
  [[nodiscard]] constexpr std::uint64_t npixels() const noexcept;

  void deallocate() noexcept;
};

class Occupancy1DLazy {
  using BufferT = std::vector<std::atomic<std::uint64_t>>;
  mutable std::optional<BufferT> _buff{};
  mutable std::once_flag _alloc_flag{};
  std::once_flag _dealloc_flag{};
  std::size_t _size{0};

 public:
  Occupancy1DLazy() = default;
  Occupancy1DLazy(bp_t length, bp_t bin_size) noexcept;

  explicit Occupancy1DLazy(BufferT buff) noexcept;
  Occupancy1DLazy(const Occupancy1DLazy& other) = delete;
  Occupancy1DLazy(Occupancy1DLazy&& other) noexcept;

  ~Occupancy1DLazy() = default;

  Occupancy1DLazy& operator=(const Occupancy1DLazy& other) = delete;
  Occupancy1DLazy& operator=(Occupancy1DLazy&& other) noexcept;

  [[nodiscard]] const std::vector<std::atomic<std::uint64_t>>& operator()() const noexcept;
  [[nodiscard]] std::vector<std::atomic<std::uint64_t>>& operator()() noexcept;
  [[nodiscard]] explicit operator bool() const noexcept;

  [[nodiscard]] constexpr std::size_t size() const noexcept;

  void deallocate() noexcept;
};
}  // namespace internal

class Chromosome {
  std::string _name{};
  std::size_t _id{(std::numeric_limits<std::size_t>::max)()};
  bp_t _size{};

 public:
  Chromosome() = default;
  Chromosome(std::size_t id, std::string name, bp_t size) noexcept;
  [[nodiscard]] constexpr explicit operator bool() const noexcept;

  [[nodiscard]] constexpr std::size_t id() const noexcept;
  [[nodiscard]] std::string_view name() const noexcept;
  [[nodiscard]] const char* name_cstr() const noexcept;
  [[nodiscard]] constexpr bp_t size() const noexcept;

  [[nodiscard]] std::uint64_t hash(XXH3_state_t& state) const;
  [[nodiscard]] std::uint64_t hash(XXH3_state_t& state, std::uint64_t seed) const;

  [[nodiscard]] constexpr bool operator==(const Chromosome& other) const noexcept;
  [[nodiscard]] constexpr bool operator!=(const Chromosome& other) const noexcept;
  [[nodiscard]] constexpr bool operator<(const Chromosome& other) const noexcept;
  [[nodiscard]] constexpr bool operator<=(const Chromosome& other) const noexcept;
  [[nodiscard]] constexpr bool operator>(const Chromosome& other) const noexcept;
  [[nodiscard]] constexpr bool operator>=(const Chromosome& other) const noexcept;
};

class GenomicInterval {
 public:
  using ContactMatrix = internal::ContactMatrixLazy::ContactMatrix;
  friend class Genome;

 private:
  std::size_t _id{(std::numeric_limits<std::size_t>::max)()};
  std::shared_ptr<const Chromosome> _chrom{};
  bp_t _start{};
  bp_t _end{};

  std::vector<ExtrusionBarrier> _barriers{};

  internal::ContactMatrixLazy _contacts{};
  internal::Occupancy1DLazy _lef_1d_occupancy{};

 public:
  GenomicInterval() = default;
  GenomicInterval(std::size_t id, const std::shared_ptr<const Chromosome>& chrom,
                  bp_t contact_matrix_resolution, bp_t diagonal_width);
  GenomicInterval(std::size_t id, std::shared_ptr<const Chromosome> chrom, bp_t start, bp_t end,
                  bp_t contact_matrix_resolution, bp_t diagonal_width);
  template <typename It>
  GenomicInterval(std::size_t id, const std::shared_ptr<const Chromosome>& chrom,
                  bp_t contact_matrix_resolution, bp_t diagonal_width, It first_barrier,
                  It last_barrier);
  template <typename It>
  GenomicInterval(std::size_t id, std::shared_ptr<const Chromosome> chrom, bp_t start, bp_t end,
                  bp_t contact_matrix_resolution, bp_t diagonal_width, It first_barrier,
                  It last_barrier);

  ~GenomicInterval() = default;

  GenomicInterval(const GenomicInterval& other) = delete;
  GenomicInterval(GenomicInterval&& other) noexcept = default;

  GenomicInterval& operator=(const GenomicInterval& other) = delete;
  GenomicInterval& operator=(GenomicInterval&& other) noexcept = default;

  [[nodiscard]] constexpr bool operator==(const GenomicInterval& other) const noexcept;
  [[nodiscard]] constexpr bool operator!=(const GenomicInterval& other) const noexcept;
  [[nodiscard]] constexpr bool operator<(const GenomicInterval& other) const noexcept;
  [[nodiscard]] constexpr bool operator<=(const GenomicInterval& other) const noexcept;
  [[nodiscard]] constexpr bool operator>(const GenomicInterval& other) const noexcept;
  [[nodiscard]] constexpr bool operator>=(const GenomicInterval& other) const noexcept;

  void add_extrusion_barrier(const bed::BED& record, double default_barrier_stp_active,
                             double default_barrier_stp_inactive);
  void add_extrusion_barriers(std::vector<ExtrusionBarrier> barriers);

  [[nodiscard]] constexpr std::size_t id() const noexcept;
  [[nodiscard]] std::uint64_t hash(XXH3_state_t& state) const;
  [[nodiscard]] std::uint64_t hash(XXH3_state_t& state, std::uint64_t seed) const;
  [[nodiscard]] const Chromosome& chrom() const noexcept;
  [[nodiscard]] constexpr bp_t start() const noexcept;
  [[nodiscard]] constexpr bp_t end() const noexcept;
  [[nodiscard]] constexpr bp_t size() const noexcept;
  [[nodiscard]] constexpr std::uint64_t npixels() const noexcept;
  [[nodiscard]] std::size_t num_barriers() const;
  [[nodiscard]] auto barriers() const noexcept -> const std::vector<ExtrusionBarrier>&;
  [[nodiscard]] auto barriers() noexcept -> std::vector<ExtrusionBarrier>&;
  [[nodiscard]] auto contacts() const noexcept -> const ContactMatrix&;
  [[nodiscard]] auto contacts() noexcept -> ContactMatrix&;
  [[nodiscard]] auto lef_1d_occupancy() const noexcept
      -> const std::vector<std::atomic<std::uint64_t>>&;
  [[nodiscard]] auto lef_1d_occupancy() noexcept -> std::vector<std::atomic<std::uint64_t>>&;
  void deallocate() noexcept;

  template <typename H>
  inline friend H AbslHashValue(H h, const GenomicInterval& c);
};

class Genome {
  std::vector<std::shared_ptr<const Chromosome>> _chroms{};
  phmap::btree_set<GenomicInterval> _intervals{};

  std::size_t _size{};
  std::size_t _simulated_size{};
  std::size_t _num_barriers{};

 public:
  Genome() = default;
  Genome(const std::filesystem::path& path_to_chrom_sizes,
         const std::filesystem::path& path_to_extr_barriers,
         const std::filesystem::path& path_to_genomic_intervals, bp_t contact_matrix_resolution,
         bp_t contact_matrix_diagonal_witdh, double default_barrier_pbb, double default_barrier_puu,
         bool interpret_name_field_as_puu);

  using iterator = decltype(_intervals)::iterator;
  using const_iterator = decltype(_intervals)::const_iterator;

  [[nodiscard]] auto begin() -> iterator;
  [[nodiscard]] auto end() -> iterator;

  [[nodiscard]] auto begin() const -> const_iterator;
  [[nodiscard]] auto end() const -> const_iterator;

  [[nodiscard]] auto cbegin() const -> const_iterator;
  [[nodiscard]] auto cend() const -> const_iterator;

  // NOTE: The find and contains overloads taking Chromosome or string_view as query can be quite
  // slow, as they are performing a linear search through the genome
  [[nodiscard]] auto find(const GenomicInterval& query) -> iterator;
  [[nodiscard]] auto find(const GenomicInterval& query) const -> const_iterator;
  [[nodiscard]] auto find(std::size_t query) -> iterator;
  [[nodiscard]] auto find(std::size_t query) const -> const_iterator;
  [[nodiscard]] auto find(const Chromosome& query) -> iterator;
  [[nodiscard]] auto find(const Chromosome& query) const -> const_iterator;
  [[nodiscard]] auto find(std::string_view query) -> iterator;
  [[nodiscard]] auto find(std::string_view query) const -> const_iterator;

  [[nodiscard]] bool contains(const GenomicInterval& query) const;
  [[nodiscard]] bool contains(std::size_t query) const;
  [[nodiscard]] bool contains(const Chromosome& query) const;
  [[nodiscard]] bool contains(std::string_view query) const;

  constexpr const std::vector<std::shared_ptr<const Chromosome>>& chromosomes() const noexcept;

  [[nodiscard]] std::size_t num_intervals() const noexcept;
  [[nodiscard]] constexpr std::size_t size() const noexcept;
  [[nodiscard]] constexpr std::size_t simulated_size() const noexcept;
  [[nodiscard]] std::size_t num_chromosomes() const noexcept;
  [[nodiscard]] constexpr std::size_t num_barriers() const noexcept;

  [[nodiscard]] const Chromosome& chromosome_with_longest_name() const noexcept;
  [[nodiscard]] const Chromosome& longest_chromosome() const noexcept;
  [[nodiscard]] const GenomicInterval& longest_interval() const noexcept;
  [[nodiscard]] const GenomicInterval& interval_with_most_barriers() const noexcept;
  [[nodiscard]] std::size_t max_target_contacts(std::size_t bin_size, std::size_t diagonal_width,
                                                double target_contact_density,
                                                std::size_t simulation_iterations,
                                                double lef_fraction_contact_sampling,
                                                double nlefs_per_mbp, std::size_t ncells) const;

 private:
  /// Import chromosomes from a chrom.sizes file
  [[nodiscard]] static std::vector<std::shared_ptr<const Chromosome>> import_chromosomes(
      const std::filesystem::path& path_to_chrom_sizes);

  /// Import genomic intervals from a BED file. If path to BED file is empty, assume entire
  /// chromosomes are to be simulated.
  phmap::btree_set<GenomicInterval> import_genomic_intervals(
      const std::filesystem::path& path_to_bed,
      const std::vector<std::shared_ptr<const Chromosome>>& chromosomes,
      bp_t contact_matrix_resolution, bp_t diagonal_width);

  /// Parse a BED file containing the genomic coordinates of extrusion barriers and add them to
  /// the Genome
  static std::size_t map_barriers_to_intervals(phmap::btree_set<GenomicInterval>& intervals,
                                               const bed::BED_tree<>& barriers_bed,
                                               double default_barrier_pbb,
                                               double default_barrier_puu,
                                               bool interpret_name_field_as_puu);
};

}  // namespace modle

template <>
struct fmt::formatter<modle::Chromosome> {
  constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin());
  template <typename FormatContext>
  auto format(const modle::Chromosome& chrom, FormatContext& ctx) const -> decltype(ctx.out());
};

template <>
struct fmt::formatter<modle::GenomicInterval> {
  constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin());
  template <typename FormatContext>
  auto format(const modle::GenomicInterval& gi, FormatContext& ctx) const -> decltype(ctx.out());
};

#include "../../genome_impl.hpp"  // IWYU pragma: export
