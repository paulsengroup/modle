// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <absl/container/btree_set.h>  // for btree_set, btree_set_container<>::const_iterator
#include <absl/types/span.h>           // for Span
#include <xxh3.h>                      // for XXH3_state_t, XXH_INLINE_XXH3_state_t

#include <boost/filesystem/path.hpp>  // for path
#include <iterator>                   // for iterator_traits
#include <limits>                     // for numeric_limits
#include <memory>                     // for shared_ptr
#include <shared_mutex>               // for shared_mutex
#include <string>                     // for string
#include <string_view>                // for string_view
#include <type_traits>                // for enable_if_t, remove_cv_t
#include <vector>                     // for vector

#include "modle/bed/bed.hpp"             // for BED (ptr only), BED_tree, BED_tree<>::value_type
#include "modle/common/common.hpp"       // for bp_t, contacts_t, u64, u32, u8
#include "modle/common/utils.hpp"        // for ndebug_defined
#include "modle/contacts.hpp"            // for ContactMatrix
#include "modle/extrusion_barriers.hpp"  // for ExtrusionBarrier
#include "modle/interval_tree.hpp"       // for IITree, IITree::IITree<I, T>

namespace modle {

class ExtrusionBarrier;

class Chromosome {
  using contact_matrix_t = ContactMatrix<contacts_t>;
  using bed_tree_value_t = bed::BED_tree<>::value_type;
  friend class Genome;

 public:
  Chromosome() = default;
  Chromosome(usize id, std::string_view chrom_name, bp_t chrom_start, bp_t chrom_end,
             bp_t chrom_size, bool ok_ = true);

  Chromosome(usize id, const bed::BED& chrom, bool ok_ = true);
  Chromosome(usize id, const bed::BED& chrom, const IITree<bp_t, ExtrusionBarrier>& barriers,
             bool ok_ = true);
  Chromosome(usize id, const bed::BED& chrom, IITree<bp_t, ExtrusionBarrier>&& barriers,
             bool ok_ = true);

  Chromosome(usize id, std::string_view chrom_name, bp_t chrom_start, bp_t chrom_end,
             bp_t chrom_size, const IITree<bp_t, ExtrusionBarrier>& barriers, bool ok_ = true);
  Chromosome(usize id, std::string_view chrom_name, bp_t chrom_start, bp_t chrom_end,
             bp_t chrom_size, IITree<bp_t, ExtrusionBarrier>&& barriers, bool ok_ = true);

  ~Chromosome() = default;

  Chromosome(const Chromosome& other);
  Chromosome(Chromosome&& other) noexcept;

  Chromosome& operator=(const Chromosome& other);
  Chromosome& operator=(Chromosome&& other) noexcept;

  [[nodiscard]] bool operator==(const Chromosome& other) const noexcept(utils::ndebug_defined());
  [[nodiscard]] bool operator==(std::string_view other_name) const noexcept;
  [[nodiscard]] bool operator<(const Chromosome& other) const noexcept(utils::ndebug_defined());

  void add_extrusion_barrier(const bed::BED& record, double ctcf_prob_occ_to_occ,
                             double ctcf_prob_nocc_to_nocc);

  template <typename Iter,
            typename = std::enable_if_t<std::is_same_v<
                std::remove_cv_t<typename std::iterator_traits<Iter>::value_type>, bed::BED>>>
  inline void add_extrusion_barrier(Iter barriers_begin, Iter barriers_end);

  [[nodiscard]] usize id() const;
  [[nodiscard]] std::string_view name() const;
  [[nodiscard]] const char* name_cstr() const;
  [[nodiscard]] constexpr bp_t start_pos() const;
  [[nodiscard]] constexpr bp_t end_pos() const;
  [[nodiscard]] constexpr bp_t size() const;
  [[nodiscard]] constexpr bp_t simulated_size() const;
  [[nodiscard]] usize npixels() const;
  [[nodiscard]] constexpr usize npixels(bp_t diagonal_width, bp_t bin_size) const;
  [[nodiscard]] bool ok() const;
  [[nodiscard]] usize num_lefs(double nlefs_per_mbp) const;
  [[nodiscard]] usize num_barriers() const;
  [[nodiscard]] const IITree<bp_t, ExtrusionBarrier>& barriers() const;
  [[nodiscard]] IITree<bp_t, ExtrusionBarrier>& barriers();
  [[nodiscard]] absl::Span<const bed_tree_value_t> get_features() const;
  bool allocate_contacts(bp_t bin_size, bp_t diagonal_width);
  bool deallocate_contacts();
  [[nodiscard]] const contact_matrix_t& contacts() const;
  [[nodiscard]] contact_matrix_t& contacts();
  [[nodiscard]] std::shared_ptr<const contact_matrix_t> contacts_ptr() const;
  [[nodiscard]] std::shared_ptr<contact_matrix_t> contacts_ptr();
  [[nodiscard]] u64 hash(XXH3_state_t* xxh_state, u64 seed, usize cell_id) const;
  [[nodiscard]] u64 hash(u64 seed, usize cell_id) const;

  template <typename H>
  inline friend H AbslHashValue(H h, const Chromosome& c);

 private:
  std::string _name{};
  bp_t _start{(std::numeric_limits<bp_t>::max)()};
  bp_t _end{(std::numeric_limits<bp_t>::max)()};
  bp_t _size{(std::numeric_limits<bp_t>::max)()};
  usize _id{(std::numeric_limits<usize>::max)()};
  IITree<bp_t, ExtrusionBarrier> _barriers{};
  // Protect _contacts from concurrent writes and allocations/deallocations
  std::shared_mutex _contacts_mtx{};
  std::shared_ptr<contact_matrix_t> _contacts{};

  std::vector<bed_tree_value_t> _features{};
  bool _ok{true};
};

class Genome {
 public:
  Genome() = default;
  Genome(const boost::filesystem::path& path_to_chrom_sizes,
         const boost::filesystem::path& path_to_extr_barriers,
         const boost::filesystem::path& path_to_chrom_subranges,
         absl::Span<const boost::filesystem::path> paths_to_extra_features,
         double ctcf_prob_occ_to_occ, double ctcf_prob_nocc_to_nocc, bool keep_all_chroms);

  using iterator = absl::btree_set<Chromosome>::iterator;
  using const_iterator = absl::btree_set<Chromosome>::const_iterator;

  [[nodiscard]] iterator begin();
  [[nodiscard]] iterator end();

  [[nodiscard]] const_iterator begin() const;
  [[nodiscard]] const_iterator end() const;

  [[nodiscard]] const_iterator cbegin() const;
  [[nodiscard]] const_iterator cend() const;

  // NOTE: The find and contains method taking a string_view as input parameter are performing a
  // linear search through the genome
  [[nodiscard]] iterator find(const Chromosome& other_chrom);
  [[nodiscard]] const_iterator find(const Chromosome& other_chrom) const;
  [[nodiscard]] iterator find(std::string_view other_chrom_name);
  [[nodiscard]] const_iterator find(std::string_view other_chrom_name) const;

  [[nodiscard]] bool contains(const Chromosome& other_chromosome) const;
  [[nodiscard]] bool contains(std::string_view other_chrom_name) const;

  [[nodiscard]] usize size() const;
  [[nodiscard]] usize number_of_chromosomes() const;
  [[nodiscard]] usize simulated_size() const;

  [[nodiscard]] const Chromosome& chromosome_with_longest_name() const;
  [[nodiscard]] const Chromosome& longest_chromosome() const;
  [[nodiscard]] const Chromosome& chromosome_with_max_nbarriers() const;
  [[nodiscard]] usize max_target_contacts(usize bin_size, usize diagonal_width,
                                          double target_contact_density,
                                          usize simulation_iterations,
                                          double lef_fraction_contact_sampling,
                                          double nlefs_per_mbp, usize ncells) const;

  /// A simple wrapper function that imports chromosomes and extrusion barriers that comprise the
  /// genome that is being simulated.
  [[nodiscard]] static absl::btree_set<Chromosome> instantiate_genome(
      const boost::filesystem::path& path_to_chrom_sizes,
      const boost::filesystem::path& path_to_extr_barriers,
      const boost::filesystem::path& path_to_chrom_subranges,
      absl::Span<const boost::filesystem::path> paths_to_extra_features,
      double ctcf_prob_occ_to_occ, double ctcf_prob_nocc_to_nocc, bool keep_all_chroms);

 private:
  absl::btree_set<Chromosome> _chromosomes{};

  /// Import chromosomes from a chrom.sizes file

  //! When \p path_to_extr_barriers is non-empty, import the intersection of the chromosomes present
  //! in the chrom.sizes and BED files. The optional BED file can be used to instruct MoDLE to
  //! simulate loop extrusion on a sub-region of a chromosome from the chrom.sizes file.
  [[nodiscard]] static absl::btree_set<Chromosome> import_chromosomes(
      const boost::filesystem::path& path_to_chrom_sizes,
      const boost::filesystem::path& path_to_chrom_subranges, bool keep_all_chroms);

  /// Parse a BED file containing the genomic coordinates of extrusion barriers and add them to the
  /// Genome
  static usize import_barriers(absl::btree_set<Chromosome>& chromosomes,
                               const boost::filesystem::path& path_to_extr_barriers,
                               double ctcf_prob_occ_to_occ, double ctcf_prob_nocc_to_nocc);

  /// Parse a BED file containing the genomic coordinates of extra features (e.g. promoters,
  /// enhancer) them to the Genome
  static usize import_extra_features(absl::btree_set<Chromosome>& chromosomes,
                                     const boost::filesystem::path& path_to_extra_features);
};

}  // namespace modle

#include "../../genome_impl.hpp"  // IWYU pragma: export
