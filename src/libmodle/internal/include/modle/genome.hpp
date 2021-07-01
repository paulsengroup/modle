#pragma once

#include <absl/container/btree_set.h>  // for btree_set
#include <absl/types/span.h>           // for Span
#include <xxh3.h>                      // for XXH3_state_t

#include <boost/type_traits.hpp>
#include <cstddef>  // IWYU pragma: keep for size_t
#include <cstdint>  // for uint32_t
#include <iterator>
#include <limits>       // for numeric_limits
#include <memory>       // for shared_ptr
#include <string>       // for string
#include <string_view>  // for string_view
#include <vector>       // for vector

#include "modle/bed.hpp"            // for BED
#include "modle/common/common.hpp"  // for bp_t, contacts_t
#include "modle/common/utils.hpp"   // for ndebug_defined, XXH3_Deleter
#include "modle/interval_tree.hpp"  // for IITree

namespace modle {

template <typename I>
class ContactMatrix;

class ExtrusionBarrier;

class Chromosome {
  using contact_matrix_t = ContactMatrix<contacts_t>;
  using bed_tree_value_t = bed::BED_tree<>::value_type;
  friend class Genome;

 public:
  Chromosome() = default;
  template <typename I, typename = std::enable_if_t<std::is_integral_v<I>>>
  inline Chromosome(size_t id, std::string_view chrom_name, I chrom_start, I chrom_end,
                    I chrom_size);

  Chromosome(size_t id, const bed::BED& chrom);
  Chromosome(size_t id, const bed::BED& chrom, const IITree<bp_t, ExtrusionBarrier>& barriers);
  Chromosome(size_t id, const bed::BED& chrom, IITree<bp_t, ExtrusionBarrier>&& barriers);

  template <typename I = uint64_t, typename = std::enable_if_t<std::is_integral_v<I>>>
  inline Chromosome(size_t id, std::string_view chrom_name, I chrom_start, I chrom_end,
                    I chrom_size, const IITree<bp_t, ExtrusionBarrier>& barriers);
  template <typename I = uint64_t, typename = std::enable_if_t<std::is_integral_v<I>>>
  inline Chromosome(size_t id, std::string_view chrom_name, I chrom_start, I chrom_end,
                    I chrom_size, IITree<bp_t, ExtrusionBarrier>&& barriers);

  ~Chromosome() = default;

  Chromosome(const Chromosome& other);
  Chromosome(Chromosome&& other) = default;

  Chromosome& operator=(const Chromosome& other);
  Chromosome& operator=(Chromosome&& other) = default;

  [[nodiscard]] bool operator==(const Chromosome& other) const noexcept(utils::ndebug_defined());
  [[nodiscard]] bool operator==(std::string_view other_name) const noexcept;
  [[nodiscard]] bool operator<(const Chromosome& other) const noexcept(utils::ndebug_defined());

  void add_extrusion_barrier(const bed::BED& record, double ctcf_prob_occ_to_occ,
                             double ctcf_prob_nocc_to_nocc);

  template <typename Iter,
            typename = std::enable_if_t<std::is_same_v<
                std::remove_cv_t<typename std::iterator_traits<Iter>::value_type>, bed::BED>>>
  inline void add_extrusion_barrier(Iter barriers_begin, Iter barriers_end);

  [[nodiscard]] size_t id() const;
  [[nodiscard]] std::string_view name() const;
  [[nodiscard]] const char* name_cstr() const;
  [[nodiscard]] constexpr bp_t start_pos() const;
  [[nodiscard]] constexpr bp_t end_pos() const;
  [[nodiscard]] constexpr bp_t size() const;
  [[nodiscard]] constexpr bp_t simulated_size() const;
  [[nodiscard]] bool ok() const;
  [[nodiscard]] size_t num_lefs(double nlefs_per_mbp) const;
  [[nodiscard]] size_t num_barriers() const;
  [[nodiscard]] const IITree<bp_t, ExtrusionBarrier>& barriers() const;
  [[nodiscard]] IITree<bp_t, ExtrusionBarrier>& barriers();
  [[nodiscard]] absl::Span<const bed_tree_value_t> get_features() const;
  void increment_contacts(bp_t pos1, bp_t pos2, bp_t bin_size);
  void increment_contacts(bp_t bin1, bp_t bin2);
  void allocate_contacts(bp_t bin_size, bp_t diagonal_width);
  void deallocate_contacts();
  [[nodiscard]] const contact_matrix_t& contacts() const;
  [[nodiscard]] contact_matrix_t& contacts();
  [[nodiscard]] uint64_t hash(uint64_t seed, size_t cell_id = 0);

  template <typename H>
  inline friend H AbslHashValue(H h, const Chromosome& c);

 private:
  std::string _name{};
  bp_t _start{std::numeric_limits<bp_t>::max()};
  bp_t _end{std::numeric_limits<bp_t>::max()};
  bp_t _size{std::numeric_limits<bp_t>::max()};
  size_t _id{std::numeric_limits<size_t>::max()};
  IITree<bp_t, ExtrusionBarrier> _barriers{};
  std::shared_ptr<contact_matrix_t> _contacts{nullptr};
  std::vector<bed_tree_value_t> _features{};

  std::unique_ptr<XXH3_state_t, utils::XXH3_Deleter> _xxh_state{XXH3_createState()};
};

class Genome {
  using Genome_t = bed::BED_tree<bed::BED, uint32_t>;

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

  [[nodiscard]] size_t size() const;
  [[nodiscard]] size_t simulated_size() const;

  [[nodiscard]] const Chromosome& chromosome_with_longest_name() const;
  [[nodiscard]] const Chromosome& longest_chromosome() const;
  [[nodiscard]] const Chromosome& chromosome_with_max_nbarriers() const;
  [[nodiscard]] size_t max_target_contacts(size_t bin_size, size_t diagonal_width,
                                           double target_contact_density,
                                           size_t simulation_iterations,
                                           double lef_fraction_contact_sampling,
                                           double nlefs_per_mbp, size_t ncells) const;

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
  //! in the chrom.sizes and BED files. The optional BED file can be used to instruct ModLE to
  //! simulate loop extrusion on a sub-region of a chromosome from the chrom.sizes file.
  [[nodiscard]] static absl::btree_set<Chromosome> import_chromosomes(
      const boost::filesystem::path& path_to_chrom_sizes,
      const boost::filesystem::path& path_to_chrom_subranges, bool keep_all_chroms);

  /// Parse a BED file containing the genomic coordinates of extrusion barriers and add them to the
  /// Genome
  static size_t import_barriers(absl::btree_set<Chromosome>& chromosomes,
                                const boost::filesystem::path& path_to_extr_barriers,
                                double ctcf_prob_occ_to_occ, double ctcf_prob_nocc_to_nocc);

  /// Parse a BED file containing the genomic coordinates of extra features (e.g. promoters,
  /// enhancer) them to the Genome
  static size_t import_extra_features(absl::btree_set<Chromosome>& chromosomes,
                                      const boost::filesystem::path& path_to_extra_features);
};

}  // namespace modle

#include "../../genome_impl.hpp"  // IWYU pragma: keep
