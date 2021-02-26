#pragma once

#include <absl/container/flat_hash_set.h>
#include <absl/container/inlined_vector.h>
#include <absl/types/span.h>

#include <cstdint>  // for uint*_t
#include <iosfwd>   // for size_t
#include <memory>   // unique_ptr
#include <string>
#include <string_view>
#include <vector>
#ifdef USE_XOSHIRO
#include <XoshiroCpp.hpp>
#endif

#include "modle/common.hpp"
#include "modle/contacts.hpp"  // for ContactMatrix
#include "modle/extr_barrier.hpp"

namespace modle {

#ifdef USE_XOSHIRO
using PRNG = XoshiroCpp::Xoshiro256PlusPlus;
using seeder = XoshiroCpp::SplitMix64;
#else
using PRNG = std::mt19937_64;
using seeder = std::seed_seq;
#endif

// Pre-declarations
class ExtrusionUnit;
class Lef;
namespace bed {
struct BED;
}

/// Class used to represent a DNA molecule as a vector Bin%s.
/** Aside from owning the \p std::vector<Bin>, DNA is also used to keep track of the
 * simulated_length and bin size of the DNA molecule.
 */
class DNA {
 public:
  class Bin;
  /** Enumerator used across ModLE's code-base to describe directionality in terms of DNA strand
   * (+/- or forward/reverse).
   */
  template <typename I1, typename I2>
  inline DNA(I1 length, I2 bin_size, bool allocate_bins = false);

  // Getters
  [[nodiscard]] inline std::size_t size() const;
  [[nodiscard]] inline uint32_t get_bin_size() const;
  [[nodiscard]] inline double get_total_lef_affinity() const;

  [[nodiscard]] inline std::size_t get_n_bins() const;
  [[nodiscard]] inline std::size_t get_n_barriers() const;

  [[nodiscard]] inline DNA::Bin* get_ptr_to_bin(std::size_t pos);
  [[nodiscard]] inline DNA::Bin& get_bin(std::size_t pos);

  [[nodiscard]] inline DNA::Bin* get_ptr_to_prev_bin(const Bin& current_bin);
  [[nodiscard]] inline DNA::Bin& get_prev_bin(const Bin& current_bin);

  [[nodiscard]] inline DNA::Bin* get_ptr_to_next_bin(const Bin& current_bin);
  [[nodiscard]] inline DNA::Bin& get_next_bin(const Bin& current_bin);

  [[nodiscard]] inline DNA::Bin& get_first_bin();
  [[nodiscard]] inline DNA::Bin* get_ptr_to_first_bin();

  [[nodiscard]] inline DNA::Bin& get_last_bin();
  [[nodiscard]] inline DNA::Bin* get_ptr_to_last_bin();

  [[nodiscard]] inline std::size_t get_bin_idx(std::size_t pos) const;
  [[nodiscard]] inline std::size_t get_bin_idx(const DNA::Bin& bin) const;

  [[nodiscard]] inline DNA::Bin& operator[](std::size_t idx);
  [[nodiscard]] inline const DNA::Bin& operator[](std::size_t idx) const;

  // Iterator stuff
  [[nodiscard]] inline std::vector<Bin>::iterator begin();
  [[nodiscard]] inline std::vector<Bin>::iterator end();
  [[nodiscard]] inline std::vector<Bin>::const_iterator cbegin() const;
  [[nodiscard]] inline std::vector<Bin>::const_iterator cend() const;

  // Modifiers
  inline ExtrusionBarrier* add_extr_barrier(ExtrusionBarrier& b, std::size_t pos);
  inline ExtrusionBarrier* add_extr_barrier(const modle::bed::BED& record);
  inline void remove_extr_barrier(std::size_t pos, dna::Direction direction);
  inline void allocate_bins();

 private:
  std::vector<Bin> _bins{};  ///< \p Bin%s belonging to this instance of DNA.
  std::size_t _length;       ///< DNA molecule simulated_length in bp.
  uint32_t _bin_size;        ///< Target bin size in bp.

  /** Initialize the <tt> std::vector<Bin> </tt> given a simulated_length and bin size.
   *
   * When <tt> simulated_length > bin_size </tt> the <tt> std::vector<Bin> </tt> will contain a
   * single DNA::Bin starting at <tt> pos == 0 </tt> and ending at <tt> pos == simulated_length
   * </tt>.
   *
   * When <tt> simulated_length % bin_size != 0 </tt> then the last bin in <tt>
   * std::vector<DNA::Bin> </tt> will be shorter than  <tt> bin_size </tt>.
   *
   * This will introduce a bias towards the last DNA::Bin when randomly selecting DNA::Bin by index.
   * This bias should be negligible for <tt> simulated_length >> bin_size </tt> (which should be the
   * case under normal circumstances).
   * @param length    Length of DNA molecule (bp).
   * @param bin_size  Bin size (bp).
   * @return          A vector of DNA::Bin with where consecutive DNA::Bin are non-overlapping (i.e.
   * given a pair of consecutive DNA::Bin%s <tt> b1 </tt> and <tt> b2 </tt>, <tt> b1::_end ==
   * b2::_start + 1 </tt>.
   */
  inline static std::vector<Bin> make_bins(std::size_t length, uint32_t bin_size);

 public:
  /// Class to represent a single Bin from a DNA molecule
  /** A DNA::Bin knows of the genomic coordinates corresponding to its DNA::Bin::_start and
   * DNA::Bin::_end positions, as well as the entities that are currently bound to it (currently
   * only ExtrusionBarrier and ExtrusionUnit).
   */
  class Bin {
    friend class DNA;

   public:
    ///@{
    /** DNA::Bin Constructor
     * @param idx        Index in the \p std::vector that owns this instance of DNA::Bin.
     * @param start      Genomic coordinate of the left-most position represented by this DNA::Bin.
     * @param end        Genomic coordinate of the right-most position represented by this DNA::Bin.
     * @param barriers   List of ExtrusionBarrier%s mapped in the genomic region represented by this
     * DNA::Bin.
     */
    inline Bin() = default;
    template <typename I>
    inline Bin(std::size_t idx, I bin_size, absl::Span<ExtrusionBarrier> barriers);
    template <typename I>
    inline Bin(std::size_t idx, I bin_size);
    ///@}

    [[nodiscard]] inline std::size_t size() const;
    [[nodiscard]] inline std::size_t get_n_extr_units() const;
    [[nodiscard]] inline std::size_t get_n_extr_barriers() const;

    [[nodiscard]] inline std::size_t get_index() const;
    [[nodiscard]] inline std::size_t get_start() const;
    [[nodiscard]] inline std::size_t get_end() const;

    [[nodiscard]] inline bool has_extr_barrier() const;
    [[nodiscard]] inline absl::Span<ExtrusionBarrier> get_all_extr_barriers();
    /** Given a ptr to an instance of ExtrusionBarrier that belongs to this instance of
     * DNA::Bin, return a ptr to the previous/next instance of ExtrusionBarrier that is blocking
     * ExtrusionUnit%s that are moving in the direction specified by \p d.
     * @param old_b Pointer to ExtrusionBarrier to use to look for the next ExtrusionBarrier (i.e.
     * DNA::Bin::get_(prev)|(next)_extr_barrier will only search for ExtrusionBarriers that come
     * before/after \p old_b). When <tt> b == nullptr </tt>, then the search will involve all the
     * ExtrusionBarrier that are bound to the current instance of DNA::Bin.
     * @param d dna::Direction of the ExtrusionBarrier to look for (the caller has to make sure that
     * <tt> d != dna::Direction::none </tt>).
     * @return Pointer to the next ExtrusionBarrier with the orientation specified by \p d (or
     * \p nullptr if there are no ExtrusionBarrier%s matching the search parameters).
     */

    [[nodiscard]] inline ExtrusionBarrier* get_ptr_to_prev_extr_barrier(
        const ExtrusionBarrier* old_b, dna::Direction d = dna::Direction::both);
    [[nodiscard]] inline ExtrusionBarrier* get_ptr_to_prev_extr_barrier(
        uint64_t pos, dna::Direction d = dna::Direction::both);
    [[nodiscard]] inline ExtrusionBarrier* get_ptr_to_next_extr_barrier(
        const ExtrusionBarrier* old_b, dna::Direction d = dna::Direction::both);
    [[nodiscard]] inline ExtrusionBarrier* get_ptr_to_next_extr_barrier(
        std::size_t pos, dna::Direction d = dna::Direction::both);

    [[nodiscard]] inline absl::Span<ExtrusionUnit*> get_extr_units();
    [[nodiscard]] inline float get_lef_affinity() const;

    inline ExtrusionBarrier* add_extr_barrier(ExtrusionBarrier b);
    inline ExtrusionBarrier* add_extr_barrier(uint64_t pos, double prob_of_barrier_block,
                                              dna::Direction direction);

    inline void remove_extr_barrier(dna::Direction d);
    inline uint64_t remove_all_extr_barriers();
    /** Register the ExtrusionUnit unit as being bound to this instance of DNA::Bin.
     *
     * This function will take care of allocating the \p std::vector<ExtrusionUnit> when
     * <tt> DNA::Bin::_extr_barriers = nullptr </tt>.
     * @param unit Pointer to the ExtrusionUnit to be registered as binding to this DNA::Bin.
     */
    inline uint64_t register_extr_unit_binding(ExtrusionUnit* unit);
    /** Remove the ExtrusionUnit \p unit from the list of ExtrusionUnit%s that are bound to this
     * instance of DNA::Bin.
     *
     * This function will take care of deallocating the \p std::vector<ExtrusionUnit> when there are
     * no more ExtrusionUnit%s bound to this DNA::Bin.
     * @param unit Pointer to the unit to be removed from the \p std::vector<ExtrusionUnit> bound to
     * this DNA::Bin.
     */
    inline uint64_t unregister_extr_unit_binding(ExtrusionUnit* unit);
    inline void sort_extr_barriers_by_pos();

   private:
    /// Index corresponding to this DNA::Bin in the owning \p std::vector<DNA::Bin> (DNA::_bins).
    std::size_t _idx{};
    uint32_t _bin_size{};
    float _lef_affinity{1.0};
    /// \brief Pointer to a \p std::vector of ExtrusionBarrier%s (or \p nullptr if there are no
    /// ExtrusionBarrier%s).
    absl::InlinedVector<ExtrusionBarrier, 1> _extr_barriers{};
    /// Ptr to a \p std::vector of ExtrusionUnit%s (or \p nullptr if there are no ExtrusionUnit%s).
    absl::InlinedVector<ExtrusionUnit*, 3> _extr_units{};

    // Make DNA::Bin hashable by absl::hash
    template <typename H>
    inline friend H AbslHashValue(H h, const DNA::Bin& b) {
      return H::combine(std::move(h), &b);
    }
  };
};

/// Struct representing a Chromosome.
/** The purpose of this struct is to keep together data regarding a single DNA molecule
 */
struct Chromosome {
  friend class Genome;
  inline Chromosome(std::string name, std::size_t length, uint32_t bin_size,
                    uint32_t diagonal_width_, uint64_t seed = 0, bool allocate = false);
  inline Chromosome(std::string chr_name, std::size_t chr_start, std::size_t chr_end,
                    std::size_t length, uint32_t bin_size, uint32_t diagonal_width_,
                    uint64_t seed = 0, bool allocate = false);

  // Getters
  [[nodiscard]] inline std::size_t simulated_length() const;
  [[nodiscard]] inline std::size_t real_length() const;
  [[nodiscard]] inline std::size_t get_start_pos() const;
  [[nodiscard]] inline std::size_t get_end_pos() const;
  [[nodiscard]] inline std::size_t get_nbins() const;
  [[nodiscard]] inline uint32_t get_bin_size() const;
  [[nodiscard]] inline std::size_t get_nbarriers() const;
  [[nodiscard]] inline std::size_t get_nlefs() const;
  [[nodiscard]] inline double get_total_lef_affinity() const;

  // inline void write_barriers_to_tsv(std::string_view path_to_outfile, bool force_overwrite)
  // const;
  inline void allocate_contacts();
  inline void allocate();
  [[nodiscard]] inline bool ok() const;

  std::string name;            // NOLINT
  std::size_t start;           // NOLINT
  std::size_t end;             // NOLINT
  std::size_t total_length;    // NOLINT
  std::size_t diagonal_width;  // NOLINT
  DNA dna;                     // NOLINT
  /// Pointers to the ExtrusionBarrier associated to \p Chromosome::dna.
  std::vector<ExtrusionBarrier*> barriers{};  // NOLINT
  std::vector<Lef*> lefs{};                   // NOLINT
  ContactMatrix<uint32_t> contacts{};         // NOLINT

 private:
  modle::PRNG _rand_eng;
  uint64_t _seed;
  bool _ok{true};
};

}  // namespace modle

#include "../../dna_impl.hpp"
