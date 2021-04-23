#pragma once

#include <absl/container/btree_set.h>  // for btree_set

#include <cstddef>      // IWYU pragma: keep for size_t
#include <cstdint>      // for uint32_t
#include <memory>       // for shared_ptr
#include <string>       // for string
#include <string_view>  // for string_view
#include <vector>       // for vector

#include "modle/bed.hpp"     // for BED
#include "modle/common.hpp"  // for bp_t, contacts_t
#include "modle/utils.hpp"   // for ndebug_defined

namespace modle {
template <typename I>
class ContactMatrix;
namespace chrom_sizes {
struct ChromSize;
}

class Chromosome {
 public:
  inline explicit Chromosome(const chrom_sizes::ChromSize& chrom,
                             const std::vector<bed::BED>& barriers = {});
  inline explicit Chromosome(chrom_sizes::ChromSize&& chrom, std::vector<bed::BED>&& barriers = {});

  [[nodiscard]] inline bool operator==(const Chromosome& other) const
      noexcept(utils::ndebug_defined());
  [[nodiscard]] inline bool operator==(std::string_view other) const
      noexcept(utils::ndebug_defined());
  [[nodiscard]] inline bool operator<(const Chromosome& other) const
      noexcept(utils::ndebug_defined());

  inline void instantiate_contact_matrix(size_t bin_size, size_t diagonal_width);
  inline void clear_contacts();
  inline void add_extrusion_barrier(bed::BED&& barrier);
  inline void add_extrusion_barrier(bed::BED barrier);

  [[nodiscard]] inline std::string_view name() const;
  [[nodiscard]] inline constexpr bp_t start_pos() const;
  [[nodiscard]] inline constexpr bp_t end_pos() const;
  [[nodiscard]] inline constexpr bp_t size() const;
  [[nodiscard]] inline constexpr bp_t simulated_size() const;
  [[nodiscard]] inline bool ok() const;
  [[nodiscard]] inline size_t nbarriers() const;
  [[nodiscard]] inline const absl::btree_set<bed::BED>& get_barriers() const;
  template <typename I>
  inline void increment_contacts(bp_t pos1, bp_t pos2, bp_t bin_size, I n = 1);
  inline void increment_contacts(bp_t pos1, bp_t pos2, bp_t bin_size);
  inline void allocate_contacts(bp_t bin_size, bp_t diagonal_width);
  inline void deallocate_contacts();
  [[nodiscard]] inline const ContactMatrix<contacts_t>& contacts() const;

  template <typename H>
  inline friend H AbslHashValue(H h, const Chromosome& c);

  struct Comparator {
    using is_transparent = void;
    [[nodiscard]] inline bool operator()(const Chromosome& c1, const Chromosome& c2) const
        noexcept(utils::ndebug_defined());
    [[nodiscard]] inline bool operator()(const Chromosome& c1, std::string_view c2) const
        noexcept(utils::ndebug_defined());
    [[nodiscard]] inline bool operator()(std::string_view c1, const Chromosome& c2) const
        noexcept(utils::ndebug_defined());
    [[nodiscard]] inline bool operator()(std::string_view c1, std::string_view c2) const
        noexcept(utils::ndebug_defined());
  };

 private:
  std::string _name;
  size_t _size;
  size_t _start;
  size_t _end;
  absl::btree_set<bed::BED> _barriers{};
  std::shared_ptr<ContactMatrix<contacts_t>> _contacts{nullptr};
};

}  // namespace modle

#include "../../dna_impl.hpp"  // IWYU pragma: keep
