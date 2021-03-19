#pragma once

#include <absl/container/btree_set.h>  // for btree_set

#include <cstddef>      // IWYU pragma: keep for size_t
#include <cstdint>      // for uint32_t
#include <memory>       // for shared_ptr
#include <string>       // for string
#include <string_view>  // for string_view
#include <vector>       // for vector

#include "modle/bed.hpp"     // for BED
#include "modle/common.hpp"  // for Bp, Contacts

namespace modle {
template <typename I>
class ContactMatrix;
namespace chr_sizes {
struct ChrSize;
}

class Chromosome {
 public:
  inline explicit Chromosome(const chr_sizes::ChrSize& chrom,
                             const std::vector<bed::BED>& barriers = {});
  inline explicit Chromosome(chr_sizes::ChrSize&& chrom, std::vector<bed::BED>&& barriers = {});

  [[nodiscard]] inline bool operator==(const Chromosome& other) const;
  [[nodiscard]] inline bool operator==(std::string_view other) const;
  [[nodiscard]] inline bool operator<(const Chromosome& other) const;

  inline void instantiate_contact_matrix(std::size_t bin_size, std::size_t diagonal_width);
  inline void clear_contacts();
  inline void add_extrusion_barrier(bed::BED&& barrier);
  inline void add_extrusion_barrier(bed::BED barrier);

  [[nodiscard]] inline std::string_view name() const;
  [[nodiscard]] inline Bp start_pos() const;
  [[nodiscard]] inline Bp end_pos() const;
  [[nodiscard]] inline Bp size() const;
  [[nodiscard]] inline Bp simulated_size() const;
  [[nodiscard]] inline bool ok() const;
  [[nodiscard]] inline std::size_t nbarriers() const;
  [[nodiscard]] inline const absl::btree_set<bed::BED>& get_barriers() const;
  template <typename I>
  inline void increment_contacts(Bp pos1, Bp pos2, Bp bin_size, I n = 1);
  inline void increment_contacts(Bp pos1, Bp pos2, Bp bin_size);
  inline void allocate_contacts(Bp bin_size, Bp diagonal_width);
  inline void deallocate_contacts();
  [[nodiscard]] inline const ContactMatrix<Contacts>& contacts() const;

  template <typename H>
  inline friend H AbslHashValue(H h, const Chromosome& c);

  struct Comparator {
    using is_transparent = void;
    [[nodiscard]] inline bool operator()(const Chromosome& c1, const Chromosome& c2) const;
    [[nodiscard]] inline bool operator()(const Chromosome& c1, std::string_view c2) const;
    [[nodiscard]] inline bool operator()(std::string_view c1, const Chromosome& c2) const;
    [[nodiscard]] inline bool operator()(std::string_view c1, std::string_view c2) const;
  };

 private:
  std::string _name;
  std::size_t _size;
  std::size_t _start;
  std::size_t _end;
  absl::btree_set<bed::BED> _barriers{};
  std::shared_ptr<ContactMatrix<Contacts>> _contacts{nullptr};
};

}  // namespace modle

#include "../../dna_impl.hpp"  // IWYU pragma: keep
