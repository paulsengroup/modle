#pragma once

#include <absl/container/btree_set.h>  // for btree_set
#include <xxh3.h>

#include <cstddef>      // IWYU pragma: keep for size_t
#include <cstdint>      // for uint32_t
#include <limits>       // for numeric_limits
#include <memory>       // for shared_ptr
#include <string>       // for string
#include <string_view>  // for string_view
#include <vector>       // for vector

#include "modle/bed.hpp"     // for BED
#include "modle/common.hpp"  // for bp_t, contacts_t
#include "modle/contacts.hpp"
#include "modle/utils.hpp"  // for ndebug_defined

namespace modle {
namespace chrom_sizes {
struct ChromSize;
}

class Chromosome {
 public:
  explicit Chromosome(const chrom_sizes::ChromSize& chrom,
                      size_t id = std::numeric_limits<size_t>::max(),
                      const std::vector<bed::BED>& barriers = {});
  explicit Chromosome(chrom_sizes::ChromSize&& chrom,
                      size_t id = std::numeric_limits<size_t>::max(),
                      std::vector<bed::BED>&& barriers = {});

  [[nodiscard]] bool operator==(const Chromosome& other) const noexcept(utils::ndebug_defined());
  [[nodiscard]] bool operator==(std::string_view other) const noexcept(utils::ndebug_defined());
  [[nodiscard]] bool operator<(const Chromosome& other) const noexcept(utils::ndebug_defined());

  void add_extrusion_barrier(bed::BED&& barrier);
  void add_extrusion_barrier(bed::BED barrier);

  [[nodiscard]] size_t id() const;
  [[nodiscard]] std::string_view name() const;
  [[nodiscard]] const char* name_cstr() const;
  [[nodiscard]] constexpr bp_t start_pos() const;
  [[nodiscard]] constexpr bp_t end_pos() const;
  [[nodiscard]] constexpr bp_t size() const;
  [[nodiscard]] constexpr bp_t simulated_size() const;
  [[nodiscard]] bool ok() const;
  [[nodiscard]] size_t nlefs(double nlefs_per_mbp) const;
  [[nodiscard]] size_t nbarriers() const;
  [[nodiscard]] size_t num_valid_barriers() const;
  [[nodiscard]] const absl::btree_set<bed::BED>& get_barriers() const;
  void increment_contacts(bp_t pos1, bp_t pos2, bp_t bin_size);
  void increment_contacts(bp_t bin1, bp_t bin2);
  void allocate_contacts(bp_t bin_size, bp_t diagonal_width);
  void deallocate_contacts();
  [[nodiscard]] const ContactMatrix<contacts_t>& contacts() const;
  [[nodiscard]] uint64_t hash(uint64_t seed, size_t cell_id = 0);

  template <typename H>
  inline friend H AbslHashValue(H h, const Chromosome& c);

  struct Comparator {
    using is_transparent = void;
    [[nodiscard]] bool operator()(const Chromosome& c1, const Chromosome& c2) const
        noexcept(utils::ndebug_defined());
    [[nodiscard]] bool operator()(const Chromosome& c1, std::string_view c2) const
        noexcept(utils::ndebug_defined());
    [[nodiscard]] bool operator()(std::string_view c1, const Chromosome& c2) const
        noexcept(utils::ndebug_defined());
    [[nodiscard]] bool operator()(std::string_view c1, std::string_view c2) const
        noexcept(utils::ndebug_defined());
  };

 private:
  size_t _id{std::numeric_limits<size_t>::max()};
  std::string _name;
  size_t _size;
  size_t _start;
  size_t _end;
  absl::btree_set<bed::BED> _barriers{};
  std::shared_ptr<ContactMatrix<contacts_t>> _contacts{nullptr};

  struct XXH3_Deleter {
    void operator()(XXH3_state_t* state) noexcept { XXH3_freeState(state); }
  };

  std::unique_ptr<XXH3_state_t, XXH3_Deleter> _xxh_state{XXH3_createState()};
};

}  // namespace modle

#include "../../dna_impl.hpp"
