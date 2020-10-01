#ifndef MODLE_DNA_HPP
#define MODLE_DNA_HPP

#include <string_view>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
//#include "modle/lefs.hpp"

namespace modle {

class ExtrusionBarrier;
class ExtrusionUnit;

class DNA {
 public:
  class Bin {
    friend class DNA;

   public:
    Bin(uint32_t idx, uint64_t start, uint64_t end, ExtrusionBarrier const* fwd_barrier,
        ExtrusionBarrier const* rev_barrier);
    Bin(uint32_t idx, uint64_t start, uint64_t end);

    [[nodiscard]] bool has_fwd_barrier() const;
    [[nodiscard]] bool has_rev_barrier() const;
    void add_fwd_barrier(const ExtrusionBarrier* barrier);
    void add_rev_barrier(const ExtrusionBarrier* barrier);
    void remove_fwd_barrier();
    void remove_rev_barrier();
    void remove_barriers();
    [[nodiscard]] uint32_t get_start() const;
    [[nodiscard]] uint32_t get_end() const;
    [[nodiscard]] uint32_t get_center() const;
    [[nodiscard]] uint32_t size() const;
    // The following two functions return the number of extr. units that are bound to a given bin
    [[nodiscard]] uint32_t add_extr_unit_binding(ExtrusionUnit* unit);
    [[nodiscard]] uint32_t remove_extr_unit_binding(ExtrusionUnit* unit);
    [[nodiscard]] bool pos_is_occupied(uint32_t pos) const;
    [[nodiscard]] uint32_t n_extruders() const;
    [[nodiscard]] absl::flat_hash_set<ExtrusionUnit*>& get_extruders();

   private:
    uint32_t _idx;
    uint32_t _start;
    uint32_t _end;
    ExtrusionBarrier const* _fwd_barrier;
    ExtrusionBarrier const* _rev_barrier;
    std::unique_ptr<absl::flat_hash_set<ExtrusionUnit*>> _extr_units{nullptr};
  };

  DNA(uint64_t length, uint32_t bin_size);
  void add_fwd_barrier(const ExtrusionBarrier& barrier, uint32_t pos);
  void add_rev_barrier(const ExtrusionBarrier& barrier, uint32_t pos);
  void remove_fwd_barrier(uint32_t pos);
  void remove_rev_barrier(uint32_t pos);
  [[nodiscard]] uint32_t length() const;
  [[nodiscard]] uint32_t n_bins() const;
  [[nodiscard]] uint32_t n_barriers() const;
  [[nodiscard]] uint32_t bin_size() const;
  [[nodiscard]] std::shared_ptr<DNA::Bin> get_ptr_to_bin_from_pos(uint32_t pos);
  [[nodiscard]] std::vector<std::shared_ptr<Bin>>::iterator begin();
  [[nodiscard]] std::vector<std::shared_ptr<Bin>>::iterator end();
  [[nodiscard]] std::vector<std::shared_ptr<Bin>>::const_iterator cbegin() const;
  [[nodiscard]] std::vector<std::shared_ptr<Bin>>::const_iterator cend() const;
  [[nodiscard]] std::shared_ptr<DNA::Bin> get_ptr_to_previous_bin(const DNA::Bin* current_bin);
  [[nodiscard]] std::shared_ptr<DNA::Bin> get_ptr_to_next_bin(const DNA::Bin* current_bin);

 private:
  std::vector<std::shared_ptr<Bin>> _bins;
  uint64_t _length;
  static std::vector<std::shared_ptr<Bin>> make_bins(uint64_t length, uint32_t bin_size);
};
};  // namespace modle

#endif  // MODLE_DNA_HPP
