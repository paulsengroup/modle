#ifndef MODLE_DNA_HPP
#define MODLE_DNA_HPP

#include <string_view>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "modle/extr_barrier.hpp"
#include "modle/parsers.hpp"
//#include "modle/lefs.hpp"

namespace modle {

class DNA {
 public:
  class Bin {
    friend class DNA;

   public:
    Bin(uint32_t idx, uint32_t start, uint32_t end, std::shared_ptr<ExtrusionBarrier> fwd_barrier,
        std::shared_ptr<ExtrusionBarrier> rev_barrier);
    Bin(uint32_t idx, uint32_t start, uint32_t end);

    [[nodiscard]] bool has_fwd_barrier() const;
    [[nodiscard]] bool has_rev_barrier() const;
    void add_fwd_barrier(std::shared_ptr<ExtrusionBarrier> barrier);
    void add_rev_barrier(std::shared_ptr<ExtrusionBarrier> barrier);
    void remove_fwd_barrier();
    void remove_rev_barrier();
    void remove_barriers();
    [[nodiscard]] uint32_t get_start() const;
    [[nodiscard]] uint32_t get_end() const;
    [[nodiscard]] uint32_t get_center() const;
    [[nodiscard]] uint32_t size() const;

   private:
    uint32_t _idx;
    uint32_t _start;
    uint32_t _end;
    std::shared_ptr<ExtrusionBarrier> _fwd_barrier;
    std::shared_ptr<ExtrusionBarrier> _rev_barrier;
    //      absl::flat_hash_set<std::shared_ptr<Lef>> _lefs{};
  };

  DNA(uint32_t length, uint32_t bin_size);
  void add_fwd_barrier(std::shared_ptr<ExtrusionBarrier> barrier, uint32_t pos);
  void add_rev_barrier(std::shared_ptr<ExtrusionBarrier> barrier, uint32_t pos);
  void remove_fwd_barrier(uint32_t pos);
  void remove_rev_barrier(uint32_t pos);
  [[nodiscard]] uint32_t length() const;
  [[nodiscard]] uint32_t nbins() const;
  [[nodiscard]] uint32_t bin_size() const;
  [[nodiscard]] std::shared_ptr<DNA::Bin> get_ptr_to_bin_from_pos(uint32_t pos);
  [[nodiscard]] std::vector<std::shared_ptr<Bin>>::iterator begin();
  [[nodiscard]] std::vector<std::shared_ptr<Bin>>::iterator end();
  [[nodiscard]] std::vector<std::shared_ptr<Bin>>::const_iterator cbegin() const;
  [[nodiscard]] std::vector<std::shared_ptr<Bin>>::const_iterator cend() const;
  [[nodiscard]] std::shared_ptr<DNA::Bin> get_ptr_to_previous_bin(const std::shared_ptr<DNA::Bin>& current_bin);
  [[nodiscard]] std::shared_ptr<DNA::Bin> get_ptr_to_next_bin(const std::shared_ptr<DNA::Bin>& current_bin);

 private:
  std::vector<std::shared_ptr<Bin>> _bins;
  uint32_t _length;
  static std::vector<std::shared_ptr<Bin>> make_bins(uint32_t length, uint32_t bin_size);
};
};  // namespace modle

#endif  // MODLE_DNA_HPP
