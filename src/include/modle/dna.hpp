#ifndef MODLE_DNA_HPP
#define MODLE_DNA_HPP

#include <random>
#include <string_view>
#include <vector>

#include "absl/container/inlined_vector.h"
#include "modle/contacts.hpp"

namespace modle {

class ExtrusionBarrier;
class ExtrusionUnit;

class DNA {
 public:
  enum Direction { none = 0, fwd = 1, rev = 2, both = 3 };

  class Bin {
    friend class DNA;

   public:
    Bin(uint32_t idx, uint64_t start, uint64_t end, std::vector<ExtrusionBarrier>& barriers);
    Bin(uint32_t idx, uint64_t start, uint64_t end);

    [[nodiscard]] uint32_t get_start() const;
    [[nodiscard]] uint32_t get_end() const;
    [[nodiscard]] uint32_t size() const;
    [[nodiscard]] uint32_t get_n_extr_units() const;
    [[nodiscard]] uint32_t get_index() const;
    [[nodiscard]] ExtrusionBarrier* get_next_extr_barrier(ExtrusionBarrier* b,
                                                          Direction d = DNA::Direction::both) const;
    [[nodiscard]] std::vector<ExtrusionBarrier>* get_all_extr_barriers() const;
    void add_extr_barrier(ExtrusionBarrier b);
    void add_extr_barrier(double prob_of_barrier_block, DNA::Direction direction);
    uint32_t add_extr_unit_binding(ExtrusionUnit* unit);
    uint32_t remove_extr_unit_binding(ExtrusionUnit* unit);
    [[nodiscard]] absl::InlinedVector<ExtrusionUnit*, 10>& get_extr_units();

   private:
    uint32_t _idx;
    uint32_t _start;
    uint32_t _end;
    std::unique_ptr<std::vector<ExtrusionBarrier>> _extr_barriers{nullptr};
    // TODO: Consider making this a vector
    std::unique_ptr<absl::InlinedVector<ExtrusionUnit*, 10>> _extr_units{nullptr};
    void add_barrier(ExtrusionBarrier& b);
    void remove_barrier(Direction d);
    void remove_barriers();
  };

  DNA(uint64_t length, uint32_t bin_size);

  // Getters
  [[nodiscard]] uint32_t length() const;
  [[nodiscard]] uint32_t get_n_bins() const;
  [[nodiscard]] uint32_t get_n_barriers() const;
  [[nodiscard]] uint32_t get_bin_size() const;
  [[nodiscard]] DNA::Bin& get_bin_from_pos(uint32_t pos);
  [[nodiscard]] DNA::Bin* get_ptr_to_bin_from_pos(uint32_t pos);
  [[nodiscard]] DNA::Bin& get_prev_bin(const Bin& current_bin);
  [[nodiscard]] DNA::Bin* get_ptr_to_prev_bin(const Bin& current_bin);
  [[nodiscard]] DNA::Bin& get_next_bin(const Bin& current_bin);
  [[nodiscard]] DNA::Bin* get_ptr_to_next_bin(const Bin& current_bin);
  [[nodiscard]] DNA::Bin& get_first_bin();
  [[nodiscard]] DNA::Bin* get_ptr_to_first_bin();
  [[nodiscard]] DNA::Bin& get_last_bin();
  [[nodiscard]] DNA::Bin* get_ptr_to_last_bin();

  // Iterator stuff
  [[nodiscard]] std::vector<Bin>::iterator begin();
  [[nodiscard]] std::vector<Bin>::iterator end();
  [[nodiscard]] std::vector<Bin>::const_iterator cbegin() const;
  [[nodiscard]] std::vector<Bin>::const_iterator cend() const;

  // Modifiers
  void add_barrier(ExtrusionBarrier& b, uint32_t pos);
  void remove_barrier(uint32_t pos, Direction direction);

 private:
  std::vector<Bin> _bins;
  uint64_t _length;
  uint32_t _bin_size;

  // Initializer
  static std::vector<Bin> make_bins(uint64_t length, uint32_t bin_size);
};

struct Chromosome {
  Chromosome(std::string name, uint64_t length, uint32_t bin_size, uint32_t avg_lef_processivity);
  [[nodiscard]] uint32_t length() const;
  [[nodiscard]] uint32_t n_bins() const;
  [[nodiscard]] uint32_t n_barriers() const;
  void write_contacts_to_tsv(const std::string& path_to_file, bool complete = false) const;
  std::string name;
  DNA dna;
  std::vector<ExtrusionBarrier*> barriers;
  ContactMatrix<uint32_t> contacts;
};

};  // namespace modle

#endif  // MODLE_DNA_HPP
