#include "modle/lefs.hpp"

#include "fmt/printf.h"

namespace modle {

ExtrusionUnit::ExtrusionUnit(Lef& lef, double prob_of_extr_unit_bypass)
    : _parent_lef(lef), _n_stall_generator(prob_of_extr_unit_bypass) {}

uint32_t ExtrusionUnit::get_pos() const {
  assert(this->_bin);  // NOLINT
  return (this->_bin->get_end() + this->_bin->get_start()) / 2;
}

DNA::Direction ExtrusionUnit::get_extr_direction() const { return this->_direction; }

bool ExtrusionUnit::is_stalled() const { return this->_stalls_left > 0; }

double ExtrusionUnit::get_prob_of_extr_unit_bypass() const {
  return this->_parent_lef.get_probability_of_extr_unit_bypass();
}

std::size_t ExtrusionUnit::get_bin_index() const { return this->_bin->get_index(); }

bool ExtrusionUnit::try_extrude(std::mt19937& rand_eng) {
  assert(this->_direction == DNA::Direction::fwd ||  // NOLINT
         this->_direction == DNA::Direction::rev);   // NOLINT
  assert(this->_bin != nullptr);                     // NOLINT
  if (this->is_stalled()) {
    this->decrement_stalls();
    if (!this->is_stalled()) {
      // TODO: Problem, we need a way to call the "extend lifetime" logic here, ideally without
      // implementing it twice
      this->check_constraints(rand_eng);
    }
    return false;
  }

  if (this->_direction == DNA::Direction::fwd) {
    return this->try_moving_to_next_bin();
  }
  return this->try_moving_to_prev_bin();
}

bool ExtrusionUnit::try_moving_to_next_bin() {
  //  fmt::fprintf(stderr, "Trying to move to next bin...");
  if (this->_bin == &this->_parent_lef.get_last_bin()) {
    this->set_stalls(UINT32_MAX);  // Stall forever
    return false;
  }
  //  fmt::fprintf(stderr, " Moving from bin #%lu", this->_bin->get_index());
  this->_bin->remove_extr_unit_binding(this);
  this->_bin = &this->_parent_lef.get_ptr_to_chr()->dna.get_next_bin(*this->_bin);
  this->_bin->add_extr_unit_binding(this);
  this->_blocking_barrier = nullptr;
  //  fmt::fprintf(stderr, " to bin #%lu!\n", this->_bin->get_index());
  return true;
}

bool ExtrusionUnit::try_moving_to_prev_bin() {
  if (this->_bin == &this->_parent_lef.get_first_bin()) {
    this->set_stalls(UINT32_MAX);  // Stall forever
    return false;
  }
  //  fmt::fprintf(stderr, " Moving from bin #%lu", this->_bin->get_index());
  this->_bin->remove_extr_unit_binding(this);
  this->_bin = &this->_parent_lef.get_ptr_to_chr()->dna.get_prev_bin(*this->_bin);
  this->_bin->add_extr_unit_binding(this);
  this->_blocking_barrier = nullptr;
  //  fmt::fprintf(stderr, " to bin #%lu!\n", this->_bin->get_index());
  return true;
}

bool ExtrusionUnit::hard_stall() const {
  return this->_blocking_barrier &&  // Not sure why clang-tidy is complaining here
         this->get_extr_direction() == this->_blocking_barrier->get_direction_of_block();
}

void ExtrusionUnit::set_stalls(uint32_t n) { this->_stalls_left = n; }
void ExtrusionUnit::increment_stalls(uint32_t n) {
  assert(UINT32_MAX - n >= this->_stalls_left);  // NOLINT
  this->_stalls_left += n;
}

void ExtrusionUnit::decrement_stalls(uint32_t n) {
  assert(this->_stalls_left >= n);  // NOLINT
  this->_stalls_left -= n;
}

void ExtrusionUnit::reset_stalls() { this->set_stalls(0); }

void ExtrusionUnit::unload() {
  this->reset_stalls();
  this->_bin->remove_extr_unit_binding(this);
  this->_bin = nullptr;
  this->_blocking_barrier = nullptr;
  this->_direction = DNA::Direction::none;
}

void ExtrusionUnit::bind(Chromosome* chr, uint32_t pos, DNA::Direction direction,
                         std::mt19937& rand_eng) {
  // TODO: We should also set a stall if another extr unit with the proper orientation is bound to
  // this bin
  assert(pos < chr->length());                                                   // NOLINT
  assert(direction == DNA::Direction::fwd || direction == DNA::Direction::rev);  // NOLINT
  this->_bin = chr->dna.get_ptr_to_bin_from_pos(pos);
  this->_bin->add_extr_unit_binding(this);
  this->_direction = direction;

  this->_blocking_barrier = nullptr;
  if (this->_bin->has_extr_barrier()) {
    auto& barriers = this->_bin->get_all_extr_barriers();
    if (direction == DNA::Direction::fwd) {
      // I think in this case ptr arithmetic is less error-prone than taking ptr to references/iters
      // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
      for (auto* b = &barriers.front(); b < &barriers.back(); ++b) {
        if (b->get_pos() >= pos) {
          this->_blocking_barrier = b;
          break;
        }
      }
    } else {
      // I think in this case ptr arithmetic is less error-prone than taking ptr to references/iters
      // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
      for (auto* b = &barriers.back(); b > &barriers.front(); --b) {
        if (b->get_pos() <= pos) {
          this->_blocking_barrier = b;
          break;
        }
      }
    }
    if (this->_blocking_barrier) {
      auto n_stalls = this->_blocking_barrier->generate_num_stalls(rand_eng);
      if (this->_blocking_barrier->get_direction_of_block() != this->get_extr_direction()) {
        n_stalls /= 2;
      }
      this->set_stalls(n_stalls);
      // fmt::fprintf(stderr, "LEF bound to a bin with a blocking extrusion barrier!\n");
      // NOLINTNEXTLINE
      assert(this->_blocking_barrier >= &this->_bin->get_all_extr_barriers().front() &&
             this->_blocking_barrier <= &this->_bin->get_all_extr_barriers().back());
    }
  }
}

bool ExtrusionUnit::is_bound() const { return this->_bin != nullptr; }

uint64_t ExtrusionUnit::check_constraints(std::mt19937& rand_eng) {
  if (this->is_stalled()) {
    return 0;
  }
  // This already applies the appropriate stall to colliding units.
  // Communicating this to Lef should not be necessary.
  this->check_for_extruder_collisions(rand_eng);
  // This already figures out if whether we should apply a "big" or "small" stall,
  // based on the direction of extrusion as well as the direction of the barrier
  // that is blocking the extrusion unit
  try {
    return this->check_for_extrusion_barrier(rand_eng);
  } catch (const std::runtime_error& err) {
    throw std::runtime_error(
        fmt::format("Exception caught while processing ExtrUnit bound at pos %lu! %s",
                    this->get_pos(), err.what()));
  }
}

uint32_t ExtrusionUnit::check_for_extruder_collisions(std::mt19937& rang_eng) {
  assert(!this->is_stalled());                       // NOLINT
  assert(this->_direction == DNA::Direction::fwd ||  // NOLINT
         this->_direction == DNA::Direction::rev);   // NOLINT
  assert(this->_bin);                                // NOLINT
  assert(this->_bin->get_n_extr_units() > 0);        // NOLINT
  uint32_t n_stalls = 0;
  // Avoid checking further if this is the only ExtrusionUnit bound to this->_bin
  if (this->_bin->get_n_extr_units() == 1) {
    return n_stalls;
  }

  for (auto& other : this->_bin->get_extr_units()) {
    assert(other->_direction == DNA::Direction::fwd ||  // NOLINT
           other->_direction == DNA::Direction::rev);   // NOLINT
    // Skip over extr. units belonging to the same LEF and apply a stall if this and other are
    // extruding in opposite directions
    if (&this->_parent_lef != &other->_parent_lef && this->_direction != other->_direction) {
      n_stalls = this->_n_stall_generator(rang_eng);
      this->increment_stalls(n_stalls);
      // In case other is not stalled, also apply the stall to that unit
      if (!other->is_stalled()) {
        other->increment_stalls(n_stalls);
      }
      break;
    }
  }
  return n_stalls;
}

uint64_t ExtrusionUnit::check_for_extrusion_barrier(std::mt19937& rang_eng) {
  uint32_t applied_stall = 0;
  if (this->get_extr_direction() == DNA::fwd) {
    this->_blocking_barrier = this->_bin->get_ptr_to_next_extr_barrier(this->_blocking_barrier);
  } else {
    this->_blocking_barrier = this->_bin->get_ptr_to_prev_extr_barrier(this->_blocking_barrier);
  }
  if (this->_blocking_barrier) {
    applied_stall = this->_blocking_barrier->generate_num_stalls(rang_eng);
    if (this->_blocking_barrier->get_direction_of_block() != this->get_extr_direction()) {
      // "Small" stall
      applied_stall /= 2;  // TODO: Make this tunable
    }
    this->increment_stalls(applied_stall);
  }
  return applied_stall;
}

Lef::Lef(uint32_t bin_size, uint32_t avg_processivity, double probability_of_extruder_bypass,
         double unloader_strength_coefficient)
    : _avg_processivity(avg_processivity),
      _probability_of_extr_unit_bypass(probability_of_extruder_bypass),
      _unloader_strength_coeff(unloader_strength_coefficient),
      _lifetime_generator(compute_prob_of_unloading(bin_size)),
      _left_unit(std::make_unique<ExtrusionUnit>(*this, probability_of_extruder_bypass)),
      _right_unit(std::make_unique<ExtrusionUnit>(*this, probability_of_extruder_bypass)) {}

bool Lef::is_bound() const {
  assert(this->_left_unit->is_bound() == this->_right_unit->is_bound());  // NOLINT
  return this->_left_unit->is_bound();
}

void Lef::randomly_bind_to_chr(Chromosome* chr, std::mt19937& rand_eng, bool register_contact) {
  std::uniform_int_distribution<uint32_t> pos(0, static_cast<uint32_t>(chr->length() - 1));
  this->bind_at_pos(chr, pos(rand_eng), rand_eng, register_contact);
}
void Lef::bind_at_pos(Chromosome* chr, uint32_t pos, std::mt19937& rand_eng,
                      bool register_contact) {
  assert(!this->_left_unit->is_bound() && !this->_right_unit->is_bound());  // NOLINT
  this->_chr = chr;
  this->_binding_pos = pos;
  this->_lifetime = this->_lifetime_generator(rand_eng);
  std::discrete_distribution<int8_t> bin_idx_offset_gen({1, 1, 0, 1, 1});
  // We assume that the left unit always travels towards the 5', while the right unit goes to the 3'
  this->_left_unit->bind(this->_chr, pos, DNA::Direction::rev, rand_eng);
  this->_right_unit->bind(this->_chr, pos, DNA::Direction::fwd, rand_eng);
  if (int8_t bin_idx_offset = bin_idx_offset_gen(rand_eng) - 2; bin_idx_offset < 0) {
    for (; bin_idx_offset < 0; ++bin_idx_offset) {
      this->_left_unit->try_moving_to_prev_bin();
    }
  } else {
    assert(bin_idx_offset != 0);  // NOLINT
    for (; bin_idx_offset < 2; ++bin_idx_offset) {
      this->_right_unit->try_moving_to_next_bin();
    }
  }
  if (register_contact) {
    const auto n = this->_chr->dna.get_ptr_to_bin_from_pos(pos)->get_index();
    this->_chr->contacts.increment(n, n);
  }
}

std::string_view Lef::get_chr_name() const { return this->_chr->name; }

std::pair<DNA::Bin*, DNA::Bin*> Lef::get_ptr_to_bins() {
  return std::make_pair(this->_left_unit->_bin, this->_left_unit->_bin);
}

void Lef::unload() {
  this->_chr = nullptr;
  this->_binding_pos = UINT64_MAX;
  this->_left_unit->unload();
  this->_right_unit->unload();
  this->_lifetime = 0;  // Probably unnecessary
}

uint32_t Lef::extrude(std::mt19937& rand_eng) {
  if (this->_lifetime-- > 0) {
    const auto bp_extruded =  // NOLINTNEXTLINE(readability-implicit-bool-conversion)
        (this->_left_unit->try_extrude(rand_eng) + this->_right_unit->try_extrude(rand_eng)) *
        this->_right_unit->_bin->size();
    this->_tot_bp_extruded += bp_extruded;
    return static_cast<uint32_t>(bp_extruded);
  }
  this->unload();
  return 0;
}

void Lef::register_contact() {
  // We don't register contacts for the first and last bins in a chromosome
  if (this->_left_unit->get_bin_index() == 0 &&
      this->_right_unit->get_bin_index() == this->get_last_bin().get_index()) {
    this->_chr->contacts.increment(this->_left_unit->get_bin_index(),
                                   this->_right_unit->get_bin_index());
  }
}

void Lef::check_constraints(std::mt19937& rang_eng) {
  if (this->_left_unit->is_stalled() && this->_right_unit->is_stalled()) {
    return;
  }
  this->_left_unit->check_constraints(rang_eng);
  this->_right_unit->check_constraints(rang_eng);
  if (this->_unloader_strength_coeff > 0 && this->_left_unit->hard_stall() &&
      this->_right_unit->hard_stall()) {
    const auto n_of_stalls_due_to_barriers =
        std::lround(this->_unloader_strength_coeff *
                    std::min(this->_left_unit->_stalls_left, this->_right_unit->_stalls_left));
    assert(UINT32_MAX - n_of_stalls_due_to_barriers > this->_lifetime);  // NOLINT
    this->_lifetime += n_of_stalls_due_to_barriers;
  }
}

uint32_t Lef::get_loop_size() const {
  assert(this->_right_unit->get_pos() >= this->_left_unit->get_pos());  // NOLINT
  return this->_right_unit->get_pos() - this->_left_unit->get_pos();
}

uint32_t Lef::get_avg_processivity() const { return this->_avg_processivity; }

bool Lef::try_rebind(Chromosome& chr, std::mt19937& rand_eng, double prob_of_rebinding,
                     bool register_contact) {
  assert(!this->is_bound());                                 // NOLINT
  assert(prob_of_rebinding >= 0 && prob_of_rebinding <= 1);  // NOLINT
  std::uniform_real_distribution<> d1(0.0, 1.0);
  if (d1(rand_eng) >= 1 - prob_of_rebinding) {
    //    fmt::fprintf(stderr, "Trying to rebind...\n");
    this->randomly_bind_to_chr(&chr, rand_eng, register_contact);
    return true;
  }
  return false;
}

const DNA::Bin& Lef::get_first_bin() const { return this->_chr->dna.get_first_bin(); }
const DNA::Bin& Lef::get_last_bin() const { return this->_chr->dna.get_last_bin(); }

Chromosome* Lef::get_ptr_to_chr() {
  assert(this->_chr);  // NOLINT
  return this->_chr;
}

std::pair<uint32_t, uint32_t> Lef::get_pos() const {
  return std::make_pair(this->_left_unit->get_pos(), this->_right_unit->get_pos());
}

double Lef::get_probability_of_extr_unit_bypass() const {
  return this->_probability_of_extr_unit_bypass;
}

uint64_t Lef::get_tot_bp_extruded() const { return this->_tot_bp_extruded; }

void Lef::reset_tot_bp_extruded() { this->_tot_bp_extruded = 0; }

double Lef::compute_prob_of_unloading(uint32_t bin_size, uint8_t n_of_active_extr_units) const {
  return n_of_active_extr_units / (static_cast<double>(this->get_avg_processivity()) / bin_size);
}

}  // namespace modle