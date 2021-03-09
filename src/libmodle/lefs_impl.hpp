#pragma once

#include <fmt/format.h>  // for format
#ifdef USE_XOSHIRO
#include <XoshiroCpp.hpp>
#endif

#include <algorithm>  // for min
#include <cassert>
#include <cmath>      // for lround
#include <stdexcept>  // for runtime_error
#include <string>
#include <vector>

#include "modle/contacts.hpp"      // for ContactMatrix
#include "modle/extr_barrier.hpp"  // for ExtrusionBarrier

namespace modle {

Lef::Lef(uint32_t bin_size, uint32_t avg_lef_lifetime, double probability_of_extruder_bypass,
         double hard_stall_multiplier, double soft_stall_multiplier)
    : _avg_lifetime(avg_lef_lifetime),
      _probability_of_extr_unit_bypass(probability_of_extruder_bypass),
      _hard_stall_multiplier(hard_stall_multiplier),
      _soft_stall_multiplier(soft_stall_multiplier),
      _lifetime_generator(compute_prob_of_unloading(bin_size)),
      _left_unit(nullptr, probability_of_extruder_bypass),
      _right_unit(nullptr, probability_of_extruder_bypass) {}

bool Lef::is_bound() const {
  assert(this->_left_unit.is_bound() == this->_right_unit.is_bound());  // NOLINT
  return this->_left_unit.is_bound();
}

void Lef::bind_chr_at_random_pos(Chromosome* chr, modle::PRNG& rand_eng, bool register_contact) {
  std::uniform_int_distribution<uint32_t> pos(0,
                                              static_cast<uint32_t>(chr->simulated_length() - 1));
  this->bind_at_pos(chr, pos(rand_eng), rand_eng, register_contact);
}

void Lef::assign_to_chr(Chromosome* chr) { this->_chr = chr; }

void Lef::bind_at_pos(Chromosome* chr, uint32_t pos, modle::PRNG& rand_eng, bool register_contact) {
  assert(!this->_left_unit.is_bound() && !this->_right_unit.is_bound());  // NOLINT
  this->_chr = chr;
  this->_binding_pos = pos;
  // TODO consider whether it is ok to assign a lifetime of 0
  this->_lifetime = this->_lifetime_generator(rand_eng);

  // We assume that the left unit always travels towards the 5', while the right unit goes to the 3'
  this->_left_unit.bind(this->_chr, pos, dna::Direction::rev);
  this->_right_unit.bind(this->_chr, pos, dna::Direction::fwd);
  if (register_contact) {
    const auto n = this->_chr->dna.get_ptr_to_bin(pos)->get_index();
    this->_chr->contacts.increment(n, n);
  }
  // Checking constraints here is necessary to avoid skipping another LEF or extr. barrier
  // already bound to the bin corresponding to pos, as later calls to
  // try_extrude_and_check_constraints first extrudes, and then checks the constraints
  this->check_constraints(rand_eng);

  // We apply an offset to avoid the "artificial" checker board pattern that arises because each
  // call to Lef::try_extrude() always extrudes an even (i.e. 2) number of bins: one left and one
  // right. Here we are calling try_extrude_and_check_constraints, so that we can deal with the
  // rare, but possible LEF-LEF/BAR collisions that occur immediately after binding
  if (const auto idx_offset = this->_bin_idx_offset_generator(rand_eng); idx_offset > 0) {
    for (auto i = 0; i < idx_offset; ++i) {
      this->_right_unit.try_extrude_and_check_constraints(rand_eng);
    }
  } else if (idx_offset < 0) {
    for (auto i = idx_offset; i < 0; ++i) {
      this->_left_unit.try_extrude_and_check_constraints(rand_eng);
    }
  }

  if (this->hard_stall()) {
    this->apply_hard_stall_and_extend_lifetime();
  }
}

std::string_view Lef::get_chr_name() const { return this->_chr->name; }

std::pair<DNA::Bin*, DNA::Bin*> Lef::get_ptr_to_bins() {
  return {&this->_left_unit.get_bin(), &this->_right_unit.get_bin()};
}

void Lef::unload() {
  this->_binding_pos = std::numeric_limits<decltype(this->_binding_pos)>::max();
  this->_left_unit.unload();
  this->_right_unit.unload();
  this->_lifetime = 0;  // Probably unnecessary
}

uint32_t Lef::try_extrude() {
  if (this->_lifetime-- > 0) {
    const auto bins_extruded =
        static_cast<uint64_t>(  // NOLINTNEXTLINE(readability-implicit-bool-conversion)
            this->_left_unit.try_extrude() + this->_right_unit.try_extrude());

    this->_tot_bp_extruded += bins_extruded * this->get_bin_size();
    return static_cast<uint32_t>(bins_extruded);
  }
  this->unload();
  return 0;
}

void Lef::register_contact() {
  // We don't register contacts for the first and last bins in a chromosome
  if (this->_left_unit.get_bin_index() != 0 &&
      this->_right_unit.get_bin_index() != this->get_last_bin().get_index()) {
    this->_chr->contacts.increment(this->_left_unit.get_bin_index(),
                                   this->_right_unit.get_bin_index());
  }
}

void Lef::check_constraints(modle::PRNG& rang_eng) {
  this->_left_unit.check_for_lef_lef_collisions(rang_eng);
  this->_right_unit.check_for_lef_lef_collisions(rang_eng);
  if (this->_left_unit._nstalls_lef_bar > 0 && this->_right_unit._nstalls_lef_bar > 0) {
    return;
  }

  if (this->_left_unit._nstalls_lef_bar == 0) {
    this->_left_unit.check_for_lef_bar_collision(rang_eng);
  }

  if (this->_right_unit._nstalls_lef_bar == 0) {
    this->_right_unit.check_for_lef_bar_collision(rang_eng);
  }

  if (this->hard_stall()) {
    this->apply_hard_stall_and_extend_lifetime();
  }
}

std::size_t Lef::get_loop_size() const {
  assert(this->_right_unit.get_pos() >= this->_left_unit.get_pos());  // NOLINT
  return this->_right_unit.get_pos() - this->_left_unit.get_pos();
}

std::size_t Lef::get_avg_lifetime() const { return static_cast<std::size_t>(this->_avg_lifetime); }

bool Lef::try_rebind(modle::PRNG& rand_eng, double prob_of_rebinding, bool register_contact) {
  assert(!this->is_bound());                                 // NOLINT
  assert(prob_of_rebinding >= 0 && prob_of_rebinding <= 1);  // NOLINT
  if (prob_of_rebinding == 1.0 || this->_prob_of_rebind_generator(rand_eng) <= prob_of_rebinding) {
    this->bind_chr_at_random_pos(this->_chr, rand_eng, register_contact);
    return true;
  }
  return false;
}

bool Lef::try_rebind(modle::PRNG& rand_eng) { return this->try_rebind(rand_eng, 1, false); }

std::size_t Lef::bind_at_random_pos(modle::PRNG& rand_eng, bool register_contact) {
  this->try_rebind(rand_eng, 1, register_contact);
  return this->get_pos().first;
}

const DNA::Bin& Lef::get_first_bin() const { return this->_chr->dna.get_first_bin(); }
const DNA::Bin& Lef::get_last_bin() const { return this->_chr->dna.get_last_bin(); }

Chromosome* Lef::get_ptr_to_chr() {
  assert(this->_chr);  // NOLINT
  return this->_chr;
}

std::pair<std::size_t, std::size_t> Lef::get_pos() const {
  return std::make_pair(this->_left_unit.get_pos(), this->_right_unit.get_pos());
}

double Lef::get_probability_of_extr_unit_bypass() const {
  return this->_probability_of_extr_unit_bypass;
}

std::size_t Lef::get_bin_size() const { return this->_chr->get_bin_size(); }
std::size_t Lef::get_nbins() const { return this->_chr->get_nbins(); }

uint64_t Lef::get_tot_bp_extruded() const { return this->_tot_bp_extruded; }

double Lef::get_hard_stall_multiplier() const { return this->_hard_stall_multiplier; }

double Lef::get_soft_stall_multiplier() const { return this->_soft_stall_multiplier; }

void Lef::reset_tot_bp_extruded() { this->_tot_bp_extruded = 0; }

double Lef::compute_prob_of_unloading(uint32_t bin_size, uint8_t n_of_active_extr_units) const {
  return n_of_active_extr_units / (static_cast<double>(this->get_avg_lifetime()) / bin_size);
}

bool Lef::hard_stall() const {
  return this->_left_unit.hard_stall() && this->_right_unit.hard_stall();
}

uint32_t Lef::apply_hard_stall_and_extend_lifetime(bool allow_lifetime_extension) {
  assert(this->hard_stall());  // NOLINT
  const auto nstalls = static_cast<uint32_t>(
      std::round(this->_hard_stall_multiplier *
                 std::min(this->_left_unit._nstalls_lef_bar, this->_right_unit._nstalls_lef_bar)));

  this->_left_unit.increment_lef_bar_stalls(nstalls);
  this->_right_unit.increment_lef_bar_stalls(nstalls);

  if (allow_lifetime_extension &&
      std::numeric_limits<decltype(this->_lifetime)>::max() - nstalls > this->_lifetime) {
    this->_lifetime += nstalls;
  }

  return nstalls;
}

void Lef::finalize_extrusion_unit_construction() {
  this->_left_unit._parent_lef = this;
  this->_right_unit._parent_lef = this;
}

ExtrusionUnit::ExtrusionUnit(Lef* lef, double prob_of_extr_unit_bypass)
    : _parent_lef(lef), _n_lef_lef_stall_generator(prob_of_extr_unit_bypass) {}

std::size_t ExtrusionUnit::get_pos() const {
  return ((this->_bin_idx + (this->_bin_idx + 1)) * this->get_bin_size()) / 2;
}

dna::Direction ExtrusionUnit::get_extr_direction() const { return this->_direction; }

bool ExtrusionUnit::is_stalled() const {
  return (this->_nstalls_lef_lef + this->_nstalls_lef_bar) > 0;
}

double ExtrusionUnit::get_prob_of_extr_unit_bypass() const {
  return this->_parent_lef->get_probability_of_extr_unit_bypass();
}

std::size_t ExtrusionUnit::get_bin_index() const { return this->_bin_idx; }

std::size_t ExtrusionUnit::get_bin_size() const { return this->_parent_lef->get_bin_size(); }

DNA::Bin& ExtrusionUnit::get_bin() {
  assert(this->_bin_idx < this->_dna->size());  // NOLINT
  return (*this->_dna)[this->_bin_idx];
}

const DNA::Bin& ExtrusionUnit::get_bin() const {
  assert(this->_bin_idx < this->_dna->size());  // NOLINT
  return (*this->_dna)[this->_bin_idx];
}

bool ExtrusionUnit::try_extrude() {
  assert(this->_direction == dna::Direction::fwd ||  // NOLINT
         this->_direction == dna::Direction::rev);   // NOLINT

  if (this->is_stalled()) {
    this->decrement_stalls();
    return false;
  }

  if (this->_direction == dna::Direction::fwd) {
    return this->try_moving_to_next_bin();
  }
  return this->try_moving_to_prev_bin();
}

bool ExtrusionUnit::try_extrude_and_check_constraints(modle::PRNG& rand_eng) {
  const auto extr_successful = this->try_extrude();
  if (extr_successful) {
    this->check_constraints(rand_eng);
  }

  return extr_successful;
}

bool ExtrusionUnit::try_moving_to_next_bin() {
  if (this->_bin_idx == this->_dna->get_n_bins() - 1) {
    this->set_lef_bar_stalls(  // Stall until lef rebind
        std::numeric_limits<decltype(this->_nstalls_lef_bar)>::max());
    return false;
  }
  this->get_bin().unregister_extr_unit_binding(this);
  ++this->_bin_idx;
  this->get_bin().register_extr_unit_binding(this);
  this->_blocking_barrier = nullptr;
  return true;
}

bool ExtrusionUnit::try_moving_to_prev_bin() {
  if (this->_bin_idx == 0) {
    this->set_lef_bar_stalls(  // Stall until lef rebind
        std::numeric_limits<decltype(this->_nstalls_lef_bar)>::max());
    return false;
  }
  this->get_bin().unregister_extr_unit_binding(this);
  --this->_bin_idx;
  this->get_bin().register_extr_unit_binding(this);
  this->_blocking_barrier = nullptr;
  return true;
}

bool ExtrusionUnit::hard_stall() const {
  return this->_blocking_barrier &&  // Not sure why clang-tidy is complaining here
         this->get_extr_direction() == this->_blocking_barrier->get_direction_of_block();
}

void ExtrusionUnit::set_lef_bar_stalls(uint32_t n) { this->_nstalls_lef_bar = n; }
void ExtrusionUnit::set_lef_lef_stalls(uint32_t n) { this->_nstalls_lef_lef = n; }
void ExtrusionUnit::increment_lef_lef_stalls(uint32_t n) {  // NOLINTNEXTLINE
  if (std::numeric_limits<decltype(this->_nstalls_lef_lef)>::max() - n >= this->_nstalls_lef_lef) {
    this->_nstalls_lef_lef += n;
  }
}

void ExtrusionUnit::increment_lef_bar_stalls(uint32_t n) {  // NOLINTNEXTLINE
  if (std::numeric_limits<decltype(this->_nstalls_lef_bar)>::max() - n >= this->_nstalls_lef_bar) {
    this->_nstalls_lef_bar += n;
  }
}

void ExtrusionUnit::decrement_stalls(uint32_t n) {
  assert(this->is_stalled());
  // Stalls due to LEF-LEF collisions have precedence over LEF-barrier stalls
  if (this->_nstalls_lef_lef > 0) {
    this->decrement_lef_lef_stalls(n);
  } else {
    this->decrement_lef_bar_stalls(n);
  }
}

void ExtrusionUnit::decrement_lef_bar_stalls(uint32_t n) {  // NOLINTNEXTLINE
  assert(this->_nstalls_lef_bar >=
         std::numeric_limits<decltype(this->_nstalls_lef_bar)>::min() + n);
  this->_nstalls_lef_bar -= n;
}

void ExtrusionUnit::decrement_lef_lef_stalls(uint32_t n) {  // NOLINTNEXTLINE
  assert(this->_nstalls_lef_lef >=
         std::numeric_limits<decltype(this->_nstalls_lef_lef)>::min() + n);
  this->_nstalls_lef_lef -= n;
}

void ExtrusionUnit::reset_lef_bar_stalls() { this->set_lef_bar_stalls(0); }

void ExtrusionUnit::reset_lef_lef_stalls() { this->set_lef_lef_stalls(0); }

void ExtrusionUnit::reset_stalls() {
  this->reset_lef_bar_stalls();
  this->reset_lef_lef_stalls();
}

void ExtrusionUnit::unload() {
  this->reset_stalls();
  this->get_bin().unregister_extr_unit_binding(this);
  this->_bin_idx = std::numeric_limits<decltype(this->_bin_idx)>::max();
  this->_blocking_barrier = nullptr;
  // this->_direction = dna::Direction::none;
}

void ExtrusionUnit::bind(Chromosome* chr, uint32_t pos, dna::Direction direction) {
  assert(pos < chr->end);                                                        // NOLINT
  assert(direction == dna::Direction::fwd || direction == dna::Direction::rev);  // NOLINT
  this->_dna = &chr->dna;
  this->_bin_idx = this->_dna->get_bin_idx(pos);
  this->get_bin().register_extr_unit_binding(this);
  this->_direction = direction;

  // Note: The callee is expected to check if all constraints are satisfied before moving to the
  // next/previous bin
}

bool ExtrusionUnit::is_bound() const {
  // NOLINTNEXTLINE
  assert(this->_bin_idx == std::numeric_limits<decltype(this->_bin_idx)>::max() ||
         this->_bin_idx < this->_dna->get_n_bins());
  return this->_bin_idx != std::numeric_limits<decltype(this->_bin_idx)>::max();
}

uint64_t ExtrusionUnit::check_constraints(modle::PRNG& rand_eng) {
  uint64_t nstalls = 0;
  if (this->is_stalled()) {
    if (this->get_bin().get_n_extr_units() == 1) {
      this->_nstalls_lef_lef = 0;
    }
    return nstalls;
  }
  // This already applies the appropriate stall to colliding units.
  // Communicating this to Lef should not be necessary.
  nstalls += this->check_for_lef_lef_collisions(rand_eng);
  // This already figures out if whether we should apply a "big" or "small" stall,
  // based on the direction of extrusion as well as the direction of the barrier
  // that is blocking the extrusion unit
#ifdef NDEBUG
  return nstalls + this->check_for_lef_bar_collision(rand_eng);
#else
  try {
    return nstalls + this->check_for_lef_bar_collision(rand_eng);
  } catch (const std::runtime_error& err) {
    throw std::runtime_error(
        fmt::format("Exception caught while processing ExtrUnit bound at pos {}! {}",
                    this->get_pos(), err.what()));
  }
#endif
}

uint32_t ExtrusionUnit::check_for_lef_lef_collisions(modle::PRNG& rang_eng) {
  assert(this->_direction == dna::Direction::fwd ||   // NOLINT
         this->_direction == dna::Direction::rev);    // NOLINT
  assert(this->_bin_idx < this->_dna->get_n_bins());  // NOLINT
  assert(this->get_bin().get_n_extr_units() > 0);     // NOLINT
  uint32_t nstalls = 0;

  // If there's a single extrusion unit bound to the current bin, reset the number of lef_lef stalls
  // and return immediately
  if (this->get_bin().get_n_extr_units() == 1) {
    this->_nstalls_lef_lef = 0;
    return nstalls;
  }

  // Avoid checking further if this lef has already been stalled due to lef_lef collisions, or if
  // the prob. of lef bypass is 1.0 (i.e. lefs are effectively "transparent" to each other)
  if (this->_nstalls_lef_lef > 0 || this->_n_lef_lef_stall_generator.p() == 1.0) {
    return nstalls;
  }

  for (auto& other : this->get_bin().get_extr_units()) {
    assert(other->_direction == dna::Direction::fwd ||  // NOLINT
           other->_direction == dna::Direction::rev);   // NOLINT
    // Skip over extr. units belonging to the same LEF and apply a stall if this and other are
    // extruding in opposite directions
    if (this->_parent_lef != other->_parent_lef && this->_direction != other->_direction) {
      if (this->_n_lef_lef_stall_generator.p() == 0.0) {
        nstalls = std::numeric_limits<decltype(this->_nstalls_lef_lef)>::max();
        this->set_lef_lef_stalls(nstalls);
        other->set_lef_lef_stalls(nstalls);
      } else {
        nstalls = this->_n_lef_lef_stall_generator(rang_eng);
        this->increment_lef_lef_stalls(nstalls);
        other->increment_lef_lef_stalls(nstalls);
      }
    }
  }
  return nstalls;
}

uint32_t ExtrusionUnit::check_for_lef_bar_collision(modle::PRNG& rang_eng) {
  uint32_t nstalls = 0;
  if (this->get_extr_direction() == dna::fwd) {
    this->_blocking_barrier = this->get_bin().get_ptr_to_next_extr_barrier(this->_blocking_barrier);
  } else {
    this->_blocking_barrier = this->get_bin().get_ptr_to_prev_extr_barrier(this->_blocking_barrier);
  }
  if (this->_blocking_barrier) {
    nstalls = this->_blocking_barrier->generate_nstalls(rang_eng);
    if (this->_blocking_barrier->get_direction_of_block() != this->get_extr_direction()) {
      // Soft stall
      nstalls = static_cast<uint32_t>(std::round(this->_parent_lef->get_soft_stall_multiplier() *
                                                 static_cast<double>(nstalls)));
    }
    this->increment_lef_bar_stalls(nstalls);
  }
  return nstalls;
}

}  // namespace modle
