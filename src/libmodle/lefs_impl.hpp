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
         double unloader_strength_coefficient)
    : _avg_lifetime(avg_lef_lifetime),
      _probability_of_extr_unit_bypass(probability_of_extruder_bypass),
      _unloader_strength_coeff(unloader_strength_coefficient),
      _lifetime_generator(compute_prob_of_unloading(bin_size)),
      _left_unit(std::make_unique<ExtrusionUnit>(*this, probability_of_extruder_bypass)),
      _right_unit(std::make_unique<ExtrusionUnit>(*this, probability_of_extruder_bypass)) {}

bool Lef::is_bound() const {
  assert(this->_left_unit->is_bound() == this->_right_unit->is_bound());  // NOLINT
  return this->_left_unit->is_bound();
}

void Lef::bind_chr_at_random_pos(Chromosome* chr, modle::PRNG& rand_eng, bool register_contact) {
  std::uniform_int_distribution<uint32_t> pos(0,
                                              static_cast<uint32_t>(chr->simulated_length() - 1));
  this->bind_at_pos(chr, pos(rand_eng), rand_eng, register_contact);
}

void Lef::assign_to_chr(Chromosome* chr) { this->_chr = chr; }

void Lef::bind_at_pos(Chromosome* chr, uint32_t pos, modle::PRNG& rand_eng, bool register_contact) {
  assert(!this->_left_unit->is_bound() && !this->_right_unit->is_bound());  // NOLINT
  this->_chr = chr;
  this->_binding_pos = pos;
  this->_lifetime = this->_lifetime_generator(rand_eng);
  const auto pos_offset =
      this->_bin_idx_offset_generator(rand_eng) * static_cast<int64_t>(this->_chr->get_bin_size());

  // It would be nice to remove some static casts here, but the first static_cast to int64_t is
  // required to deal with the possibility that pos - offset overflows
  const auto pos1 = static_cast<uint32_t>(
      std::clamp(pos_offset < 0 ? static_cast<int64_t>(pos) + pos_offset : pos, 0L,
                 static_cast<int64_t>(this->_chr->simulated_length() - 1)));
  const auto pos2 = static_cast<uint32_t>(
      std::clamp(pos_offset > 0 ? static_cast<int64_t>(pos) + pos_offset : pos, 0L,
                 static_cast<int64_t>(this->_chr->simulated_length() - 1)));
  // We assume that the left unit always travels towards the 5', while the right unit goes to the 3'
  this->_left_unit->bind(this->_chr, pos1, dna::Direction::rev, rand_eng);
  this->_right_unit->bind(this->_chr, pos2, dna::Direction::fwd, rand_eng);
  if (register_contact) {
    const auto n = this->_chr->dna.get_ptr_to_bin(pos)->get_index();
    this->_chr->contacts.increment(n, n);
  }
}

std::string_view Lef::get_chr_name() const { return this->_chr->name; }

std::pair<DNA::Bin*, DNA::Bin*> Lef::get_ptr_to_bins() {
  return {&this->_left_unit->get_bin(), &this->_right_unit->get_bin()};
}

void Lef::unload() {
  this->_binding_pos = std::numeric_limits<decltype(this->_binding_pos)>::max();
  this->_left_unit->unload();
  this->_right_unit->unload();
  this->_lifetime = 0;  // Probably unnecessary
}

uint32_t Lef::extrude(modle::PRNG& rand_eng) {
  if (this->_lifetime-- > 0) {
    const auto bp_extruded =
        static_cast<uint64_t>(
            this->_left_unit->try_extrude(  // NOLINT(readability-implicit-bool-conversion)
                rand_eng) +
            this->_right_unit->try_extrude(  // NOLINT(readability-implicit-bool-conversion)
                rand_eng)) *
        this->get_bin_size();

    this->_tot_bp_extruded += bp_extruded;
    return static_cast<uint32_t>(bp_extruded);
  }
  this->unload();
  return 0;
}

void Lef::register_contact() {
  // We don't register contacts for the first and last bins in a chromosome
  if (this->_left_unit->get_bin_index() != 0 &&
      this->_right_unit->get_bin_index() != this->get_last_bin().get_index()) {
    this->_chr->contacts.increment(this->_left_unit->get_bin_index(),
                                   this->_right_unit->get_bin_index());
  }
}

void Lef::check_constraints(modle::PRNG& rang_eng) {
  if (this->_left_unit->is_stalled() && this->_right_unit->is_stalled()) {
    return;
  }
  this->_left_unit->check_constraints(rang_eng);
  this->_right_unit->check_constraints(rang_eng);
  if (this->_unloader_strength_coeff > 0 && this->_left_unit->hard_stall() &&
      this->_right_unit->hard_stall()) {
    const auto n_of_stalls_due_to_barriers = static_cast<uint32_t>(
        std::lround(this->_unloader_strength_coeff *
                    std::min(this->_left_unit->_stalls_left, this->_right_unit->_stalls_left)));
    assert(std::numeric_limits<decltype(this->_lifetime)>::max() - n_of_stalls_due_to_barriers >
           this->_lifetime);  // NOLINT
    this->_lifetime += n_of_stalls_due_to_barriers;
  }
}

std::size_t Lef::get_loop_size() const {
  assert(this->_right_unit->get_pos() >= this->_left_unit->get_pos());  // NOLINT
  return this->_right_unit->get_pos() - this->_left_unit->get_pos();
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
  return std::make_pair(this->_left_unit->get_pos(), this->_right_unit->get_pos());
}

double Lef::get_probability_of_extr_unit_bypass() const {
  return this->_probability_of_extr_unit_bypass;
}

std::size_t Lef::get_bin_size() const { return this->_chr->get_bin_size(); }
std::size_t Lef::get_nbins() const { return this->_chr->get_nbins(); }

uint64_t Lef::get_tot_bp_extruded() const { return this->_tot_bp_extruded; }

void Lef::reset_tot_bp_extruded() { this->_tot_bp_extruded = 0; }

double Lef::compute_prob_of_unloading(uint32_t bin_size, uint8_t n_of_active_extr_units) const {
  return n_of_active_extr_units / (static_cast<double>(this->get_avg_lifetime()) / bin_size);
}

ExtrusionUnit::ExtrusionUnit(Lef& lef, double prob_of_extr_unit_bypass)
    : _parent_lef(lef), _n_stall_generator(prob_of_extr_unit_bypass) {}

std::size_t ExtrusionUnit::get_pos() const {
  return ((this->_bin_idx + (this->_bin_idx + 1)) * this->get_bin_size()) / 2;
}

dna::Direction ExtrusionUnit::get_extr_direction() const { return this->_direction; }

bool ExtrusionUnit::is_stalled() const { return this->_stalls_left > 0; }

double ExtrusionUnit::get_prob_of_extr_unit_bypass() const {
  return this->_parent_lef.get_probability_of_extr_unit_bypass();
}

std::size_t ExtrusionUnit::get_bin_index() const { return this->_bin_idx; }

std::size_t ExtrusionUnit::get_bin_size() const { return this->_parent_lef.get_bin_size(); }

DNA::Bin& ExtrusionUnit::get_bin() {
  assert(this->_bin_idx != std::numeric_limits<decltype(this->_bin_idx)>::max());  // NOLINT
  return (*this->_dna)[this->_bin_idx];
}

const DNA::Bin& ExtrusionUnit::get_bin() const {
  assert(this->_bin_idx != std::numeric_limits<decltype(this->_bin_idx)>::max());  // NOLINT
  return (*this->_dna)[this->_bin_idx];
}

bool ExtrusionUnit::try_extrude(modle::PRNG& rand_eng) {
  assert(this->_direction == dna::Direction::fwd ||  // NOLINT
         this->_direction == dna::Direction::rev);   // NOLINT
  if (this->is_stalled()) {
    this->decrement_stalls();
    if (!this->is_stalled()) {
      // TODO: Problem, we need a way to call the "extend lifetime" logic here, ideally without
      // implementing it twice
      this->check_constraints(rand_eng);
    }
    return false;
  }

  if (this->_direction == dna::Direction::fwd) {
    return this->try_moving_to_next_bin();
  }
  return this->try_moving_to_prev_bin();
}

bool ExtrusionUnit::try_moving_to_next_bin() {
  if (this->_bin_idx == this->_dna->get_n_bins() - 1) {
    this->set_stalls(std::numeric_limits<decltype(this->_stalls_left)>::max());  // Stall forever
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
    this->set_stalls(std::numeric_limits<decltype(this->_stalls_left)>::max());  // Stall forever
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
  this->get_bin().unregister_extr_unit_binding(this);
  this->_bin_idx = std::numeric_limits<decltype(this->_bin_idx)>::max();
  this->_blocking_barrier = nullptr;
  this->_direction = dna::Direction::none;
}

void ExtrusionUnit::bind(Chromosome* chr, uint32_t pos, dna::Direction direction,
                         modle::PRNG& rand_eng) {
  // TODO: We should also set a stall if another extr unit with the proper orientation is bound to
  // this bin
  assert(pos < chr->end);                                                        // NOLINT
  assert(direction == dna::Direction::fwd || direction == dna::Direction::rev);  // NOLINT
  this->_dna = &chr->dna;
  this->_bin_idx = this->_dna->get_bin_idx(pos);
  this->get_bin().register_extr_unit_binding(this);
  this->_direction = direction;

  if (this->_direction == dna::fwd) {
    this->_blocking_barrier = this->get_bin().get_ptr_to_next_extr_barrier(pos);
  } else {
    this->_blocking_barrier = this->get_bin().get_ptr_to_prev_extr_barrier(pos);
  }
  if (this->_blocking_barrier) {
    auto n_stalls = this->_blocking_barrier->generate_num_stalls(rand_eng);
    if (this->_blocking_barrier->get_direction_of_block() != this->get_extr_direction()) {
      n_stalls /= 2;
    }
    this->set_stalls(n_stalls);
    // NOLINTNEXTLINE
    // assert(this->_blocking_barrier >= &this->_bin->get_all_extr_barriers().front() &&
    //       this->_blocking_barrier <= &this->_bin->get_all_extr_barriers().back());
  }
}

bool ExtrusionUnit::is_bound() const {
  // NOLINTNEXTLINE
  assert(this->_bin_idx == std::numeric_limits<decltype(this->_bin_idx)>::max() ||
         this->_bin_idx < this->_dna->get_n_bins());
  return this->_bin_idx != std::numeric_limits<decltype(this->_bin_idx)>::max();
}

uint64_t ExtrusionUnit::check_constraints(modle::PRNG& rand_eng) {
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
        fmt::format("Exception caught while processing ExtrUnit bound at pos {}! {}",
                    this->get_pos(), err.what()));
  }
}

uint32_t ExtrusionUnit::check_for_extruder_collisions(modle::PRNG& rang_eng) {
  assert(!this->is_stalled());                        // NOLINT
  assert(this->_direction == dna::Direction::fwd ||   // NOLINT
         this->_direction == dna::Direction::rev);    // NOLINT
  assert(this->_bin_idx < this->_dna->get_n_bins());  // NOLINT
  assert(this->get_bin().get_n_extr_units() > 0);     // NOLINT
  uint32_t n_stalls = 0;
  // Avoid checking further if this is the only ExtrusionUnit bound to this->_bin
  if (this->get_bin().get_n_extr_units() == 1) {
    return n_stalls;
  }

  for (auto& other : this->get_bin().get_extr_units()) {
    assert(other->_direction == dna::Direction::fwd ||  // NOLINT
           other->_direction == dna::Direction::rev);   // NOLINT
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

uint64_t ExtrusionUnit::check_for_extrusion_barrier(modle::PRNG& rang_eng) {
  uint32_t applied_stall = 0;
  if (this->get_extr_direction() == dna::fwd) {
    this->_blocking_barrier = this->get_bin().get_ptr_to_next_extr_barrier(this->_blocking_barrier);
  } else {
    this->_blocking_barrier = this->get_bin().get_ptr_to_prev_extr_barrier(this->_blocking_barrier);
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

}  // namespace modle
