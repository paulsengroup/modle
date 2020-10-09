#include "modle/lefs.hpp"

#include "absl/strings/str_format.h"

namespace modle {

ExtrusionUnit::ExtrusionUnit(Lef& lef, double prob_of_extr_unit_bypass)
    : _parent_lef(lef), _n_of_stalls_gen(prob_of_extr_unit_bypass) {}

uint32_t ExtrusionUnit::get_pos() const {
  assert(this->_bin);
  return (this->_bin->get_end() + this->_bin->get_start()) / 2;
}

DNA::Direction ExtrusionUnit::get_extr_direction() const { return this->_direction; }

bool ExtrusionUnit::is_stalled() const { return this->_stalls_left > 0; }

double ExtrusionUnit::get_prob_of_extr_unit_bypass() const {
  return this->_parent_lef.get_probability_of_extr_unit_bypass();
}

bool ExtrusionUnit::try_extrude() {
  assert(this->_direction == DNA::Direction::fwd || this->_direction == DNA::Direction::rev);
  if (this->is_stalled()) {
    this->decrement_stalls();
    return false;
  }

  if (this->_direction == DNA::Direction::fwd)
    return this->try_moving_to_next_bin();
  else
    return this->try_moving_to_prev_bin();
}

bool ExtrusionUnit::try_moving_to_next_bin() {
  //  absl::FPrintF(stderr, "Trying to move to next bin...");
  if (this->_bin == &this->_parent_lef.get_last_bin()) {
    this->set_stalls(UINT32_MAX);  // Stall forever
    return false;
  }
  //  absl::FPrintF(stderr, " Moving from bin #%lu", this->_bin->get_index());
  this->_bin->remove_extr_unit_binding(this);
  this->_bin = &this->_parent_lef.get_ptr_to_chr()->dna.get_next_bin(*this->_bin);
  this->_bin->add_extr_unit_binding(this);
  //  absl::FPrintF(stderr, " to bin #%lu!\n", this->_bin->get_index());
  return true;
}

bool ExtrusionUnit::try_moving_to_prev_bin() {
  if (this->_bin == &this->_parent_lef.get_first_bin()) {
    this->set_stalls(UINT32_MAX);  // Stall forever
    return false;
  }
  //  absl::FPrintF(stderr, " Moving from bin #%lu", this->_bin->get_index());
  this->_bin->remove_extr_unit_binding(this);
  this->_bin = &this->_parent_lef.get_ptr_to_chr()->dna.get_prev_bin(*this->_bin);
  this->_bin->add_extr_unit_binding(this);
  //  absl::FPrintF(stderr, " to bin #%lu!\n", this->_bin->get_index());
  return true;
}

bool ExtrusionUnit::hard_stall() const {
  return this->_blocking_barrier &&
         this->get_extr_direction() == this->_blocking_barrier->get_direction();
}

void ExtrusionUnit::set_stalls(uint32_t n) { this->_stalls_left = n; }
void ExtrusionUnit::increment_stalls(uint32_t n) {
  assert(UINT32_MAX - n >= this->_stalls_left);
  this->_stalls_left += n;
}

void ExtrusionUnit::decrement_stalls(uint32_t n) {
  assert(this->_stalls_left >= n);
  this->_stalls_left -= n;
}

void ExtrusionUnit::reset_stalls() { this->set_stalls(0); }

void ExtrusionUnit::unload() {
  this->reset_stalls();
  this->_bin->remove_extr_unit_binding(this);
  this->_bin = nullptr;
  this->_blocking_barrier = nullptr;
}

void ExtrusionUnit::bind(Chromosome* chr, uint32_t pos, DNA::Direction direction) {
  assert(pos < chr->length());
  assert(direction == DNA::Direction::fwd || direction == DNA::Direction::rev);
  this->_bin = chr->dna.get_ptr_to_bin_from_pos(pos);
  this->_bin->add_extr_unit_binding(this);
  this->_direction = direction;
}

bool ExtrusionUnit::is_bound() const { return this->_bin != nullptr; }

uint64_t ExtrusionUnit::check_constraints(std::mt19937& rand_gen) {
  if (this->is_stalled()) return 0;
  // This already applies the appropriate stall to colliding units.
  // Communicating this to Lef should not be necessary.
  if (const auto stalls = this->check_for_extruder_collisions(rand_gen); stalls > 0) return stalls;
  // This already figures out if whether we should apply a "big" or "small" stall,
  // based on the direction of extrusion as well as the direction of the barrier
  // that is blocking the extrusion unit
  return this->check_for_extrusion_barrier(rand_gen);
}

uint32_t ExtrusionUnit::check_for_extruder_collisions(std::mt19937& rand_gen) {
  assert(!this->is_stalled());
  assert(this->_direction == DNA::Direction::fwd || this->_direction == DNA::Direction::rev);
  assert(this->_bin->get_n_extr_units() > 0);
  uint32_t n_stalls = 0;
  if (this->_bin->get_n_extr_units() == 1) return n_stalls;

  for (auto& other : this->_bin->get_extr_units()) {
    assert(other->_direction == DNA::Direction::fwd || other->_direction == DNA::Direction::rev);
    if (&this->_parent_lef != &other->_parent_lef && this->_direction != other->_direction) {
      std::geometric_distribution<uint32_t> dist(this->get_prob_of_extr_unit_bypass());
      n_stalls = dist(rand_gen);
      this->increment_stalls(n_stalls);
      if (!other->is_stalled()) other->increment_stalls(n_stalls);
      break;
    }
  }
  return n_stalls;
}

uint64_t ExtrusionUnit::check_for_extrusion_barrier(std::mt19937& rand_gen) {
  assert(!this->is_stalled());
  uint64_t applied_stall = 0;
  this->_blocking_barrier = this->_bin->get_next_extr_barrier(this->_blocking_barrier);
  if (this->_blocking_barrier) {
    applied_stall = this->_blocking_barrier->generate_num_of_blocking_events(rand_gen);
    if (this->_blocking_barrier->get_direction() != this->get_extr_direction()) {
      // "Small" stall
      applied_stall /= 2;  // TODO: Make this tunable
    }
    this->increment_stalls(applied_stall);
    //    if (applied_stall > 0) {
    //      absl::FPrintF(
    //          stderr, "Applying a %s stall of %lu. # of stalls left: %lu\n",
    //          this->_blocking_barrier->get_direction() == this->get_extr_direction() ? "big" :
    //          "small", applied_stall, this->_stalls_left + applied_stall);
    //    }
  }
  return applied_stall;
}

Lef::Lef(uint32_t bin_size, uint32_t avg_processivity, double probability_of_extruder_bypass)
    : _avg_processivity(avg_processivity),
      _probability_of_extr_unit_bypass(probability_of_extruder_bypass),
      _lifetime_generator(compute_prob_of_unloading(bin_size)) {}

bool Lef::is_bound() const {
  assert((this->_left_unit->is_bound() && this->_right_unit->is_bound()) ||
         (!this->_left_unit->is_bound() && !this->_right_unit->is_bound()));
  return this->_left_unit->is_bound();
}

void Lef::bind_at_pos(Chromosome& chr, uint32_t pos, std::mt19937& rand_eng) {
  this->_chr = &chr;
  this->_lifetime = this->_lifetime_generator(rand_eng);
  // We assume that the left unit always travels towards the 5', while the right unit goes to the 3'
  if (!this->_left_unit)
    this->_left_unit =
        std::make_unique<ExtrusionUnit>(*this, this->_probability_of_extr_unit_bypass);
  if (!this->_right_unit)
    this->_right_unit =
        std::make_unique<ExtrusionUnit>(*this, this->_probability_of_extr_unit_bypass);
  this->_left_unit->bind(this->_chr, pos, DNA::Direction::rev);
  this->_right_unit->bind(this->_chr, pos, DNA::Direction::fwd);
}

std::string_view Lef::get_chr_name() const { return this->_chr->name; }

std::pair<DNA::Bin*, DNA::Bin*> Lef::get_ptr_to_bins() {
  return std::make_pair(this->_left_unit->_bin, this->_left_unit->_bin);
}

void Lef::unload() {
  //  absl::FPrintF(stderr, "Unloading LEF at bin %lu-%lu pos %lu-%lu...\n",
  //                this->_left_unit->_bin->get_index(), this->_right_unit->_bin->get_index(),
  //                this->_left_unit->get_pos(), this->_right_unit->get_pos());
  this->_chr = nullptr;
  //  (void)this->_left_unit->_bin->remove_extr_unit_binding(&this->_left_unit);
  //  (void)this->_right_unit->_bin->remove_extr_unit_binding(&this->_right_unit);
  this->_left_unit->unload();
  this->_right_unit->unload();
  this->_lifetime = 0;  // Probably unnecessary
}

uint32_t Lef::extrude() {
  //  if (this->get_loop_size() > this->_avg_processivity * 2) {
  //    absl::FPrintF(stderr, "The loop generated by lef %p looks suspiciously large: %lu\n", this,
  //                  this->get_loop_size());
  //  }

  if (this->_lifetime-- > 0) {
    return this->_left_unit->try_extrude() + this->_right_unit->try_extrude();
  }
  this->unload();
  return 0;
}

void Lef::register_contact() {
  this->_chr->contacts.increment(this->_left_unit->_bin->get_index(),
                                 this->_right_unit->_bin->get_index());
}

void Lef::check_constraints(std::mt19937& rand_gen) {
  if (this->_left_unit->is_stalled() && this->_right_unit->is_stalled()) {
    //    absl::FPrintF(stderr, "Both extrusion units of lef %p are stalled!\n", this);
    return;
  }
  uint32_t n_of_stalls_due_to_barriers = std::max(this->_left_unit->check_constraints(rand_gen),
                                                  this->_right_unit->check_constraints(rand_gen));
  //  if (n_of_stalls_due_to_barriers > 0) {
  //    absl::FPrintF(stderr, "n_of_stalls_due_to_barriers=%lu\n", n_of_stalls_due_to_barriers);
  //  }
  if (this->_left_unit->is_stalled() && this->_right_unit->is_stalled()) {
    n_of_stalls_due_to_barriers *= 2;
    if (this->_left_unit->hard_stall() && this->_right_unit->hard_stall()) {
      n_of_stalls_due_to_barriers *= 2;
      assert(UINT32_MAX - n_of_stalls_due_to_barriers > this->_lifetime);
      this->_lifetime += n_of_stalls_due_to_barriers;
    }

    //    absl::FPrintF(stderr, "left_stall %lu -> %lu\n", this->_left_unit->_stalls_left,
    //                  std::max(n_of_stalls_due_to_barriers, this->_left_unit->_stalls_left));
    //    absl::FPrintF(stderr, "right_stall %lu -> %lu\n", this->_right_unit->_stalls_left,
    //                  std::max(n_of_stalls_due_to_barriers, this->_right_unit->_stalls_left));

    this->_left_unit->_stalls_left =
        std::max(n_of_stalls_due_to_barriers, this->_left_unit->_stalls_left);
    this->_right_unit->_stalls_left =
        std::max(n_of_stalls_due_to_barriers, this->_right_unit->_stalls_left);
    //    usleep(100000);
  }
}

uint32_t Lef::get_loop_size() const {
  assert(this->_right_unit->get_pos() >= this->_left_unit->get_pos());
  return this->_right_unit->get_pos() - this->_left_unit->get_pos();
}

uint32_t Lef::get_avg_processivity() const { return this->_avg_processivity; }

bool Lef::try_rebind(Chromosome& chr, std::mt19937& rand_eng, double prob_of_rebinding) {
  assert(!this->is_bound());
  assert(prob_of_rebinding >= 0 && prob_of_rebinding <= 1);
  std::uniform_real_distribution<> d1(0.0, 1.0);
  if (d1(rand_eng) >= 1 - prob_of_rebinding) {
    //    absl::FPrintF(stderr, "Trying to rebind...\n");
    std::uniform_int_distribution<uint32_t> d2(0, chr.length() - 1);
    this->bind_at_pos(chr, d2(rand_eng), rand_eng);
    return true;
  }
  return false;
}

const DNA::Bin& Lef::get_first_bin() const { return this->_chr->dna.get_first_bin(); }
const DNA::Bin& Lef::get_last_bin() const { return this->_chr->dna.get_last_bin(); }

Chromosome* Lef::get_ptr_to_chr() {
  assert(this->_chr);
  return this->_chr;
}

std::pair<uint32_t, uint32_t> Lef::get_pos() const {
  return std::make_pair(this->_left_unit->get_pos(), this->_right_unit->get_pos());
}

double Lef::get_probability_of_extr_unit_bypass() const {
  return this->_probability_of_extr_unit_bypass;
}

double Lef::compute_prob_of_unloading(uint32_t bin_size, uint8_t n_of_active_extr_units) const {
  /*
  const double avg_number_of_extrusion_events =
      this->_avg_processivity / static_cast<double>(n_of_active_extr_units);
  absl::FPrintF(stderr, "Prob of unloading = %.6f\n", 1.0 / avg_number_of_extrusion_events);
  return 1.0 / avg_number_of_extrusion_events;
   */
  return n_of_active_extr_units / (static_cast<double>(this->get_avg_processivity()) / bin_size);
}

}  // namespace modle