#include "modle/lefs.hpp"

#include "absl/strings/str_format.h"
#include "modle/contacts.hpp"

namespace modle {

ExtrusionUnit::ExtrusionUnit(Lef* lef) : _lef_ptr(lef) {}

ExtrusionUnit::ExtrusionUnit(Lef* lef, uint32_t pos, DNA* dna, std::shared_ptr<DNA::Bin> bin,
                             uint32_t extrusion_speed, ExtrusionUnit::Direction direction,
                             bool stall)
    : _lef_ptr(lef),
      _pos(pos),
      _dna(dna),
      _bin(std::move(bin)),
      _extrusion_speed(extrusion_speed),
      _direction(direction),
      _stalled(stall) {}

uint32_t ExtrusionUnit::extrude() {
  if (this->_direction == ExtrusionUnit::Direction::fwd) {
    if (this->_pos < (this->_dna->length() - this->get_extrusion_speed()) &&
        (this->_pos += this->get_extrusion_speed()) > this->_bin->get_end()) {
      //      absl::FPrintF(stderr, "Switching to next bin at %lu!\n", this->_pos);
      (void)this->_bin->remove_extr_unit_binding(this);
      this->_bin = this->_dna->get_ptr_to_next_bin(this->_bin.get());
      (void)this->_bin->add_extr_unit_binding(this);
    }
  } else {
    if (this->_pos > this->get_extrusion_speed() &&
        (this->_pos -= this->get_extrusion_speed()) < this->_bin->get_start()) {
      //      absl::FPrintF(stderr, "Switching to previous bin at %lu!\n", this->_pos);
      (void)this->_bin->remove_extr_unit_binding(this);
      this->_bin = this->_dna->get_ptr_to_previous_bin(this->_bin.get());
      (void)this->_bin->add_extr_unit_binding(this);
    }
  }
  return this->get_extrusion_speed();
}

uint32_t ExtrusionUnit::get_pos() const { return this->_pos; }

void ExtrusionUnit::set_stall() {
  //  absl::FPrintF(stderr, "STALLING LEF at pos %lu dir %s...\n", this->get_pos(),
  //                this->_direction ? "FWD" : "REV");
  this->_stalled = true;
}
void ExtrusionUnit::remove_stall() { this->_stalled = false; }
void ExtrusionUnit::unload() {
  this->_dna = nullptr;
  this->_pos = -1;
  (void)this->_bin->remove_extr_unit_binding(this);
  this->_bin = nullptr;
  this->remove_stall();
}

bool ExtrusionUnit::bind(DNA* dna, uint32_t pos, std::shared_ptr<DNA::Bin> bin,
                         Direction direction) {
  this->_dna = dna;
  this->_pos = pos;
  this->_direction = direction;
  this->_bin = std::move(bin);
  (void)this->_bin->add_extr_unit_binding(this);

  return true;
}

bool ExtrusionUnit::is_extruding_fwd() const {
  return this->_direction == ExtrusionUnit::Direction::fwd;
}
bool ExtrusionUnit::is_extruding_rev() const { return !this->is_extruding_fwd(); }

bool ExtrusionUnit::is_bound() const { return this->_bin != nullptr; }

uint32_t ExtrusionUnit::get_extrusion_speed() const {
  return !this->_stalled * this->_extrusion_speed;
}
void ExtrusionUnit::check_for_extruder_collisions() {
  assert(this->_bin->n_extruders() > 1);
  for (auto& other : this->_bin->get_extruders()) {
    if (this->_lef_ptr != other->_lef_ptr && this->_direction != other->_direction &&
        std::abs(static_cast<int64_t>(this->_pos) - other->_pos) <=
            (this->get_extrusion_speed() + other->get_extrusion_speed())) {
      this->set_stall();
      other->set_stall();
      return;
    }
  }
}

Lef::Lef(uint32_t avg_processivity)
    : _avg_processivity(avg_processivity), _bernoulli_dist(compute_prob_of_unloading()) {}

uint32_t Lef::get_left_extrusion_speed() const { return this->_left_unit.get_extrusion_speed(); }

uint32_t Lef::get_right_extrusion_speed() const { return this->_right_unit.get_extrusion_speed(); }

bool Lef::left_is_stalled() const { return this->_left_unit._stalled; }
bool Lef::right_is_stalled() const { return this->_right_unit._stalled; }
bool Lef::is_bound() const {
  assert((this->_left_unit.is_bound() && this->_right_unit.is_bound()) ||
         (!this->_left_unit.is_bound() && !this->_right_unit.is_bound()));
  return this->_left_unit.is_bound();
}

bool Lef::bind_at_pos(std::string_view chr_name, DNA& dna, modle::ContactMatrix& contacts,
                      uint32_t pos) {
  auto bin = dna.get_ptr_to_bin_from_pos(pos);
  // 1st ref. is for std::vector, 2nd ref. is for the instance passed to this func (and any
  // additional ref. is due to a lef binding to this bin)
  //  if (bin.use_count() - 1 > 1) return false;
  this->_chr = chr_name;
  this->_dna = &dna;
  this->_contacts = &contacts;
  // We assume that the left unit always travels towards the 5', while the right unit goes to the 3'
  this->_left_unit.bind(&dna, pos, bin, ExtrusionUnit::Direction::rev);
  this->_right_unit.bind(&dna, pos, std::move(bin), ExtrusionUnit::Direction::fwd);
  //  absl::FPrintF(stderr, "LEF bound at bin %lu pos %lu...\n", this->_left_unit._bin->get_index(),
  //  pos);

  return true;
}

std::string_view Lef::get_chr_name() const { return this->_chr; }

std::pair<DNA::Bin*, DNA::Bin*> Lef::get_ptr_to_bins() {
  return std::make_pair(this->_left_unit._bin.get(), this->_left_unit._bin.get());
}

void Lef::unload() {
  //  absl::FPrintF(stderr, "Unloading LEF at bin %lu-%lu pos %lu-%lu...\n",
  //                this->_left_unit._bin->get_index(), this->_right_unit._bin->get_index(),
  //                this->_left_unit.get_pos(), this->_right_unit.get_pos());
  this->_chr = "";
  this->_dna = nullptr;
  this->_contacts = nullptr;
  //  (void)this->_left_unit._bin->remove_extr_unit_binding(&this->_left_unit);
  //  (void)this->_right_unit._bin->remove_extr_unit_binding(&this->_right_unit);
  this->_left_unit.unload();
  this->_right_unit.unload();
  this->_bernoulli_dist.reset();  // TODO: check if this makes sense:
  // https://en.cppreference.com/w/cpp/numeric/random/bernoulli_distribution/reset
}

 uint32_t Lef::extrude() { return this->_left_unit.extrude() + this->_right_unit.extrude(); }
 /*
uint32_t Lef::extrude() {
  auto l = this->_left_unit.extrude();
  auto r = this->_right_unit.extrude();
  if (l + r == 0) absl::FPrintF(stderr, "LEF stalled left and right.\n");
  else if (l == 0) absl::FPrintF(stderr, "LEF stalled left.\n");
  else if (r == 0) absl::FPrintF(stderr, "LEF stalled right.\n");
  return l + r;
}
  */

void Lef::register_contact() {
  this->_contacts->increment(this->_left_unit._bin->get_index(),
                             this->_right_unit._bin->get_index());
}

void Lef::check_constraints() {
  // Check for clashes with other LEFs
  // Check for boundaries
  if (this->_left_unit._pos <=
          (this->_left_unit._bin->get_start() + this->_left_unit.get_extrusion_speed()) &&
      this->_left_unit._bin->blocking_rev()) {
    //    absl::FPrintF(stderr, "Stalling left because of an extr barrier....\n");
    //    sleep(1);
    this->_left_unit.set_stall();
    //  } else {
    //    this->remove_left_stall();
  }
  if (this->_right_unit._pos >=
          (this->_right_unit._bin->get_end() - this->_right_unit.get_extrusion_speed()) &&
      this->_right_unit._bin->blocking_fwd()) {
    //    absl::FPrintF(stderr, "Stalling right because of an extr barrier....\n");
    //    sleep(1);
    this->_right_unit.set_stall();
    //  } else {
    //    this->remove_right_stall();
  }
  if (const auto n_extr = this->_left_unit._bin->n_extruders(); n_extr > 1) {
    this->_left_unit.check_for_extruder_collisions();
  } else if (n_extr == 1 && !this->_left_unit._bin->blocking_rev()) {
    this->_left_unit.remove_stall();
  }
  if (const auto n_extr = this->_right_unit._bin->n_extruders(); n_extr > 1) {
    this->_right_unit.check_for_extruder_collisions();
  } else if (n_extr == 1 && !this->_right_unit._bin->blocking_fwd()) {
    this->_right_unit.remove_stall();
  }
}

uint32_t Lef::get_loop_size() const {
  assert(this->_right_unit._pos >= this->_left_unit._pos);
  return this->_right_unit._pos - this->_left_unit._pos;
}

uint32_t Lef::get_avg_processivity() const { return this->_avg_processivity; }

void Lef::set_avg_processivity(uint32_t avg_proc) {
  assert(avg_proc % 2 == 0);
  this->_avg_processivity = avg_proc;
  this->_left_unit._extrusion_speed = this->_avg_processivity / 2;
  this->_right_unit._extrusion_speed = this->_avg_processivity / 2;
}

bool Lef::try_unload(std::default_random_engine& rng) {
  if (this->_bernoulli_dist(rng)) {
    //    absl::FPrintF(stderr, "Unloading!\n");
    this->unload();
    return true;
  }
  return false;
}

bool Lef::try_rebind(std::string_view chr_name, DNA& chr, modle::ContactMatrix& contacts,
                     std::default_random_engine& rng, double prob_of_rebinding) {
  assert(!this->is_bound());
  assert(prob_of_rebinding >= 0 && prob_of_rebinding <= 1);
  std::uniform_int_distribution<uint32_t> d1(0, chr.length());
  std::uniform_real_distribution<> d2(0.0, 1.0);
  if (d2(rng) <= prob_of_rebinding) {
    //    absl::FPrintF(stderr, "Trying to rebind...\n");
    return this->bind_at_pos(chr_name, chr, contacts, d1(rng));
  }
  return false;
}

uint32_t Lef::get_left_pos() const { return this->_left_unit.get_pos(); }
uint32_t Lef::get_right_pos() const { return this->_right_unit.get_pos(); }

double Lef::compute_prob_of_unloading() const {
  const double bp_per_extr_event =
      this->get_left_extrusion_speed() + this->get_right_extrusion_speed();
  const double mean_extr_events = this->get_avg_processivity() / bp_per_extr_event;
  absl::FPrintF(
      stderr,
      "avg_processivity=%lu; bp_extruded=%.0f; mean_extr_events=%.0f; prob_unloading=%.10f\n",
      this->get_avg_processivity(), bp_per_extr_event, mean_extr_events, 1.0 / mean_extr_events);
  return 1.0 / mean_extr_events;
}

}  // namespace modle