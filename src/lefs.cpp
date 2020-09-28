#include "modle/lefs.hpp"

#include <absl/strings/str_format.h>

namespace modle {
Lef::Lef(uint32_t binding_pos, uint32_t avg_processivity, DNA* dna, std::string_view chr)
    : _chr(chr),
      _dna(dna),
      _left_pos(binding_pos),
      _left_bin(),
      _right_pos(binding_pos),
      _avg_processivity(avg_processivity),
      _dist(this->_avg_processivity) {}

Lef::Lef(uint32_t binding_pos, uint32_t avg_processivity, uint8_t left_extr_speed,
         uint8_t right_extr_speed, DNA* dna, std::string_view chr)
    : _chr(chr),
      _dna(dna),
      _left_pos(binding_pos),
      _right_pos(binding_pos),
      _left_bin(),
      _left_extrusion_speed(left_extr_speed),
      _right_extrusion_speed(right_extr_speed),
      _avg_processivity(avg_processivity),
      _dist(this->_avg_processivity) {}

uint32_t Lef::get_left_extrusion_speed() const {
  return !this->left_is_stalled() * this->_left_extrusion_speed;
}

uint32_t Lef::get_right_extrusion_speed() const {
  return !this->right_is_stalled() * this->_right_extrusion_speed;
}

bool Lef::left_is_stalled() const { return this->_stall_left; }
bool Lef::right_is_stalled() const { return this->_stall_right; }

bool Lef::is_bound() const { return this->_left_bin && this->_right_bin; }

bool Lef::bind_at_pos(std::string_view chr_name, DNA& dna, uint32_t pos) {
  auto bin = dna.get_ptr_to_bin_from_pos(pos);
  // 1st ref. is for std::vector, 2nd ref. is for the instance passed to this func (and any
  // additional ref. is due to a lef binding to this bin)
  if (bin.use_count() - 1 > 1) return false;
  this->_chr = chr_name;
  this->_dna = &dna;
  this->_left_pos = pos;
  this->_right_pos = pos;
  this->_left_bin = bin;
  this->_right_bin = std::move(bin);

  return true;
}

std::string_view Lef::get_chr_name() const { return this->_chr; }

std::pair<std::shared_ptr<DNA::Bin>, std::shared_ptr<DNA::Bin>> Lef::get_bins() {
  return std::make_pair(this->_left_bin, this->_right_bin);
}

void Lef::unload() {
  this->_chr = "";
  this->_dna = nullptr;
  this->_left_pos = -1;
  this->_right_pos = -1;
  this->_left_bin = nullptr;
  this->_right_bin = nullptr;
  this->remove_left_stall();
  this->remove_right_stall();
  this->_dist
      .reset();  // TODO: check if this makes sense:
                 // https://en.cppreference.com/w/cpp/numeric/random/bernoulli_distribution/reset
}

void Lef::stall_left() { this->_stall_left = true; }
void Lef::stall_right() { this->_stall_right = true; }

void Lef::remove_left_stall() { this->_stall_left = false; }
void Lef::remove_right_stall() { this->_stall_right = false; }

bool Lef::extrude(std::default_random_engine& rng) {
  if (this->_dist(rng)) {
    if (this->_left_pos > 0 &&
        (this->_left_pos -= this->get_left_extrusion_speed()) < this->_left_bin->get_start()) {
      this->_left_bin = this->_dna->get_ptr_to_previous_bin(this->_left_bin);
    }
    if (this->_right_pos < this->_dna->length() &&
        (this->_right_pos += this->get_right_extrusion_speed()) > this->_right_bin->get_end()) {
      this->_right_bin = this->_dna->get_ptr_to_next_bin(this->_right_bin);
    }
    return true;
  }
  this->unload();
  return false;
}

void Lef::check_constraints() {
  // Check for clashes with other LEFs
  // Check for boundaries
  if (this->_left_pos == this->_left_bin->get_start() && this->_left_bin->has_rev_barrier()) {
    this->stall_left();
  }
  if (this->_right_pos == this->_right_bin->get_end() && this->_right_bin->has_fwd_barrier()) {
    this->stall_right();
  }
  // Instead of always unloading, I should incorporate somewhere the probability of unloading (based on process. I guess)
  if (this->left_is_stalled() && this->right_is_stalled()) this->unload();
}

uint32_t Lef::get_loop_size() const { return this->_right_pos - this->_left_pos; }

Lsf::Lsf(uint32_t left_pos, uint32_t right_pos, uint32_t lifetime)
    : _left_pos(left_pos), _right_pos(right_pos), _lifetime(lifetime) {}

bool Lsf::next() {
  this->_lifetime -= this->_lifetime > 0;
  return this->_lifetime;
}

}  // namespace modle