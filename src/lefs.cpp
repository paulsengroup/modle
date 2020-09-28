#include "modle/lefs.hpp"

#include <absl/strings/str_format.h>

namespace modle {
Lef::Lef(uint32_t binding_pos, uint32_t avg_processivity)
    : _left_pos(binding_pos),
      _left_bin(),
      _right_pos(binding_pos),
      _avg_processivity(avg_processivity),
      _dist(this->_avg_processivity) {}

Lef::Lef(uint32_t binding_pos, uint32_t avg_processivity, uint8_t left_extr_speed,
         uint8_t right_extr_speed)
    : _left_pos(binding_pos),
      _right_pos(binding_pos),
      _left_bin(),
      _left_extrusion_speed(left_extr_speed),
      _right_extrusion_speed(right_extr_speed),
      _avg_processivity(avg_processivity),
      _dist(this->_avg_processivity) {}

bool Lef::is_bound() const { return this->_left_bin && this->_right_bin; }

bool Lef::bind_at_pos(std::string_view chr_name, DNA &chr, uint32_t pos) {
  auto bin = chr.get_ptr_to_bin_from_pos(pos);
  // 1st ref. is for std::vector, 2nd ref. is for the instance passed to this func (and any
  // additional ref. is due to a lef binding to this bin)
  if (bin.use_count() - 1 > 1) return false;
  this->_chr = chr_name;
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

Lsf::Lsf(uint32_t left_pos, uint32_t right_pos, uint32_t lifetime)
    : _left_pos(left_pos), _right_pos(right_pos), _lifetime(lifetime) {}

bool Lsf::next() {
  this->_lifetime -= this->_lifetime > 0;
  return this->_lifetime;
}

}  // namespace modle