#pragma once

#include "modle/common.hpp"

namespace modle {

constexpr bp_t Chromosome::start_pos() const { return this->_start; }

constexpr bp_t Chromosome::end_pos() const { return this->_end; }

constexpr bp_t Chromosome::size() const { return this->_size; }
constexpr bp_t Chromosome::simulated_size() const { return this->_end - this->_start; }
}  // namespace modle
