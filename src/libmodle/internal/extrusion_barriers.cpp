// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "modle/extrusion_barriers.hpp"

#include <cpp-sort/sorters/pdq_sorter.h>  // for pdq_sort

#include <algorithm>  // for clamp
#include <cassert>    // for assert
#include <limits>     // for numeric_limits
#include <vector>     // for vector

#include "modle/common/common.hpp"  // for bp_t, Direction, fwd, none
#include "modle/common/dna.hpp"     // for dna::Direction, dna::NONE
#include "modle/common/random.hpp"  // for PRNG_t

namespace modle {

ExtrusionBarriers::ExtrusionBarriers(usize size)
    : _pos(size),
      _direction(size, dna::NONE),
      _stp_active(size),
      _stp_inactive(size),
      _state(size) {}

ExtrusionBarriers::ExtrusionBarriers(std::initializer_list<bp_t> pos,
                                     std::initializer_list<dna::Direction> direction,
                                     std::initializer_list<TP> stp_active,    // NOLINT
                                     std::initializer_list<TP> stp_inactive,  // NOLINT
                                     std::initializer_list<State> state, bool sort)
    : _pos(pos),
      _direction(direction),
      _stp_active(stp_active),
      _stp_inactive(stp_inactive),
      _state(state) {
  this->assert_buffer_sizes_are_equal();
  assert(std::all_of(this->_direction.begin(), this->_direction.end(),
                     [](const auto d) { return d != dna::NONE; }));

  if (sort) {
    this->sort();
  }
}

ExtrusionBarriers::ExtrusionBarriers(std::initializer_list<ExtrusionBarrier> barriers,
                                     std::initializer_list<State> states, bool sort)
    : ExtrusionBarriers(barriers.size()) {
  assert(this->size() == states.size());
  this->_state = states;

  for (usize i = 0; i < barriers.size(); ++i) {
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    const auto& b = *(std::begin(barriers) + static_cast<isize>(i));

    this->_pos[i] = b.pos;
    this->_direction[i] = b.blocking_direction;
    this->_stp_active[i] = b.stp_active;
    this->_stp_inactive[i] = b.stp_inactive;
  }

  if (sort) {
    this->sort();
  }
}

usize ExtrusionBarriers::size() const noexcept {
  this->assert_buffer_sizes_are_equal();
  return this->_pos.size();
}

bool ExtrusionBarriers::empty() const noexcept { return this->size() == 0; }

void ExtrusionBarriers::resize(usize new_size) {
  this->_pos.resize(new_size);
  this->_direction.resize(new_size, dna::NONE);
  this->_stp_active.resize(new_size);
  this->_stp_inactive.resize(new_size);
  this->_state.resize(new_size);

  this->assert_buffer_sizes_are_equal();
}

void ExtrusionBarriers::clear() noexcept {
  this->_pos.clear();
  this->_direction.clear();
  this->_stp_active.clear();
  this->_stp_inactive.clear();
  this->_state.clear();

  this->assert_buffer_sizes_are_equal();
}

bp_t ExtrusionBarriers::pos(usize i) const noexcept {
  this->assert_index_within_bounds(i);
  return this->_pos[i];
}
const std::vector<bp_t>& ExtrusionBarriers::pos() const noexcept { return this->_pos; }

dna::Direction ExtrusionBarriers::direction(usize i) const noexcept {
  this->assert_index_within_bounds(i);
  return this->_direction[i];
}
const std::vector<dna::Direction>& ExtrusionBarriers::direction() const noexcept {
  return this->_direction;
}

double ExtrusionBarriers::stp_active(usize i) const noexcept {
  this->assert_index_within_bounds(i);
  return this->_stp_active[i]();
}
auto ExtrusionBarriers::stp_active() const noexcept -> const std::vector<TP>& {
  return this->_stp_active;
}

double ExtrusionBarriers::stp_inactive(usize i) const noexcept {
  this->assert_index_within_bounds(i);
  return this->_stp_inactive[i]();
}
auto ExtrusionBarriers::stp_inactive() const noexcept -> const std::vector<TP>& {
  return this->_stp_inactive;
}

auto ExtrusionBarriers::state(usize i) const noexcept -> State {
  this->assert_index_within_bounds(i);
  return this->_state[i];
}
auto ExtrusionBarriers::state() const noexcept -> const std::vector<State>& { return this->_state; }

bool ExtrusionBarriers::is_active(usize i) const noexcept {
  this->assert_index_within_bounds(i);
  return this->_state[i] == State::ACTIVE;
}

bool ExtrusionBarriers::is_not_active(usize i) const noexcept {
  this->assert_index_within_bounds(i);
  return this->_state[i] == State::INACTIVE;
}

double ExtrusionBarriers::occupancy(usize i) const noexcept {
  this->assert_index_within_bounds(i);
  return ExtrusionBarrier::compute_occupancy_from_stp(this->stp_active(i), this->stp_inactive(i));
}

auto ExtrusionBarriers::next_state(usize i, random::PRNG_t& rand_eng) noexcept -> State {
  const auto p = random::generate_canonical<double, std::numeric_limits<double>::digits>(rand_eng);

  if (this->is_not_active(i) && p > this->stp_inactive(i)) {
    this->_state[i] = State::ACTIVE;
  } else if (this->is_active(i) && p > this->stp_active(i)) {
    this->_state[i] = State::INACTIVE;
  }

  return this->_state[i];
}

void ExtrusionBarriers::next_state(random::PRNG_t& rand_eng) noexcept {
  for (usize i = 0; i < this->size(); ++i) {
    this->next_state(i, rand_eng);
  }
}

usize ExtrusionBarriers::count_active() const noexcept {
  return static_cast<usize>(std::count(this->_state.begin(), this->_state.end(), State::ACTIVE));
}
usize ExtrusionBarriers::count_inactive() const noexcept {
  const auto num_active =
      static_cast<usize>(std::count(this->_state.begin(), this->_state.end(), State::ACTIVE));
  assert(num_active <= this->size());
  return this->size() - num_active;
}

void ExtrusionBarriers::assert_buffer_sizes_are_equal() const noexcept {
  assert(this->_pos.size() == this->_direction.size());
  assert(this->_pos.size() == this->_stp_active.size());
  assert(this->_pos.size() == this->_stp_inactive.size());
  assert(this->_pos.size() == this->_state.size());
}

void ExtrusionBarriers::assert_index_within_bounds([[maybe_unused]] usize i) const noexcept {
  assert(i < this->size());
}

void ExtrusionBarriers::set(usize i, bp_t pos, dna::Direction direction, TP stp_active,
                            TP stp_inactive, State state) noexcept {
  this->assert_index_within_bounds(i);
  assert(direction != dna::NONE);
  this->_pos[i] = pos;
  this->_direction[i] = direction;
  this->_stp_active[i] = stp_active;
  this->_stp_inactive[i] = stp_inactive;
  this->_state[i] = state;
}

void ExtrusionBarriers::push_back(bp_t pos, dna::Direction direction, TP stp_active,
                                  TP stp_inactive, State state) {
  assert(direction != dna::NONE);
  this->_pos.push_back(pos);
  this->_direction.push_back(direction);
  this->_stp_active.push_back(stp_active);
  this->_stp_inactive.push_back(stp_inactive);
  this->_state.push_back(state);
}

void ExtrusionBarriers::set(usize i, const ExtrusionBarrier& barrier, State state) noexcept {
  this->set(i, barrier.pos, barrier.blocking_direction, barrier.stp_active, barrier.stp_inactive,
            state);
}

void ExtrusionBarriers::push_back(const ExtrusionBarrier& barrier, State state) {
  this->push_back(barrier.pos, barrier.blocking_direction, barrier.stp_active, barrier.stp_inactive,
                  state);
}

void ExtrusionBarriers::set(usize i, State state) noexcept {
  this->assert_index_within_bounds(i);
  this->_state[i] = state;
}

auto ExtrusionBarriers::init_state(usize i, random::PRNG_t& rand_eng) noexcept -> State {
  this->assert_index_within_bounds(i);
  this->_state[i] =
      random::bernoulli_trial{this->occupancy(i)}(rand_eng) ? State::ACTIVE : State::INACTIVE;

  return this->_state[i];
}

void ExtrusionBarriers::init_states(random::PRNG_t& rand_eng) noexcept {
  for (usize i = 0; i < this->size(); ++i) {
    this->init_state(i, rand_eng);
  }
}

void ExtrusionBarriers::sort() {
  std::vector<usize> buff(this->size());
  this->sort(buff);
}

void ExtrusionBarriers::sort(std::vector<usize>& idx_buff) {
  idx_buff.resize(this->size());
  std::iota(idx_buff.begin(), idx_buff.end(), usize(0));

  cppsort::pdq_sort(idx_buff.begin(), idx_buff.end(),
                    [&](const auto i1, const auto i2) { return this->pos(i1) < this->pos(i2); });
  for (usize i = 0; i < this->size(); ++i) {
    // https://stackoverflow.com/a/22218699
    while (idx_buff[i] != i) {
      const auto j = idx_buff[i];
      const auto k = idx_buff[j];

      std::swap(this->_pos[j], this->_pos[k]);
      std::swap(this->_direction[j], this->_direction[k]);
      std::swap(this->_stp_active[j], this->_stp_active[k]);
      std::swap(this->_stp_inactive[j], this->_stp_inactive[k]);
      std::swap(this->_state[j], this->_state[k]);
      std::swap(idx_buff[i], idx_buff[j]);
    }
  }
}

}  // namespace modle
