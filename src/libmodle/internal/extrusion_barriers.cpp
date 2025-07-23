// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "modle/extrusion_barriers.hpp"

#include <cpp-sort/sorters/pdq_sorter.h>

#include <algorithm>
#include <cassert>
#include <limits>
#include <vector>

#include "modle/common/common.hpp"
#include "modle/common/dna.hpp"
#include "modle/common/random.hpp"

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
                                     std::initializer_list<State> state, bool sort_barriers)
    : _pos(pos),
      _direction(direction),
      _stp_active(stp_active),
      _stp_inactive(stp_inactive),
      _state(state) {
  assert_buffer_sizes_are_equal();
  assert(std::all_of(_direction.begin(), _direction.end(),
                     [](const auto d) { return d != dna::NONE; }));

  if (sort_barriers) {
    sort();
  }
}

ExtrusionBarriers::ExtrusionBarriers(std::initializer_list<ExtrusionBarrier> barriers,
                                     std::initializer_list<State> states, bool sort_barriers)
    : ExtrusionBarriers(barriers.size()) {
  assert(size() == states.size());
  _state = states;

  for (usize i = 0; i < barriers.size(); ++i) {
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    const auto& b = *(std::begin(barriers) + static_cast<isize>(i));

    _pos[i] = b.pos;
    _direction[i] = b.blocking_direction;
    _stp_active[i] = b.stp_active;
    _stp_inactive[i] = b.stp_inactive;
  }

  if (sort_barriers) {
    sort();
  }
}

usize ExtrusionBarriers::size() const noexcept {
  assert_buffer_sizes_are_equal();
  return _pos.size();
}

bool ExtrusionBarriers::empty() const noexcept { return size() == 0; }

void ExtrusionBarriers::resize(usize new_size) {
  _pos.resize(new_size);
  _direction.resize(new_size, dna::NONE);
  _stp_active.resize(new_size);
  _stp_inactive.resize(new_size);
  _state.resize(new_size);

  assert_buffer_sizes_are_equal();
}

void ExtrusionBarriers::clear() noexcept {
  _pos.clear();
  _direction.clear();
  _stp_active.clear();
  _stp_inactive.clear();
  _state.clear();

  assert_buffer_sizes_are_equal();
}

bp_t ExtrusionBarriers::pos(usize i) const noexcept {
  assert_index_within_bounds(i);
  return _pos[i];
}
const std::vector<bp_t>& ExtrusionBarriers::pos() const noexcept { return _pos; }

dna::Direction ExtrusionBarriers::direction(usize i) const noexcept {
  assert_index_within_bounds(i);
  return _direction[i];
}
const std::vector<dna::Direction>& ExtrusionBarriers::direction() const noexcept {
  return _direction;
}

double ExtrusionBarriers::stp_active(usize i) const noexcept {
  assert_index_within_bounds(i);
  return _stp_active[i]();
}
auto ExtrusionBarriers::stp_active() const noexcept -> const std::vector<TP>& {
  return _stp_active;
}

double ExtrusionBarriers::stp_inactive(usize i) const noexcept {
  assert_index_within_bounds(i);
  return _stp_inactive[i]();
}
auto ExtrusionBarriers::stp_inactive() const noexcept -> const std::vector<TP>& {
  return _stp_inactive;
}

auto ExtrusionBarriers::state(usize i) const noexcept -> State {
  assert_index_within_bounds(i);
  return _state[i];
}
auto ExtrusionBarriers::state() const noexcept -> const std::vector<State>& { return _state; }

bool ExtrusionBarriers::is_active(usize i) const noexcept {
  assert_index_within_bounds(i);
  return _state[i] == State::ACTIVE;
}

bool ExtrusionBarriers::is_not_active(usize i) const noexcept {
  assert_index_within_bounds(i);
  return _state[i] == State::INACTIVE;
}

double ExtrusionBarriers::occupancy(usize i) const noexcept {
  assert_index_within_bounds(i);
  return ExtrusionBarrier::compute_occupancy_from_stp(stp_active(i), stp_inactive(i));
}

auto ExtrusionBarriers::next_state(usize i, random::PRNG_t& rand_eng) noexcept -> State {
  const auto p = random::generate_canonical<double, std::numeric_limits<double>::digits>(rand_eng);

  if (is_not_active(i) && p > stp_inactive(i)) {
    _state[i] = State::ACTIVE;
  } else if (is_active(i) && p > stp_active(i)) {
    _state[i] = State::INACTIVE;
  }

  return _state[i];
}

void ExtrusionBarriers::next_state(random::PRNG_t& rand_eng) noexcept {
  for (usize i = 0; i < size(); ++i) {
    next_state(i, rand_eng);
  }
}

usize ExtrusionBarriers::count_active() const noexcept {
  return static_cast<usize>(std::count(_state.begin(), _state.end(), State::ACTIVE));
}
usize ExtrusionBarriers::count_inactive() const noexcept {
  const auto num_active =
      static_cast<usize>(std::count(_state.begin(), _state.end(), State::ACTIVE));
  assert(num_active <= size());
  return size() - num_active;
}

void ExtrusionBarriers::assert_buffer_sizes_are_equal() const noexcept {
  assert(_pos.size() == _direction.size());
  assert(_pos.size() == _stp_active.size());
  assert(_pos.size() == _stp_inactive.size());
  assert(_pos.size() == _state.size());
}

void ExtrusionBarriers::assert_index_within_bounds([[maybe_unused]] usize i) const noexcept {
  assert(i < size());
}

void ExtrusionBarriers::set(usize i, bp_t pos, dna::Direction direction, TP stp_active,
                            TP stp_inactive, State state) noexcept {
  assert_index_within_bounds(i);
  assert(direction != dna::NONE);
  _pos[i] = pos;
  _direction[i] = direction;
  _stp_active[i] = stp_active;
  _stp_inactive[i] = stp_inactive;
  _state[i] = state;
}

void ExtrusionBarriers::push_back(bp_t pos, dna::Direction direction, TP stp_active,
                                  TP stp_inactive, State state) {
  assert(direction != dna::NONE);
  _pos.push_back(pos);
  _direction.push_back(direction);
  _stp_active.push_back(stp_active);
  _stp_inactive.push_back(stp_inactive);
  _state.push_back(state);
}

void ExtrusionBarriers::set(usize i, const ExtrusionBarrier& barrier, State state) noexcept {
  set(i, barrier.pos, barrier.blocking_direction, barrier.stp_active, barrier.stp_inactive, state);
}

void ExtrusionBarriers::push_back(const ExtrusionBarrier& barrier, State state) {
  push_back(barrier.pos, barrier.blocking_direction, barrier.stp_active, barrier.stp_inactive,
            state);
}

void ExtrusionBarriers::set(usize i, State state) noexcept {
  assert_index_within_bounds(i);
  _state[i] = state;
}

auto ExtrusionBarriers::init_state(usize i, random::PRNG_t& rand_eng) noexcept -> State {
  assert_index_within_bounds(i);
  _state[i] = random::bernoulli_trial{occupancy(i)}(rand_eng) ? State::ACTIVE : State::INACTIVE;

  return _state[i];
}

void ExtrusionBarriers::init_states(random::PRNG_t& rand_eng) noexcept {
  for (usize i = 0; i < size(); ++i) {
    init_state(i, rand_eng);
  }
}

void ExtrusionBarriers::sort() {
  std::vector<usize> buff(size());
  sort(buff);
}

void ExtrusionBarriers::sort(std::vector<usize>& idx_buff) {
  idx_buff.resize(size());
  std::iota(idx_buff.begin(), idx_buff.end(), usize(0));

  cppsort::pdq_sort(idx_buff.begin(), idx_buff.end(),
                    [&](const auto i1, const auto i2) { return pos(i1) < pos(i2); });
  for (usize i = 0; i < size(); ++i) {
    // https://stackoverflow.com/a/22218699
    while (idx_buff[i] != i) {
      const auto j = idx_buff[i];
      const auto k = idx_buff[j];

      std::swap(_pos[j], _pos[k]);
      std::swap(_direction[j], _direction[k]);
      std::swap(_stp_active[j], _stp_active[k]);
      std::swap(_stp_inactive[j], _stp_inactive[k]);
      std::swap(_state[j], _state[k]);
      std::swap(idx_buff[i], idx_buff[j]);
    }
  }
}

}  // namespace modle
