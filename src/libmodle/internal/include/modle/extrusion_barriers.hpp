// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <cassert>  // for assert
#include <initializer_list>
#include <vector>  // for vector

#include "modle/common/common.hpp"  // for bp_t
#include "modle/common/dna.hpp"     // for Direction
#include "modle/common/random.hpp"  // for PRNG_t
namespace modle {

namespace internal {
class TransitionProbability {
  double _p{0.0};

 public:
  constexpr TransitionProbability() = default;
  // NOLINTNEXTLINE(hicpp-explicit-conversions)
  constexpr TransitionProbability(double p) noexcept;
  [[nodiscard]] constexpr double operator()() const noexcept;
};
}  // namespace internal

struct ExtrusionBarrier;

class ExtrusionBarriers {
 public:
  enum class State : u8f { INACTIVE = 0, ACTIVE = 1 };
  using TP = internal::TransitionProbability;

 private:
  std::vector<bp_t> _pos{};
  std::vector<dna::Direction> _direction{};
  std::vector<TP> _stp_active{};
  std::vector<TP> _stp_inactive{};
  std::vector<State> _state{};

 public:
  ExtrusionBarriers() = default;
  explicit ExtrusionBarriers(usize size);
  template <class BarrierIt, class StateIt>
  ExtrusionBarriers(BarrierIt first_barrier, BarrierIt last_barrier, StateIt first_state,
                    bool sort = true);
  template <class BarrierIt>
  ExtrusionBarriers(BarrierIt first_barrier, BarrierIt last_barrier, State state = State::INACTIVE,
                    bool sort = true);
  ExtrusionBarriers(std::initializer_list<bp_t> pos,
                    std::initializer_list<dna::Direction> direction,
                    std::initializer_list<TP> stp_active, std::initializer_list<TP> stp_inactive,
                    std::initializer_list<State> state, bool sort = true);
  ExtrusionBarriers(std::initializer_list<ExtrusionBarrier> barriers,
                    std::initializer_list<State> states, bool sort = true);

  void set(usize i, bp_t pos, dna::Direction direction, TP stp_active, TP stp_inactive,
           State state = State::INACTIVE) noexcept;
  void push_back(bp_t pos, dna::Direction direction, TP stp_active, TP stp_inactive,
                 State state = State::INACTIVE);
  void set(usize i, const ExtrusionBarrier& barrier, State state = State::INACTIVE) noexcept;
  void push_back(const ExtrusionBarrier& barrier, State state = State::INACTIVE);
  void set(usize i, State state) noexcept;
  auto init_state(usize i, random::PRNG_t& rand_eng) noexcept -> State;
  void init_states(random::PRNG_t& rand_eng) noexcept;

  [[nodiscard]] usize size() const noexcept;
  [[nodiscard]] bool empty() const noexcept;
  void resize(usize new_size);
  void clear() noexcept;

  [[nodiscard]] bp_t pos(usize i) const noexcept;
  [[nodiscard]] dna::Direction block_direction(usize i) const noexcept;
  [[nodiscard]] double stp_active(usize i) const noexcept;
  [[nodiscard]] double stp_inactive(usize i) const noexcept;
  [[nodiscard]] auto state(usize i) const noexcept -> State;
  [[nodiscard]] bool is_active(usize i) const noexcept;
  [[nodiscard]] bool is_not_active(usize i) const noexcept;
  [[nodiscard]] double occupancy(usize i) const noexcept;

  [[nodiscard]] const std::vector<bp_t>& pos() const noexcept;
  [[nodiscard]] const std::vector<dna::Direction>& block_direction() const noexcept;
  [[nodiscard]] auto stp_active() const noexcept -> const std::vector<TP>&;
  [[nodiscard]] auto stp_inactive() const noexcept -> const std::vector<TP>&;
  [[nodiscard]] auto state() const noexcept -> const std::vector<State>&;

  auto next_state(usize i, random::PRNG_t& rand_eng) noexcept -> State;
  void next_state(random::PRNG_t& rand_eng) noexcept;

  [[nodiscard]] usize count_active() const noexcept;
  [[nodiscard]] usize count_inactive() const noexcept;

  void sort();
  void sort(std::vector<usize>& idx_buff);

 private:
  void assert_buffer_sizes_are_equal() const noexcept;
  void assert_index_within_bounds(usize i) const noexcept;
};

struct ExtrusionBarrier {
  using TP = internal::TransitionProbability;
  bp_t pos{(std::numeric_limits<bp_t>::max)()};
  TP stp_active{};
  TP stp_inactive{};
  dna::Direction blocking_direction{dna::NONE};

  constexpr ExtrusionBarrier(bp_t pos_, TP transition_prob_active_to_active,
                             TP transition_prob_inactive_to_inactive,
                             dna::Direction motif_direction);
  constexpr ExtrusionBarrier(bp_t pos_, TP transition_prob_active_to_active,
                             TP transition_prob_inactive_to_inactive, char motif_direction);

  [[nodiscard]] constexpr bool operator==(const ExtrusionBarrier& other) const noexcept;
  [[nodiscard]] constexpr bool operator!=(const ExtrusionBarrier& other) const noexcept;
  [[nodiscard]] constexpr bool operator<(const ExtrusionBarrier& other) const noexcept;
  [[nodiscard]] constexpr bool operator>(const ExtrusionBarrier& other) const noexcept;
  [[nodiscard]] constexpr bool operator>=(const ExtrusionBarrier& other) const noexcept;
  [[nodiscard]] constexpr bool operator<=(const ExtrusionBarrier& other) const noexcept;

  static constexpr double compute_stp_active_from_occupancy(TP stp_inactive,
                                                            double occupancy) noexcept;
  static constexpr double compute_occupancy_from_stp(TP stp_active, TP stp_inactive) noexcept;
};

namespace internal {
[[nodiscard]] constexpr dna::Direction char_to_strand(char d);
}  // namespace internal

}  // namespace modle

template <>
struct fmt::formatter<modle::ExtrusionBarrier> {
  char presentation = 's';  // s == short, f == full
  constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin());
  template <typename FormatContext>
  inline auto format(const modle::ExtrusionBarrier& b, FormatContext& ctx) -> decltype(ctx.out());
};

#include "../../extrusion_barriers_impl.hpp"
