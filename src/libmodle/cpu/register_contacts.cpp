// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

// clang-format off
#include "modle/simulation.hpp"
// clang-format on

#include <absl/types/span.h>                               // for Span

#include <algorithm>                                       // for min, all_of
#include <cassert>                                         // for assert

#include "modle/common/common.hpp"                         // for usize, bp_t, contacts_t, MODLE...
#include "modle/common/genextreme_value_distribution.hpp"  // for genextreme_value_distribution
#include "modle/common/random.hpp"                         // for PRNG_t, uniform_int_distribution
#include "modle/contact_matrix_dense.hpp"                  // for ContactMatrixDense
#include "modle/extrusion_factors.hpp"                     // for Lef, ExtrusionUnit
#include "modle/genome.hpp"                                // for GenomicInterval

namespace modle {

/// Check whether a LEF extrusion units fall in the interval (start_pos and end_pos)
[[nodiscard]] static constexpr bool lef_within_bound(const Lef& lef, const bp_t start_pos,
                                                     const bp_t end_pos) {
  assert(lef.is_bound());
  assert(end_pos >= start_pos);
  return MODLE_LIKELY(lef.rev_unit.pos() > start_pos && lef.rev_unit.pos() < end_pos &&
                      lef.fwd_unit.pos() > start_pos && lef.fwd_unit.pos() < end_pos);
}

/// Check whether a pait of genomic coordinates (FP) fall in the interval [start_pos and end_pos)
[[nodiscard]] static constexpr bool pos_within_bound(const double p1, const double p2,
                                                     const bp_t start_pos, const bp_t end_pos) {
  assert(start_pos <= end_pos);

  const auto start_pos_ = static_cast<double>(start_pos);
  const auto end_pos_ = static_cast<double>(end_pos);

  return MODLE_LIKELY(p1 >= start_pos_ && p2 >= start_pos_ && p1 < end_pos_ && p2 < end_pos_);
}

struct PosPair {
  double p1{};
  double p2{};
};

[[nodiscard]] static constexpr PosPair randomize_extrusion_unit_positions(
    const Lef& lef, const double mu, const double sigma, const double xi,
    const bool randomize_contacts, random::PRNG_t& rand_eng) {
  auto noise_gen = [&]() constexpr noexcept {
    if (randomize_contacts) {
      return genextreme_value_distribution<double>{mu, sigma, xi}(rand_eng);
    }
    return 0.0;
  };

  // We are performing most operations using double to deal with the possibility that the
  // noise generated to compute p1 is larger than the pos of the rev unit
  const auto [p1, p2] = std::minmax({static_cast<double>(lef.rev_unit.pos()) - noise_gen(),
                                     static_cast<double>(lef.fwd_unit.pos()) + noise_gen()});

  return {p1, p2};
}

[[nodiscard]] static const Lef& sample_lef_with_replacement(const absl::Span<const Lef> lefs,
                                                            random::PRNG_t& rand_eng) noexcept {
  assert(!lefs.empty());
  const auto i = random::uniform_int_distribution<usize>{0, lefs.size() - 1}(rand_eng);
  return lefs[i];
}

[[nodiscard]] static usize compute_num_contacts_loop(const usize num_contacts,
                                                     const double tad_to_loop_contact_ratio,
                                                     random::PRNG_t& rand_eng) {
  // Handle special case where TAD contact sampling has been disabled
  if (tad_to_loop_contact_ratio == 0) {
    return num_contacts;
  }

  // Handle special case where loop contact sampling has been disabled
  if (!std::isfinite(tad_to_loop_contact_ratio)) {
    assert(!std::isnan(tad_to_loop_contact_ratio));
    return 0;
  }

  const auto prob_loop_contact = 1.0 / (tad_to_loop_contact_ratio + 1.0);
  // Using usize as template param the distribution causes an ambiguous call to abs() inside
  // boost::random (boost v1.79)
  return static_cast<usize>(
      random::binomial_distribution<isize>{isize(num_contacts), prob_loop_contact}(rand_eng));
}

void Simulation::sample_and_register_contacts(State& s, usize num_sampling_events) const {
  assert(s.num_active_lefs == s.num_lefs);

  // Ensure we do not overshoot the target contact density
  if (this->target_contact_density > 0.0) {
    num_sampling_events = std::min(num_sampling_events, s.num_target_contacts - s.num_contacts);
  }

  // This can happen when simulating a huge number of cells with a low target contact density
  if (num_sampling_events == 0) {
    return;
  }

  const auto num_loop_contacts =
      compute_num_contacts_loop(num_sampling_events, this->tad_to_loop_contact_ratio, s.rand_eng);
  const auto num_tad_contacts = num_sampling_events - num_loop_contacts;

  assert(s.interval);
  s.num_contacts +=
      this->register_contacts_loop(*s.interval, s.get_lefs(), num_loop_contacts, s.rand_eng);
  s.num_contacts +=
      this->register_contacts_tad(*s.interval, s.get_lefs(), num_tad_contacts, s.rand_eng);

  if (this->track_1d_lef_position) {
    this->register_1d_lef_occupancy(*s.interval, s.get_lefs(), num_sampling_events, s.rand_eng);
  }

  assert(s.num_contacts <= s.num_target_contacts);
}

usize Simulation::register_contacts_loop(const bp_t start_pos, const bp_t end_pos,
                                         ContactMatrixDense<contacts_t>& contacts,
                                         const absl::Span<const Lef> lefs,
                                         usize num_sampling_events,
                                         random::PRNG_t& rand_eng) const {
  if (num_sampling_events == 0) {
    return 0;
  }

  using CS = ContactSamplingStrategy;
  assert(this->contact_sampling_strategy & CS::loop);

  // NOLINTNEXTLINE(readability-implicit-bool-conversion)
  const bool noisify_contacts = this->contact_sampling_strategy & CS::noisify;

  usize num_contacts_registered = 0;
  while (num_contacts_registered != num_sampling_events) {
    const auto& lef = sample_lef_with_replacement(lefs, rand_eng);
    if (MODLE_LIKELY(lef.is_bound() && lef_within_bound(lef, start_pos, end_pos))) {
      const auto [p1, p2] =
          randomize_extrusion_unit_positions(lef, this->genextreme_mu, this->genextreme_sigma,
                                             this->genextreme_xi, noisify_contacts, rand_eng);

      if (!pos_within_bound(p1, p2, start_pos, end_pos)) {
        continue;
      }

      const auto pos1 = static_cast<bp_t>(p1) - start_pos;
      const auto pos2 = static_cast<bp_t>(p2) - start_pos;
      contacts.increment(pos1 / this->bin_size, pos2 / this->bin_size);
      ++num_contacts_registered;
    }
  }

  return num_contacts_registered;
}

usize Simulation::register_contacts_tad(bp_t start_pos, bp_t end_pos,
                                        ContactMatrixDense<contacts_t>& contacts,
                                        absl::Span<const Lef> lefs, usize num_sampling_events,
                                        random::PRNG_t& rand_eng) const {
  if (num_sampling_events == 0) {
    return 0;
  }

  using CS = ContactSamplingStrategy;
  assert(this->contact_sampling_strategy & CS::tad);

  // NOLINTNEXTLINE(readability-implicit-bool-conversion)
  const bool noisify_contacts = this->contact_sampling_strategy & CS::noisify;

  usize num_contacts_registered = 0;
  while (num_contacts_registered != num_sampling_events) {
    const auto& lef = sample_lef_with_replacement(lefs, rand_eng);
    if (MODLE_LIKELY(lef.is_bound() && lef_within_bound(lef, start_pos, end_pos))) {
      const auto [p1, p2] =
          randomize_extrusion_unit_positions(lef, this->genextreme_mu, this->genextreme_sigma,
                                             this->genextreme_xi, noisify_contacts, rand_eng);

      if (!pos_within_bound(p1, p2, start_pos, end_pos)) {
        continue;
      }
      const auto p11 = random::uniform_int_distribution<bp_t>{static_cast<bp_t>(p1),
                                                              static_cast<bp_t>(p2)}(rand_eng);
      const auto p22 = random::uniform_int_distribution<bp_t>{static_cast<bp_t>(p1),
                                                              static_cast<bp_t>(p2)}(rand_eng);

      const auto pos1 = static_cast<bp_t>(p11) - start_pos;
      const auto pos2 = static_cast<bp_t>(p22) - start_pos;
      contacts.increment(pos1 / this->bin_size, pos2 / this->bin_size);
      ++num_contacts_registered;
    }
  }
  return num_contacts_registered;
}

usize Simulation::register_1d_lef_occupancy(bp_t start_pos, bp_t end_pos,
                                            std::vector<std::atomic<u64>>& occupancy_buff,
                                            absl::Span<const Lef> lefs, usize num_sampling_events,
                                            random::PRNG_t& rand_eng) const {
  if (num_sampling_events == 0) {
    return 0;
  }

  using CS = ContactSamplingStrategy;
  // NOLINTNEXTLINE(readability-implicit-bool-conversion)
  const bool noisify_positions = this->contact_sampling_strategy & CS::noisify;

  usize num_successful_sampling_events = 0;
  while (num_successful_sampling_events != num_sampling_events) {
    const auto& lef = sample_lef_with_replacement(lefs, rand_eng);
    if (MODLE_LIKELY(lef.is_bound() && lef_within_bound(lef, start_pos, end_pos))) {
      const auto [p1, p2] =
          randomize_extrusion_unit_positions(lef, this->genextreme_mu, this->genextreme_sigma,
                                             this->genextreme_xi, noisify_positions, rand_eng);

      if (!pos_within_bound(p1, p2, start_pos, end_pos)) {
        continue;
      }

      const auto i1 = utils::conditional_static_cast<usize>((static_cast<bp_t>(p1) - start_pos) /
                                                            this->bin_size);
      const auto i2 = utils::conditional_static_cast<usize>((static_cast<bp_t>(p2) - start_pos) /
                                                            this->bin_size);
      occupancy_buff[i1]++;
      occupancy_buff[i2]++;
      ++num_successful_sampling_events;
    }
  }
  return num_successful_sampling_events;
}

usize Simulation::register_contacts_loop(GenomicInterval& interval,
                                         const absl::Span<const Lef> lefs,
                                         usize num_sampling_events,
                                         random::PRNG_t& rand_eng) const {
  return this->register_contacts_loop(interval.start() + 1, interval.end() - 1, interval.contacts(),
                                      lefs, num_sampling_events, rand_eng);
}

usize Simulation::register_contacts_tad(GenomicInterval& interval, absl::Span<const Lef> lefs,
                                        usize num_sampling_events, random::PRNG_t& rand_eng) const {
  return this->register_contacts_tad(interval.start() + 1, interval.end() - 1, interval.contacts(),
                                     lefs, num_sampling_events, rand_eng);
}

usize Simulation::register_1d_lef_occupancy(GenomicInterval& interval, absl::Span<const Lef> lefs,
                                            usize num_sampling_events,
                                            random::PRNG_t& rand_eng) const {
  return this->register_1d_lef_occupancy(interval.start() + 1, interval.end() - 1,
                                         interval.lef_1d_occupancy(), lefs, num_sampling_events,
                                         rand_eng);
}

}  // namespace modle
