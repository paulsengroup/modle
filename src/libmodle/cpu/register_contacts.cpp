// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

// clang-format off
#include "modle/simulation.hpp"
// clang-format on

#include <absl/types/span.h>  // for Span

#include <algorithm>  // for min, all_of
#include <cassert>    // for assert

#include "modle/common/common.hpp"                         // for usize, bp_t, contacts_t, MODLE...
#include "modle/common/genextreme_value_distribution.hpp"  // for genextreme_value_distribution
#include "modle/common/random.hpp"                         // for PRNG_t, uniform_int_distribution
#include "modle/contacts.hpp"                              // for ContactMatrix
#include "modle/extrusion_factors.hpp"                     // for Lef, ExtrusionUnit
#include "modle/genome.hpp"                                // for Chromosome

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

[[nodiscard]] static constexpr usize compute_num_contacts_loop(
    const usize num_contacts, const double tad_to_loop_contact_ratio) {
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
  const auto num_loop_contacts = std::round(prob_loop_contact * static_cast<double>(num_contacts));

  assert(static_cast<usize>(num_loop_contacts) <= num_contacts);
  return static_cast<usize>(num_loop_contacts);
}

void Simulation::sample_and_register_contacts(State& s, usize num_contacts_to_sample) const {
  assert(s.num_active_lefs == s.num_lefs);

  // Ensure we do not overshoot the target contact density
  if (this->target_contact_density > 0.0) {
    num_contacts_to_sample =
        std::min(num_contacts_to_sample, s.num_target_contacts - s.num_contacts);
  }

  // This can happen when simulating a huge number of cells with a low target contact density
  if (num_contacts_to_sample == 0) {
    return;
  }

  const auto start_pos = s.is_modle_sim_state() ? s.chrom->start_pos() : s.window_start;
  const auto end_pos = s.is_modle_sim_state() ? s.chrom->end_pos() : s.window_end;

  const auto num_loop_contacts =
      compute_num_contacts_loop(num_contacts_to_sample, this->tad_to_loop_contact_ratio);
  const auto num_tad_contacts = num_contacts_to_sample - num_loop_contacts;

  assert(s.contacts);
  this->register_contacts_loop(start_pos + 1, end_pos - 1, *s.contacts, s.get_lefs(),
                               num_loop_contacts, s.rand_eng);
  this->register_contacts_tad(start_pos + 1, end_pos - 1, *s.contacts, s.get_lefs(),
                              num_tad_contacts, s.rand_eng);

  s.num_contacts += num_contacts_to_sample;
  assert(s.num_contacts <= s.num_target_contacts);
}

void Simulation::register_contacts_loop(const bp_t start_pos, const bp_t end_pos,
                                        ContactMatrix<contacts_t>& contacts,
                                        const absl::Span<const Lef> lefs,
                                        usize num_contacts_to_register,
                                        random::PRNG_t& rand_eng) const {
  if (num_contacts_to_register == 0) {
    return;
  }

  using CS = ContactSamplingStrategy;
  assert(this->contact_sampling_strategy & CS::loop);

  // NOLINTNEXTLINE(readability-implicit-bool-conversion)
  const bool noisify_contacts = this->contact_sampling_strategy & CS::noisify;

  while (num_contacts_to_register != 0) {
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
      --num_contacts_to_register;
    }
  }
}

void Simulation::register_contacts_tad(bp_t start_pos, bp_t end_pos,
                                       ContactMatrix<contacts_t>& contacts,
                                       absl::Span<const Lef> lefs, usize num_contacts_to_register,
                                       random::PRNG_t& rand_eng) const {
  if (num_contacts_to_register == 0) {
    return;
  }

  using CS = ContactSamplingStrategy;
  assert(this->contact_sampling_strategy & CS::tad);

  // NOLINTNEXTLINE(readability-implicit-bool-conversion)
  const bool noisify_contacts = this->contact_sampling_strategy & CS::noisify;

  while (num_contacts_to_register != 0) {
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
      --num_contacts_to_register;
    }
  }
}

void Simulation::register_contacts_loop(Chromosome& chrom, const absl::Span<const Lef> lefs,
                                        usize num_contacts_to_register,
                                        random::PRNG_t& rand_eng) const {
  this->register_contacts_loop(chrom.start_pos() + 1, chrom.end_pos() - 1, chrom.contacts(), lefs,
                               num_contacts_to_register, rand_eng);
}

void Simulation::register_contacts_tad(Chromosome& chrom, absl::Span<const Lef> lefs,
                                       usize num_contacts_to_register,
                                       random::PRNG_t& rand_eng) const {
  this->register_contacts_tad(chrom.start_pos() + 1, chrom.end_pos() - 1, chrom.contacts(), lefs,
                              num_contacts_to_register, rand_eng);
}

}  // namespace modle
