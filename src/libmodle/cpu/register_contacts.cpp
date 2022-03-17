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
#include "modle/common/random_sampling.hpp"                // for random_sample
#include "modle/contacts.hpp"                              // for ContactMatrix
#include "modle/extrusion_factors.hpp"                     // for Lef, ExtrusionUnit
#include "modle/genome.hpp"                                // for Chromosome

namespace modle {

void Simulation::sample_and_register_contacts(State& s, const double avg_nlefs_to_sample) const {
  assert(s.num_active_lefs == s.num_lefs);
  auto nlefs_to_sample =
      std::min(s.num_lefs, random::poisson_distribution<usize>{avg_nlefs_to_sample}(s.rand_eng));

  // Select LEFs to be used for contact registration.
  // We are sampling from fwd ranks to avoid having to allocate a vector of indices just to
  // do this sampling
  const auto lef_idx = s.get_idx_buff(nlefs_to_sample);
  random_sample(s.get_fwd_ranks().begin(), s.get_fwd_ranks().end(), lef_idx.begin(), lef_idx.size(),
                s.rand_eng);
#ifndef NDEBUG  // GCC 9.5.0 chokes if we use if constexpr (utils::ndebug_not_defined()) here
  if (!std::all_of(lef_idx.begin(), lef_idx.end(), [&](const auto i) { return i < s.num_lefs; })) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("lef_idx.size()={}; num_lefs={};\nlef_idx=[{}]\n"), lef_idx.size(),
                    s.num_lefs, fmt::join(lef_idx, ", ")));
  }
#endif
  const auto start_pos = s.is_modle_sim_state() ? s.chrom->start_pos() : s.window_start;
  const auto end_pos = s.is_modle_sim_state() ? s.chrom->end_pos() : s.window_end;

  assert(s.contacts);
  s.num_contacts +=
      this->register_contacts_loop(start_pos + 1, end_pos - 1, *s.contacts, s.get_lefs(), lef_idx,
                                   s.num_target_contacts - s.num_contacts, s.rand_eng);
  s.num_contacts +=
      this->register_contacts_tad(start_pos + 1, end_pos - 1, *s.contacts, s.get_lefs(), lef_idx,
                                  s.num_target_contacts - s.num_contacts, s.rand_eng);
  assert(s.num_contacts <= s.num_target_contacts);
}

usize Simulation::register_contacts_loop(const bp_t start_pos, const bp_t end_pos,
                                         ContactMatrix<contacts_t>& contacts,
                                         const absl::Span<const Lef> lefs,
                                         const absl::Span<const usize> lef_idx,
                                         const usize max_num_contacts_to_register,
                                         random::PRNG_t& rand_eng) const {
  // Register contacts for the selected LEFs (excluding LEFs that have one of their units at the
  // beginning/end of a chromosome)

  using CS = ContactSamplingStrategy;
  // NOLINTNEXTLINE(readability-implicit-bool-conversion)
  if (!(this->contact_sampling_strategy & CS::loop) ||
      this->number_of_loop_contacts_per_sampling_event == 0) {
    return 0;
  }

  // NOLINTNEXTLINE(readability-implicit-bool-conversion)
  const bool noisify_contacts = this->contact_sampling_strategy & CS::noisify;

  auto noise_gen = [&]() {
    if (noisify_contacts) {
      return genextreme_value_distribution<double>{this->genextreme_mu, this->genextreme_sigma,
                                                   this->genextreme_xi}(rand_eng);
    }
    return 0.0;
  };

  const auto start_pos_dbl = static_cast<double>(start_pos);
  const auto end_pos_dbl = static_cast<double>(end_pos);

  usize new_contacts = 0;
  for (const auto i : lef_idx) {
    assert(i < lefs.size());
    const auto& lef = lefs[i];
    if (MODLE_LIKELY(lef.is_bound() && lef.rev_unit.pos() > start_pos &&
                     lef.rev_unit.pos() < end_pos && lef.fwd_unit.pos() > start_pos &&
                     lef.fwd_unit.pos() < end_pos)) {
      for (usize j = 0; j < this->number_of_loop_contacts_per_sampling_event; ++j) {
        if (MODLE_UNLIKELY(new_contacts == max_num_contacts_to_register)) {
          return new_contacts;
        }
        // We are performing most operations using double to deal with the possibility that the
        // noise generated to compute p1 is larger than the pos of the rev unit
        const auto [p1, p2] = std::minmax({static_cast<double>(lef.rev_unit.pos()) - noise_gen(),
                                           static_cast<double>(lef.fwd_unit.pos()) + noise_gen()});

        if (MODLE_UNLIKELY(p1 < start_pos_dbl || p2 < start_pos_dbl || p1 >= end_pos_dbl ||
                           p2 >= end_pos_dbl)) {
          continue;
        }

        const auto pos1 = static_cast<bp_t>(p1) - start_pos;
        const auto pos2 = static_cast<bp_t>(p2) - start_pos;
        contacts.increment(pos1 / this->bin_size, pos2 / this->bin_size);
        ++new_contacts;
      }
    }
  }
  return new_contacts;
}

usize Simulation::register_contacts_tad(bp_t start_pos, bp_t end_pos,
                                        ContactMatrix<contacts_t>& contacts,
                                        absl::Span<const Lef> lefs,
                                        const absl::Span<const usize> lef_idx,
                                        const usize max_num_contacts_to_register,
                                        random::PRNG_t& rand_eng) const {
  // Register contacts for the selected LEFs (excluding LEFs that have one of their units at the
  // beginning/end of a chromosome)

  using CS = ContactSamplingStrategy;
  // NOLINTNEXTLINE(readability-implicit-bool-conversion)
  if (!(this->contact_sampling_strategy & CS::tad) ||
      this->number_of_loop_contacts_per_sampling_event == 0) {
    return 0;
  }

  // NOLINTNEXTLINE(readability-implicit-bool-conversion)
  const bool noisify_contacts = this->contact_sampling_strategy & CS::noisify;

  auto noise_gen = [&]() {
    if (noisify_contacts) {
      return genextreme_value_distribution<double>{this->genextreme_mu, this->genextreme_sigma,
                                                   this->genextreme_xi}(rand_eng);
    }
    return 0.0;
  };

  const auto start_pos_dbl = static_cast<double>(start_pos);
  const auto end_pos_dbl = static_cast<double>(end_pos);

  usize new_contacts = 0;
  for (const auto i : lef_idx) {
    assert(i < lefs.size());
    const auto& lef = lefs[i];
    if (MODLE_LIKELY(lef.is_bound() && lef.rev_unit.pos() > start_pos &&
                     lef.rev_unit.pos() < end_pos && lef.fwd_unit.pos() > start_pos &&
                     lef.fwd_unit.pos() < end_pos)) {
      // We are performing most operations using double to deal with the possibility that the
      // noise generated to compute p1 is larger than the pos of the rev unit
      const auto [p1, p2] = std::minmax({static_cast<double>(lef.rev_unit.pos()) - noise_gen(),
                                         static_cast<double>(lef.fwd_unit.pos()) + noise_gen()});
      if (MODLE_UNLIKELY(p1 < start_pos_dbl || p2 < start_pos_dbl || p1 >= end_pos_dbl ||
                         p2 >= end_pos_dbl)) {
        continue;
      }
      for (usize j = 0; j < this->number_of_tad_contacts_per_sampling_event; ++j) {
        if (MODLE_UNLIKELY(new_contacts == max_num_contacts_to_register)) {
          return new_contacts;
        }
        const auto p11 = random::uniform_int_distribution<bp_t>{static_cast<bp_t>(p1),
                                                                static_cast<bp_t>(p2)}(rand_eng);
        const auto p22 = random::uniform_int_distribution<bp_t>{static_cast<bp_t>(p1),
                                                                static_cast<bp_t>(p2)}(rand_eng);

        const auto pos1 = static_cast<bp_t>(p11) - start_pos;
        const auto pos2 = static_cast<bp_t>(p22) - start_pos;
        contacts.increment(pos1 / this->bin_size, pos2 / this->bin_size);
        ++new_contacts;
      }
    }
  }
  return new_contacts;
}

usize Simulation::register_contacts_loop(Chromosome& chrom, const absl::Span<const Lef> lefs,
                                         const absl::Span<const usize> selected_lef_idx,
                                         const usize max_num_contacts_to_register,
                                         random::PRNG_t& rand_eng) const {
  return this->register_contacts_loop(chrom.start_pos() + 1, chrom.end_pos() - 1, chrom.contacts(),
                                      lefs, selected_lef_idx, max_num_contacts_to_register,
                                      rand_eng);
}

usize Simulation::register_contacts_tad(Chromosome& chrom, absl::Span<const Lef> lefs,
                                        absl::Span<const usize> selected_lef_idx,
                                        const usize max_num_contacts_to_register,
                                        random::PRNG_t& rand_eng) const {
  return this->register_contacts_tad(chrom.start_pos() + 1, chrom.end_pos() - 1, chrom.contacts(),
                                     lefs, selected_lef_idx, max_num_contacts_to_register,
                                     rand_eng);
}

}  // namespace modle
