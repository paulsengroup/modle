#include "modle/genome.hpp"

#include <cstdlib>
#include <execution>
#include <filesystem>
#include <functional>

#include "absl/strings/str_format.h"
#include "modle/extr_barrier.hpp"
#include "modle/parsers.hpp"

namespace modle {

Genome::Genome(const config& c)
    : _seed(c.seed),
      _rand_eng(std::mt19937(_seed)),
      _path_to_chr_size_file(c.path_to_chr_sizes),
      _bin_size(c.bin_size),
      _avg_lef_processivity(c.average_lef_processivity),
      _probability_of_barrier_block(c.probability_of_extrusion_barrier_block),
      _probability_of_lef_rebind(c.probability_of_lef_rebind),
      _probability_of_extr_unit_bypass(c.probability_of_extrusion_unit_bypass),
      _lef_unloader_strength_coeff(c.lef_unloader_strength),
      _lefs(generate_lefs(c.number_of_lefs)),
      _chromosomes(init_chromosomes_from_file(c.diagonal_width)),
      _sampling_interval(c.contact_sampling_interval),
      _randomize_contact_sampling(c.randomize_contact_sampling),
      _sample_contacts(1.0 / _sampling_interval) {}

std::vector<uint32_t> Genome::get_chromosome_lengths() const {
  std::vector<uint32_t> lengths;
  lengths.reserve(this->get_n_chromosomes());
  for (const auto& chr : this->_chromosomes) lengths.push_back(chr.length());
  return lengths;
}

std::vector<std::string_view> Genome::get_chromosome_names() const {
  std::vector<std::string_view> names;
  names.reserve(this->get_n_chromosomes());
  for (const auto& chr : this->_chromosomes) names.push_back(chr.name);
  return names;
}

uint32_t Genome::get_n_lefs() const { return this->_lefs.size(); }

uint32_t Genome::get_n_of_free_lefs() const {
  return std::accumulate(
      this->_lefs.begin(), this->_lefs.end(), 0UL,
      [](uint32_t accumulator, const Lef& lef) { return accumulator + !lef.is_bound(); });
}

uint32_t Genome::get_n_of_busy_lefs() const {
  return this->get_n_lefs() - this->get_n_of_free_lefs();
}

void Genome::write_contacts_to_file(const std::string& output_dir, bool force_overwrite) const {
  std::filesystem::create_directories(output_dir);
  std::for_each(
      std::execution::par, this->_chromosomes.begin(), this->_chromosomes.end(),
      [&](const Chromosome& chr) {
        auto t0 = absl::Now();
        auto path_to_outfile = std::filesystem::weakly_canonical(
            absl::StrFormat("%s/%s.tsv.bz2", output_dir, chr.name));
        if (!force_overwrite && std::filesystem::exists(path_to_outfile)) {
          absl::FPrintF(stderr,
                        "File '%s' already exists. Pass --force to overwrite... SKIPPING.\n",
                        path_to_outfile);
          return;
        }
        absl::FPrintF(stderr, "Writing full contact matrix for '%s' to file '%s'...\n", chr.name,
                      path_to_outfile);
        {
          auto [bytes_in, bytes_out] = chr.contacts.write_full_matrix_to_tsv(path_to_outfile);
          absl::FPrintF(
              stderr,
              "DONE writing '%s' in %s! Compressed size: %.2f MB (compression ratio %.2fx)\n",
              path_to_outfile, absl::FormatDuration(absl::Now() - t0), bytes_out / 1.0e6,
              static_cast<double>(bytes_in) / bytes_out);
          t0 = absl::Now();
        }
        path_to_outfile = std::filesystem::weakly_canonical(
            absl::StrFormat("%s/%s_raw.tsv.bz2", output_dir, chr.name));
        absl::FPrintF(stderr, "Writing raw contact matrix for '%s' to file '%s'...\n", chr.name,
                      path_to_outfile);
        auto [bytes_in, bytes_out] = chr.contacts.write_to_tsv(path_to_outfile);
        absl::FPrintF(
            stderr, "DONE writing '%s' in %s! Compressed size: %.2f MB (compression ratio %.2fx)\n",
            path_to_outfile, absl::FormatDuration(absl::Now() - t0), bytes_out / 1.0e6,
            static_cast<double>(bytes_in) / bytes_out);
      });
}

void Genome::make_heatmaps(std::string_view output_dir, bool force_overwrite,
                           const std::string& script) const {
  const auto t0 = absl::Now();
  absl::FPrintF(stderr, "Generating heatmaps for %lu chromosomes...", this->get_n_chromosomes());
  const auto cmd =
      absl::StrFormat("%s --input-dir %s --output-dir %s --bin-size %lu%s", script, output_dir,
                      output_dir, this->_bin_size, force_overwrite ? " --force" : "");
  if (auto status = std::system(cmd.c_str()); status != 0) {
    if (WIFEXITED(status) == 0) {
      if (auto ec = WEXITSTATUS(status); ec != 0) {
        throw std::runtime_error(
            absl::StrFormat("Command '%s' terminated with exit code %d.", cmd, ec));
      }
    }
  }
  absl::FPrintF(stderr, "DONE! Plotting took %s.\n", absl::FormatDuration(absl::Now() - t0));
}

std::vector<Chromosome> Genome::init_chromosomes_from_file(uint32_t diagonal_width) const {
  std::vector<Chromosome> chromosomes;
  auto parser = ChrSizeParser(this->_path_to_chr_size_file);
  for (const auto& chr : parser.parse()) {
    chromosomes.emplace_back(chr.name, chr.size, this->_bin_size, diagonal_width);
  }
  return chromosomes;
}

std::vector<Lef> Genome::generate_lefs(uint32_t n) {
  std::vector<Lef> v;
  v.reserve(n);
  for (auto i = 0UL; i < n; ++i)
    v.emplace_back(this->_bin_size, this->_avg_lef_processivity,
                   this->_probability_of_extr_unit_bypass, this->_lef_unloader_strength_coeff);
  return v;
}

uint64_t Genome::n50() const {
  auto chr_lengths = this->get_chromosome_lengths();
  uint64_t n50_thresh = this->size() / 2;
  std::sort(chr_lengths.rbegin(), chr_lengths.rend(), std::greater<>());
  uint64_t tot = 0;
  for (const auto& len : chr_lengths) {
    if ((tot += len) >= n50_thresh) return len;
  }
  return chr_lengths.back();
}

void Genome::randomly_generate_extrusion_barriers(uint32_t n_barriers) {
  const auto& weights = this->get_chromosome_lengths();
  std::discrete_distribution<> chr_idx(weights.begin(), weights.end());
  std::bernoulli_distribution strand_selector(0.5);

  for (auto i = 0UL; i < n_barriers; ++i) {
    // Randomly select a chromosome, barrier binding pos and direction
    auto& chr = this->_chromosomes[chr_idx(this->_rand_eng)];
    std::uniform_int_distribution<uint64_t> uniform_rng(0, chr.length());
    const auto barrier_position = uniform_rng(this->_rand_eng);
    const DNA::Direction direction =
        strand_selector(this->_rand_eng) ? DNA::Direction::rev : DNA::Direction::fwd;

    // Add the new extrusion barrier to the appropriate bin
    auto& bin = chr.dna.get_bin_from_pos(barrier_position);
    bin.add_extr_barrier(barrier_position, this->_probability_of_barrier_block, direction);
  }
  for (auto& chr : this->get_chromosomes()) chr.sort_barriers_by_pos();
}

std::pair<uint64_t, uint64_t> Genome::import_extrusion_barriers_from_bed(
    std::string_view path_to_bed, double probability_of_block) {
  auto p = modle::BEDParser(path_to_bed, BED::Standard::BED6);
  uint64_t nrecords = 0;
  uint64_t nrecords_ignored = 0;
  absl::flat_hash_map<std::string_view, DNA*> chromosomes;
  chromosomes.reserve(this->get_n_chromosomes());
  for (auto& chr : this->get_chromosomes()) chromosomes.emplace(chr.name, &chr.dna);
  for (auto& record : p.parse_all()) {
    ++nrecords;
    if (!chromosomes.contains(record.chrom)) {
      ++nrecords_ignored;
      continue;
    }
    if (record.strand == '.') continue;  // TODO: Figure out if this is an OK thing to do
    if (probability_of_block != 0) record.score = probability_of_block;
    if (record.score < 0 || record.score > 1) {
      throw std::runtime_error(
          absl::StrFormat("Invalid score field detected for record %s[%lu-%lu]: expected a score "
                          "between 0 and 1, got %.4g.",
                          record.name, record.chrom_start, record.chrom_end, record.score));
    }
    chromosomes.at(record.chrom)->add_extr_barrier(record);
  }
  return {nrecords, nrecords_ignored};
}

const std::vector<Chromosome>& Genome::get_chromosomes() const { return this->_chromosomes; }
std::vector<Chromosome>& Genome::get_chromosomes() { return this->_chromosomes; }

std::vector<Lef>& Genome::get_lefs() { return this->_lefs; }
const std::vector<Lef>& Genome::get_lefs() const { return this->_lefs; }

uint32_t Genome::get_n_chromosomes() const { return this->_chromosomes.size(); }

uint64_t Genome::size() const {
  return std::accumulate(
      this->_chromosomes.begin(), this->_chromosomes.end(), 0ULL,
      [&](uint64_t accumulator, const Chromosome& chr) { return accumulator + chr.length(); });
}

uint32_t Genome::get_n_bins() const {
  return std::accumulate(
      this->_chromosomes.begin(), this->_chromosomes.end(), 0UL,
      [](uint64_t accumulator, const Chromosome& chr) { return accumulator + chr.get_n_bins(); });
}

uint32_t Genome::get_n_barriers() const {
  return std::accumulate(this->_chromosomes.begin(), this->_chromosomes.end(), 0UL,
                         [](uint32_t accumulator, const Chromosome& chr) {
                           return accumulator + chr.get_n_barriers();
                         });
}

void Genome::randomly_bind_lefs() {
  const auto& weights = this->get_chromosome_lengths();
  std::discrete_distribution<> chr_idx(weights.begin(), weights.end());

  for (auto& lef : this->_lefs) {
    // Randomly select a chromosome, barrier binding pos and direction
    auto& chr = this->_chromosomes[chr_idx(this->_rand_eng)];
    lef.randomly_bind_to_chr(chr, this->_rand_eng);  // And bind to it
  }
}

uint32_t Genome::run_burnin(double prob_of_rebinding, uint16_t target_n_of_unload_events,
                            uint64_t min_extr_rounds) {
  const auto t0 = absl::Now();
  double avg_num_of_extr_events_per_bind =
      (static_cast<double>(this->_avg_lef_processivity) / this->_bin_size) /
      2 /* N of active extr. unit */;
  uint32_t n_of_lefs_to_bind_each_round =
      std::floor(this->get_n_lefs() / avg_num_of_extr_events_per_bind);
  n_of_lefs_to_bind_each_round += n_of_lefs_to_bind_each_round == 0;

  std::vector<uint16_t> unload_events(this->get_n_lefs(), 0);
  const auto& weights = this->get_chromosome_lengths();
  std::discrete_distribution<> chr_idx(weights.begin(), weights.end());

  uint32_t start_idx = 0, end_idx = n_of_lefs_to_bind_each_round;
  for (uint32_t rounds = 0;; ++rounds) {
    for (; start_idx < end_idx && start_idx < this->get_n_lefs(); ++start_idx) {
      auto& chr = this->_chromosomes[chr_idx(this->_rand_eng)];
      this->_lefs[start_idx].randomly_bind_to_chr(chr, this->_rand_eng);
    }
    for (uint32_t i = 0; i < start_idx; ++i) {
      auto& lef = this->_lefs[i];

      if (lef.is_bound()) {
        lef.extrude();
        unload_events[i] += !lef.is_bound();
      }
    }
    for (uint32_t i = 0; i < start_idx; ++i) {
      if (auto& lef = this->_lefs[i]; lef.is_bound()) {
        lef.check_constraints(this->_rand_eng);
      } else {
        auto& chr = this->_chromosomes[chr_idx(this->_rand_eng)];
        lef.try_rebind(chr, this->_rand_eng, prob_of_rebinding, false);
      }
    }

    start_idx = std::min(end_idx, this->get_n_lefs());
    end_idx += n_of_lefs_to_bind_each_round;
    //    absl::FPrintF(stderr, "%s\n", absl::StrJoin(unload_events.begin(), unload_events.end(), ",
    //    "));

    if (rounds >= min_extr_rounds &&
        std::all_of(unload_events.begin(), unload_events.end(),
                    [&](uint16_t n) { return n >= target_n_of_unload_events; })) {
      absl::FPrintF(stderr, "Burnin completed in %s! (%lu rounds).\n",
                    absl::FormatDuration(absl::Now() - t0), rounds);
      return rounds;
    }
  }
}

void Genome::simulate_extrusion(uint32_t iterations) {
  const auto step = iterations / 500;
  const auto& weights = this->get_chromosome_lengths();
  std::discrete_distribution<> chr_idx(weights.begin(), weights.end());

  auto t0 = absl::Now();
  for (auto i = 1UL; i <= iterations; ++i) {
    bool register_contacts = this->_randomize_contact_sampling
                                 ? this->_sample_contacts(this->_rand_eng)
                                 : i % this->_sampling_interval == 0;

    for (auto& lef : this->_lefs) {
      if (lef.is_bound()) {  // Register contact and extrude if LEF is bound
        if (register_contacts) lef.register_contact();
        lef.extrude();
      }
    }

    for (auto& lef : this->_lefs) {
      // Check whether the last round of extrusions caused a collision with another LEF or an
      // extr. boundary, apply the stall and/or increase LEF lifetime where appropriate
      if (lef.is_bound()) {
        lef.check_constraints(this->_rand_eng);
      } else {
        // Randomly select a chromosome and try to bind one of the free LEFs to it
        auto& chr = this->_chromosomes[chr_idx(this->_rand_eng)];
        lef.try_rebind(chr, this->_rand_eng, this->_probability_of_lef_rebind, register_contacts);
      }
    }

    if (i % step == 0) {
      absl::FPrintF(stderr, "Running iteration %lu/%lu (%.2f iterations/s)\n", i, iterations,
                    step / absl::ToDoubleSeconds(absl::Now() - t0));
      t0 = absl::Now();
    }
  }
}
}  // namespace modle