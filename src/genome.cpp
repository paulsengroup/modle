#include "modle/genome.hpp"

#include <cstdlib>
#include <execution>
#include <filesystem>
#include <functional>

#include "absl/strings/str_format.h"
#include "modle/extr_barrier.hpp"
#include "modle/parsers.hpp"

namespace modle {

Genome::Genome(std::string_view path_to_bed, uint32_t bin_size, uint32_t n_lefs,
               uint32_t avg_lef_processivity, double probability_of_barrier_block,
               double probability_of_lef_rebind, double probability_of_extr_unit_bypass,
               uint64_t seed)
    : _path_to_bed(path_to_bed),
      _bin_size(bin_size),
      _avg_lef_processivity(avg_lef_processivity),
      _probability_of_barrier_block(probability_of_barrier_block),
      _probability_of_lef_rebind(probability_of_lef_rebind),
      _probability_of_extr_unit_bypass(probability_of_extr_unit_bypass),
      _lefs(generate_lefs(n_lefs)),
      _chromosomes(init_chromosomes_from_bed()),
      _seed(seed),
      _rand_gen(std::mt19937(_seed)) {}

std::vector<uint32_t> Genome::get_chromosome_lengths() const {
  std::vector<uint32_t> lengths;
  lengths.reserve(this->get_n_chromosomes());
  for (const auto& chr : this->_chromosomes) lengths.push_back(chr.length());
  return lengths;
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
        auto path_to_outfile = absl::StrFormat("%s/%s.tsv.bz2", output_dir, chr.name);
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
        path_to_outfile = absl::StrFormat("%s/%s_raw.tsv.bz2", output_dir, chr.name);
        absl::FPrintF(stderr, "Writing raw contact matrix for '%s' to file '%s'...\n", chr.name,
                      path_to_outfile);
        auto [bytes_in, bytes_out] = chr.contacts.write_to_tsv(path_to_outfile);
        absl::FPrintF(
            stderr, "DONE writing '%s' in %s! Compressed size: %.2f MB (compression ratio %.2fx)\n",
            path_to_outfile, absl::FormatDuration(absl::Now() - t0), bytes_out / 1.0e6,
            static_cast<double>(bytes_in) / bytes_out);
      });
}

void Genome::make_heatmaps(std::string_view output_dir, bool force_overwrite) const {
  const auto t0 = absl::Now();
  absl::FPrintF(stderr, "Generating heatmaps for %lu chromosomes...", this->get_n_chromosomes());
  std::string script = std::filesystem::is_regular_file("./utils/plot_contact_matrix.py")
                           ? "./utils/plot_contact_matrix.py"
                           : "plot_contact_matrix.py";
  auto cmd =
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

std::vector<Chromosome> Genome::init_chromosomes_from_bed() const {
  std::vector<Chromosome> chromosomes;
  SimpleBEDParser parser(this->_path_to_bed);
  SimpleBED record = parser.parse_next();

  do {
    chromosomes.emplace_back(record.chr, record.chr_end - record.chr_start, this->_bin_size,
                             this->_avg_lef_processivity);
    record = parser.parse_next();
  } while (!record.empty());

  return chromosomes;
}

std::vector<Lef> Genome::generate_lefs(uint32_t n) {
  std::vector<Lef> v;
  v.reserve(n);
  for (auto i = 0UL; i < n; ++i)
    v.emplace_back(this->_bin_size, this->_avg_lef_processivity,
                   this->_probability_of_extr_unit_bypass);
  return v;
}

void Genome::randomly_generate_barriers(uint32_t n_barriers) {
  const auto& weights = this->get_chromosome_lengths();
  std::discrete_distribution<> chr_idx(weights.begin(), weights.end());
  std::bernoulli_distribution strand_selector(0.5);

  for (auto i = 0UL; i < n_barriers; ++i) {
    // Randomly select a chromosome, barrier binding pos and direction
    auto& chr = this->_chromosomes[chr_idx(this->_rand_gen)];
    std::uniform_int_distribution<uint64_t> uniform_rng(0, chr.length());
    const auto barrier_position = uniform_rng(this->_rand_gen);
    const DNA::Direction direction =
        strand_selector(this->_rand_gen) ? DNA::Direction::rev : DNA::Direction::fwd;

    // Add the new extrusion barrier to the appropriate bin
    auto& bin = chr.dna.get_bin_from_pos(barrier_position);
    bin.add_extr_barrier(barrier_position, this->_probability_of_barrier_block, direction);
    chr.barriers.emplace_back(&bin.get_all_extr_barriers()->back());
  }
  for (auto& chr : this->get_chromosomes()) chr.sort_barriers_by_pos();
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

  for (auto& lef :
       this->_lefs) {  // Randomly select a chromosome, barrier binding pos and direction
    auto& chr = this->_chromosomes[chr_idx(this->_rand_gen)];
    std::uniform_int_distribution<uint32_t> uniform_rng(0, chr.length());
    auto pos = uniform_rng(this->_rand_gen);

    lef.bind_at_pos(chr, pos, this->_rand_gen);  // And bind to it
  }
}

void Genome::simulate_extrusion(uint32_t iterations) {
  const auto step = iterations / 500;
  const auto& weights = this->get_chromosome_lengths();
  std::discrete_distribution<> chr_idx(weights.begin(), weights.end());

  auto t0 = absl::Now();
  for (auto i = 1UL; i <= iterations; ++i) {
    for (auto& lef : this->_lefs) {
      if (lef.is_bound()) {  // Register contact and extrude if LEF is bound
        lef.register_contact();
        lef.extrude();
      }
    }

    for (auto& lef : this->_lefs) {
      // Check whether the last round of extrusions caused a collision with another LEF or an extr.
      // boundary, apply the stall and/or increase LEF lifetime where appropriate
      if (lef.is_bound()) {
        lef.check_constraints(this->_rand_gen);
      } else {
        // Randomly select a chromosome and try to bind one of the free LEFs to it
        auto& chr = this->_chromosomes[chr_idx(this->_rand_gen)];
        lef.try_rebind(chr, this->_rand_gen, this->_probability_of_lef_rebind);
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