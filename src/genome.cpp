#include "modle/genome.hpp"

#include <absl/strings/str_format.h>

#include <functional>

namespace modle {

Genome::Genome(std::string_view path_to_bed, uint32_t bin_size, uint32_t n_lefs,
               uint32_t n_barriers, uint32_t avg_lef_processivity,
               double probability_of_barrier_block, uint64_t seed)
    : _path_to_bed(path_to_bed),
      _bin_size(bin_size),
      _chromosomes(init_chromosomes_from_bed()),
      _rng_dev(std::default_random_engine(seed)),
      _lefs(generate_lefs(n_lefs, avg_lef_processivity)),
      _avg_lef_processivity(avg_lef_processivity),
      _barriers(generate_barriers(n_barriers, probability_of_barrier_block)),
      _probability_of_barrier_block(probability_of_barrier_block),
      _seed(seed) {}

absl::flat_hash_map<const std::string, DNA> Genome::init_chromosomes_from_bed() {
  absl::flat_hash_map<const std::string, DNA> chromosomes;
  SimpleBEDParser parser(this->_path_to_bed);
  SimpleBED record = parser.parse_next();
  do {
    chromosomes.emplace(record.chr, DNA{record.chr_end - record.chr_start, this->_bin_size});
    record = parser.parse_next();
  } while (!record.empty());

  return chromosomes;
}

std::vector<Lef> Genome::generate_lefs(uint32_t n, uint32_t avg_processivity) {
  return std::vector<Lef>{n, Lef{0, avg_processivity}};
}

absl::btree_multimap<std::string_view, ExtrusionBarrier> Genome::generate_barriers(
    uint32_t n, double prob_of_block) {
  std::hash<std::string_view> str_hasher;
  absl::btree_multimap<std::string_view, ExtrusionBarrier> barriers;
  const uint32_t avg_tad_size = this->size() / n;
  for (const auto& [chr, dna] : this->_chromosomes) {
    const uint32_t n_barriers = dna.length() / avg_tad_size;
    //    absl::FPrintF(stderr, "n_barriers = %lu\n", n_barriers);
    std::uniform_int_distribution<uint32_t> rng(0, dna.length());
    absl::flat_hash_set<uint32_t> barrier_positions;
    while (barrier_positions.size() < n_barriers) {
      barrier_positions.insert(rng(this->_rng_dev));
    }
    for (const auto& pos : barrier_positions) {
      // We are using the hash of the chr name + the position as seed for the rng that controls
      // whether a boundary will block or not
      barriers.emplace(
          chr, ExtrusionBarrier{pos, this->_bin_size, prob_of_block, str_hasher(chr) + pos});
    }
  }
  return barriers;
}

const absl::flat_hash_map<const std::string, DNA>& Genome::get_chromosomes() const {
  return this->_chromosomes;
}

absl::flat_hash_map<const std::string, DNA>& Genome::get_chromosomes() {
  return this->_chromosomes;
}

std::vector<Lef>& Genome::get_lefs() { return this->_lefs; }
const std::vector<Lef>& Genome::get_lefs() const { return this->_lefs; }

absl::btree_multimap<std::string_view, ExtrusionBarrier>& Genome::get_barriers() {
  return this->_barriers;
}
const absl::btree_multimap<std::string_view, ExtrusionBarrier>& Genome::get_barriers() const {
  return this->_barriers;
}

uint32_t Genome::n_chromosomes() const { return this->_chromosomes.size(); }

uint64_t Genome::size(uint32_t min_length) const {
  return std::accumulate(this->_chromosomes.begin(), this->_chromosomes.end(), 0ULL,
                         [&](uint64_t accumulator, const std::pair<const std::string, DNA>& chr) {
                           return accumulator +
                                  chr.second.length() * (chr.second.length() > min_length);
                         });
}

uint32_t Genome::n_bins() const {
  return std::accumulate(this->_chromosomes.begin(), this->_chromosomes.end(), 0UL,
                         [](uint64_t accumulator, const std::pair<const std::string, DNA>& chr) {
                           return accumulator + chr.second.nbins();
                         });
}

uint32_t Genome::tot_n_lefs() const { return this->_lefs.size(); }

uint32_t Genome::n_barriers() const { return this->_barriers.size(); }

void Genome::randomly_bind_lefs(uint32_t min_length) {
  const auto genome_size = this->size(min_length);
  auto lef_it = this->_lefs.begin();
  std::vector<DNA::Bin*> sampled_bins;

  for (auto& [chr, dna] : this->_chromosomes) {
    if (dna.length() < min_length) continue;
    const uint32_t n_lefs = (static_cast<double>(dna.length()) / genome_size) * this->tot_n_lefs();
    {  // This ugly trick is required to avoid the increase in ref. count for shared_ptr<Bin>, as
       // this will cause bind_at_pos to return false (and bind nothing)
      // A cleaner solution would be to somehow customize back_inserted to push_back a raw_ptr
      // instead of a smart_ptr
      std::vector<std::shared_ptr<DNA::Bin>> tmp_vect;
      tmp_vect.reserve(n_lefs);
      std::sample(dna.begin(), dna.end(), std::back_inserter(tmp_vect), n_lefs, this->_rng_dev);
      sampled_bins.clear();
      sampled_bins.reserve(n_lefs);
      for (const auto& ptr : tmp_vect) {
        sampled_bins.push_back(ptr.get());
      }
    }
    for (const auto& bin : sampled_bins) {
      if (!lef_it++->bind_at_pos(chr, dna, bin->get_center())) {
        absl::FPrintF(
            stderr, "%s:%lu is already occupied (bin.use_count()=%lu)\n", chr, bin->get_center(),
            this->_chromosomes.at(chr).get_ptr_to_bin_from_pos(bin->get_center()).use_count() - 1);
        --lef_it;
      }
    }
  }
  //  absl::FPrintF(stderr, "diff=%lu\n", this->_lefs.end() - lef_it);
  //  assert(lef_it == this->_lefs.end());
}

void Genome::uniformly_bind_lefs(uint32_t min_length) {
  const auto genome_size = this->size(min_length);
  const auto avg_lef_distance = genome_size / this->tot_n_lefs();
  auto lef_it = this->_lefs.begin();

  for (auto& [chr, dna] : this->_chromosomes) {
    if (dna.length() < min_length) continue;
    uint32_t pos = 0;
    do {
      absl::FPrintF(stderr, "chr=%s; pos=%lu; avg_lef_dist=%lu; dna.length()=%lu;\n", chr, pos,
                    avg_lef_distance, dna.length());
      if (!lef_it++->bind_at_pos(chr, dna, pos)) {
        absl::FPrintF(stderr,
                      "This shouldn't happen... %s:%lu is already occupied (bin.use_count()=%lu)\n",
                      chr, pos, dna.get_ptr_to_bin_from_pos(pos).use_count() - 1);
        --lef_it;
      }
      pos += avg_lef_distance;
    } while (pos < dna.length() && lef_it + 1 != this->_lefs.end());
  }
}

}  // namespace modle