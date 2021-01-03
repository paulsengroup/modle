#pragma once

#include <cstdint>      // for uint*_t
#include <memory>       // for _Destroy, allocator
#include <string>       // for basic_string, string
#include <string_view>  // for string_view
#include <utility>      // for pair
#include <vector>       // for vector

#include "modle/dna.hpp"   // for Chromosome
#include "modle/lefs.hpp"  // for Lef

namespace modle {

struct config;

/** \brief The Genome class acts as the entrypoint to simulate DNA loop extrusion.
 *
 * This class is used to store and access most of the simulation data and operations, such as the
 * chromosomes that make up the Genome, the loop extrusion factors (Lef) and the extrusion barriers
 * (ExtrusionBarrier).
 *
 * The simulation can be executed by calling Genome::simulate_extrusion with the number of
 * iterations to simulate.
 */
class Genome {
 public:
  inline explicit Genome(const config& c);

  // Getters
  [[nodiscard]] inline const std::vector<Chromosome>& get_chromosomes() const;
  [[nodiscard]] inline std::vector<Chromosome>& get_chromosomes();
  [[nodiscard]] inline std::vector<Lef>& get_lefs();
  [[nodiscard]] inline const std::vector<Lef>& get_lefs() const;
  [[nodiscard]] inline uint32_t get_n_chromosomes() const;
  [[nodiscard]] inline uint64_t size() const;
  [[nodiscard]] inline uint64_t n50() const;
  [[nodiscard]] inline uint64_t get_n_bins() const;
  [[nodiscard]] inline uint64_t get_n_lefs() const;
  [[nodiscard]] inline uint64_t get_n_barriers() const;
  [[nodiscard]] inline std::vector<std::string_view> get_chromosome_names() const;
  [[nodiscard]] inline std::vector<uint64_t> get_chromosome_lengths() const;
  [[nodiscard]] inline std::vector<double> get_chromosome_lef_affinities() const;
  [[nodiscard]] inline uint64_t get_n_of_free_lefs() const;
  [[nodiscard]] inline uint64_t get_n_of_busy_lefs() const;
  inline void write_contacts_to_file(std::string_view output_file);
  inline void write_extrusion_barriers_to_file(std::string_view output_dir,
                                               bool force_overwrite) const;

  /** Randomly generate and bind \p n_barriers ExtrusionBarrier%s
   *
   * This function randomly selects Chromosome%s using \p std::discrete_distribution and
   * Chromosome%s lengths as weight.
   *
   * After selecting a chromosome the binding position for an ExtrusionBarrier is randomly drawn
   * from an uniform distribution [0, Chromosome::simulated_length).
   *
   * The direction where the ExtrusionBarrier exerts a strong block is randomly drawn using an
   * unbiased bernoulli trial.
   *
   * ExtrusionBarrier%s are owned by the DNA::Bin to which they bound, while Chromosome::barriers
   * only stores ptrs to the actual ExtrusionBarrier%s.
   * @param n_barriers The number of barriers to randomly generate and bind.
   */
  template <typename I>
  inline void randomly_generate_extrusion_barriers(I n_barriers, uint64_t seed = 0);

  /** This function is meant to be used before calling Genome::simulate_extrusion to randomly bind
   * Lef%s.
   *
   * The Chromosome and binding position are selected following the same procedure described for
   * Chromosome::randomly_generate_extrusion_barriers.
   */

  inline std::pair<uint64_t, uint64_t> import_extrusion_barriers_from_bed(
      std::string_view path_to_bed, double probability_of_block);
  inline void assign_lefs(bool bind_lefs_after_assignment);

  inline uint64_t remove_chromosomes_wo_extr_barriers();

  inline std::pair<double, double> run_burnin(double prob_of_rebinding,
                                              uint32_t target_n_of_unload_events,
                                              uint64_t min_extr_rounds);

  /** This function is the one that actually runs the simulation
   *
   * Pseudocode for one simulation round:
   \verbatim
     For each lef in Genome::lefs:
        If lef is bound to DNA:
           Register contact
           Extrude
        If lef is bound to DNA:
           Check constraints and apply stalls where appropriate
        else:
           Try to rebind the lef
   \endverbatim
   *
   * @param iterations Number of loop extrusion events to simulate
   */
  void simulate_extrusion();
  void simulate_extrusion(uint32_t iterations);
  void simulate_extrusion(double target_contact_density);

 private:
  uint64_t _seed;
  // std::mt19937 _rand_eng;
  std::string _path_to_chrom_sizes_file;
  std::string _path_to_chr_subranges_file;
  uint32_t _bin_size;
  uint32_t _avg_lef_processivity;  ///< Average loop size if loop extrusion takes place unobstructed
  double _probability_of_barrier_block;
  double _probability_of_lef_rebind;
  double _probability_of_extr_unit_bypass;
  double _lef_unloader_strength_coeff;
  std::vector<Lef> _lefs;
  std::vector<Chromosome> _chromosomes;
  uint32_t _sampling_interval;
  bool _randomize_contact_sampling;
  uint32_t _nthreads;

  inline void simulate_extrusion(uint32_t iterations, double target_contact_density);

  // Private initializers
  [[nodiscard]] inline std::vector<Chromosome> init_chromosomes_from_file(
      uint32_t diagonal_width) const;
  [[nodiscard]] inline std::vector<Lef> generate_lefs(uint32_t n);
};
}  // namespace modle

#include "../../genome_impl.hpp"