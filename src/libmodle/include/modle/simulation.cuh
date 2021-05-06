#pragma once
#include <cstdint>
#include <cuda/std/atomic>
#include <vector>

#include "modle/common.hpp"
#include "modle/contacts.cuh"
#include "modle/extrusion_barriers.cuh"

namespace modle::cu::Simulation {
__global__ void mykernel();
__global__ void mykernel2(modle::cu::ContactMatrix<uint32_t>* m, cuda::std::atomic<uint32_t>* buff,
                          size_t nrows, size_t ncols, size_t& missed_updates, size_t& tot_contacts);

__global__ void mykernel3(const ExtrusionBarrier* barriers, size_t nbarriers, CTCF::State* mask,
                          curandStatePhilox4_32_10_t* rng_state);
[[nodiscard]] std::vector<uint32_t> run_mykernel2(size_t nrows, size_t ncols,
                                                  size_t& missed_updates, size_t& tot_contacts);

[[nodiscard]] std::vector<CTCF::State> run_mykernel3(
    const std::vector<ExtrusionBarrier>& host_barriers, uint64_t seed = 123456789);
/*
struct State {
  inline State() = default;
  size_t id{};
  char* chrom_name{nullptr};
  size_t chrom_start{};
  size_t chrom_end{};
  size_t chrom_size{};
  size_t cell_id{};
  size_t n_target_epochs{};
  size_t n_target_contacts{};
  size_t nlefs{};
  size_t nbarriers{};
  size_t* barrier_pos{};
  size_t* barrier_directions{};
  std::vector<Lef> lef_buff{};
  std::vector<double> lef_unloader_affinity{};
  std::vector<size_t> rank_buff1{};
  std::vector<size_t> rank_buff2{};
  boost::dynamic_bitset<> barrier_mask{};
  std::vector<bp_t> moves_buff1{};
  std::vector<bp_t> moves_buff2{};
  std::vector<size_t> idx_buff1{};
  std::vector<size_t> idx_buff2{};
  std::vector<size_t> epoch_buff{};
  modle::PRNG_t rand_eng{};
  uint64_t seed{};

  inline void resize(size_t size = std::numeric_limits<size_t>::max());
  inline void reset();
  [[nodiscard]] inline std::string to_string() const noexcept;
  // NOLINTNEXTLINE(readability-magic-numbers, cppcoreguidelines-avoid-magic-numbers)
};
*/
}  // namespace modle::cu::Simulation
