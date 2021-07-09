#pragma once
#include <cuda_runtime_api.h>
#include <curand_kernel.h>

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <cuda/runtime_api.hpp>
#include <deque>
#include <mutex>
#include <utility>
#include <vector>

#include "modle/config.hpp"
#include "modle/cu/common.hpp"
#include "modle/setup.hpp"

namespace modle::cu {

__constant__ extern struct Config config;  // NOLINT

namespace CTCF {
enum State : uint_fast8_t { NOT_OCCUPIED = 0, OCCUPIED = 1 };
}

class Simulation;

struct Task {  // NOLINT
  __host__ __device__ Task() = default;
  uint32_t id{};

  uint32_t chrom_id{};
  bp_t chrom_start{};
  bp_t chrom_end{};

  uint32_t cell_id{};
  uint32_t n_target_epochs{};
  uint64_t n_target_contacts{};

  uint32_t nlefs{};
  uint32_t nbarriers{};
  uint64_t seed{};
};

struct BlockState {  // NOLINT
  __host__ __device__ BlockState() = default;
  __host__ __device__ BlockState(const BlockState& other) = delete;
  __host__ __device__ BlockState(BlockState&& other) = delete;
  __host__ __device__ ~BlockState() = default;

  __host__ __device__ BlockState& operator=(const BlockState& other) = delete;
  __host__ __device__ BlockState& operator=(BlockState&& other) = delete;

  CTCF::State* barrier_mask{nullptr};

  bp_t* rev_unit_pos{nullptr};
  bp_t* fwd_unit_pos{nullptr};
  uint32_t* lef_rev_unit_idx{nullptr};
  uint32_t* lef_fwd_unit_idx{nullptr};

  uint32_t* rev_moves_buff{nullptr};
  uint32_t* fwd_moves_buff{nullptr};
  collision_t* rev_collision_mask{nullptr};
  collision_t* fwd_collision_mask{nullptr};

  uint32_t num_rev_units_at_5prime{0};
  uint32_t num_fwd_units_at_3prime{0};

  float* lef_unloader_affinities{nullptr};
  float* lef_unloader_affinities_prefix_sum{nullptr};
  uint2* contact_local_buff{nullptr};
  uint32_t contact_local_buff_size{};
  uint32_t contact_local_buff_capacity{};

  uint32_t* loading_epochs{nullptr};
  uint32_t* lefs_to_load_per_epoch{nullptr};
  uint32_t num_unique_loading_epochs{0};
  uint32_t epoch_idx{0};

  curandStatePhilox4_32_10_t* rng_state{nullptr};

  uint32_t* tmp_lef_buff1{nullptr};
  uint32_t* tmp_lef_buff2{nullptr};
  uint32_t* tmp_lef_buff3{nullptr};
  uint32_t* tmp_lef_buff4{nullptr};

  void* cub_tmp_storage{nullptr};
  size_t cub_tmp_storage_bytes{0};

  uint32_t num_active_lefs{0};
  uint32_t num_lefs_to_release{0};
  bool burnin_completed{false};
  bool simulation_completed{false};
};

class GlobalStateHost;

struct GlobalStateDev {  // NOLINT
  __host__ __device__ GlobalStateDev() = default;
  __host__ __device__ GlobalStateDev(const GlobalStateDev& other) = default;
  __host__ __device__ GlobalStateDev(GlobalStateDev&& other) = default;

  Config* _config{nullptr};  // This field points to a region in __constant__ memory

  size_t grid_size{};
  size_t block_size{};
  uint32_t shared_memory_bytes_per_block{0};

  BlockState* block_states{nullptr};
  Task* tasks{nullptr};

  uint32_t current_epoch{0};
  uint32_t nblock_states{};
  uint32_t ntasks{};
  uint32_t ntasks_completed{};

  bp_t* barrier_pos{nullptr};
  dna::Direction* barrier_directions{nullptr};
  float* barrier_probs_occ_to_occ{nullptr};
  float* barrier_probs_nocc_to_nocc{nullptr};
};

class GlobalStateHost {  // NOLINT
  friend Simulation;

 public:
  GlobalStateHost() = default;
  GlobalStateHost(size_t grid_size_, size_t block_size_, size_t device_heap_size,
                  size_t max_grid_size_, const Config& c, size_t max_nlefs, size_t max_nbarriers,
                  size_t max_ncontacts_per_block,
                  cuda::device_t dev = cuda::device::current::get());
  GlobalStateHost(const GlobalStateHost& other) = delete;
  GlobalStateHost(GlobalStateHost&& other) = delete;
  ~GlobalStateHost();

  [[nodiscard]] size_t get_grid_size() const;
  [[nodiscard]] size_t get_block_size() const;

  [[nodiscard]] GlobalStateDev* get_ptr_to_dev_instance();
  [[nodiscard]] GlobalStateDev get_copy_of_device_instance();
  void write_tasks_to_device(const std::vector<Task>& new_tasks);
  void write_barriers_to_device(const std::vector<uint32_t>& new_barrier_positions,
                                const std::vector<dna::Direction>& new_barrier_directions,
                                const std::vector<float>& new_barrier_probs_occ_to_occ,
                                const std::vector<float>& new_barrier_probs_nocc_to_nocc);

  void sync_state_with_device();
  void sync_state_with_device(const GlobalStateDev& state);
  void init(size_t grid_size_, size_t block_size_, const std::vector<Task>& tasks_host = {},
            const std::vector<uint32_t>& barrier_pos_host = {},
            const std::vector<dna::Direction>& barrier_dir_host = {},
            const std::vector<float>& barrier_probs_occ_to_occ_host = {},
            const std::vector<float>& barrier_probs_nocc_to_nocc_host = {});

  Config* _config;  // This field points to a region in __constant__ memory

  size_t grid_size;
  size_t block_size;
  uint32_t shared_memory_bytes_per_block{0};
  size_t max_grid_size;

  BlockState* block_states{nullptr};
  Task* tasks{nullptr};

  cuda::memory::device::unique_ptr<bp_t[]> barrier_pos{nullptr};                   // NOLINT
  cuda::memory::device::unique_ptr<dna::Direction[]> barrier_directions{nullptr};  // NOLINT
  cuda::memory::device::unique_ptr<float[]> barrier_probs_occ_to_occ{nullptr};     // NOLINT
  cuda::memory::device::unique_ptr<float[]> barrier_probs_nocc_to_nocc{nullptr};   // NOLINT

  uint32_t current_epoch{0};
  uint32_t nblock_states{};
  uint32_t ntasks{};
  uint32_t ntasks_completed{};

 private:
  cuda::memory::region_t _global_state_dev{};
  cuda::device_t _device{cuda::device::current::get()};

  std::vector<uint32_t> _buff1;
  uint32_t _max_ncontacts_per_block;

  [[nodiscard]] static BlockState* allocate_block_states(const cuda::device_t& dev,
                                                         size_t max_nbarriers, size_t max_nlefs,
                                                         size_t max_ncontacts, size_t max_grid_size,
                                                         size_t max_block_size);
  static void deallocate_block_states(BlockState* block_states_dev, size_t nstates);
};

class Simulation : modle::Config {
 public:  // NOLINTNEXTLINE
  Simulation(const modle::Config& config_, size_t grid_size = 208'896 / 192,
             size_t block_size = 192, size_t max_grid_size = std::numeric_limits<size_t>::max(),
             size_t device_heap_size = 1024ULL * 1024ULL * 2048LL);
  ~Simulation() = default;
  Simulation(const Simulation& other) = delete;
  Simulation(const Simulation&& other) = delete;
  inline std::string to_string() = delete;
  inline void print() = delete;
  Simulation& operator=(const Simulation& other) = delete;
  Simulation& operator=(Simulation&& other) = delete;

  static constexpr auto EXTR_UNIT_IS_IDLE = static_cast<bp_t>(-1);
  static constexpr auto EXTR_UNIT_AT_CHROM_BOUNDARY = EXTR_UNIT_IS_IDLE - 1;
  static constexpr auto NO_COLLISIONS = static_cast<collision_t>(-1);

  void run();
  void sync();

 private:
  const modle::Config* _config;
  modle::Genome _genome{};

  size_t _grid_size;
  size_t _block_size;
  GlobalStateHost _global_state_host;
  std::vector<Task> _tasks{};
  std::vector<uint32_t> _barrier_positions{};
  std::vector<cu::dna::Direction> _barrier_directions{};
  std::vector<float> _barrier_probs_occ_to_occ{};
  std::vector<float> _barrier_probs_nocc_to_nocc{};

  size_t _current_epoch{0};

  [[nodiscard]] size_t compute_shared_memory_per_block() const;

  void run_batch(const std::vector<Task>& new_tasks, const std::vector<uint32_t>& barrier_pos,
                 const std::vector<dna::Direction>& barrier_dir,
                 const std::vector<float>& barrier_probs_occ_to_occ,
                 const std::vector<float>& barrier_probs_nocc_to_nocc, size_t batch_number);

  void write_contacts_to_disk(std::deque<std::pair<Chromosome*, size_t>>& progress_queue,
                              std::mutex& progress_queue_mutex,
                              std::atomic<bool>& end_of_simulation);
  void update_contacts_for_chrom(Chromosome& chrom);
};

[[nodiscard]] Config* write_config_to_device(const Config& c);

void throw_on_cuda_error(const cudaError status, std::string_view exception_prefix_message = "");

}  // namespace modle::cu