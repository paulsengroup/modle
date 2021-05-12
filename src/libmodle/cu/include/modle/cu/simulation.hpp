#pragma once
#include <cuda_runtime_api.h>
#include <curand_kernel.h>

#include <cstddef>
#include <cstdint>
#include <cuda/runtime_api.hpp>
#include <vector>

#include "modle/config.hpp"
#include "modle/cu/common.hpp"
#include "modle/libmodle_io.hpp"

namespace modle::cu {
// using curandStatePhilox4_32_10_t = struct curandStatePhilox4_32_10;
class GlobalStateHost;

class Simulation : modle::Config {
 public:  // NOLINTNEXTLINE
  Simulation(const modle::Config& config, size_t grid_size = 208'896 / 192, size_t block_size = 192,
             size_t max_grid_size = std::numeric_limits<size_t>::max(),
             size_t device_heap_size = 1024ULL * 1024ULL * 2048LL);
  Simulation(const Simulation& other) = delete;
  inline std::string to_string() = delete;
  inline void print() = delete;
  Simulation& operator=(const Simulation& other) = delete;
  /*
    void init_global_state(size_t grid_size, size_t block_size, size_t max_grid_size,
                           size_t max_nlefs, size_t max_nbarriers,
                           cuda::device_t dev = cuda::device::current::get());
                           */
  static constexpr auto LEF_IS_IDLE = static_cast<uint32_t>(-1);

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

  void run();

 private:
  const modle::Config* _config;
  modle::io::Genome _genome{};

  size_t _grid_size;
  size_t _block_size;
  std::unique_ptr<GlobalStateHost> _global_state_host;
  std::vector<Task> _tasks{};
  std::vector<uint32_t> _barrier_positions{};
  std::vector<cu::dna::Direction> _barrier_directions{};

  size_t _current_epoch{0};

  void run_batch(const std::vector<Task>& new_tasks, const std::vector<uint32_t>& barrier_pos,
                 std::vector<dna::Direction>& barrier_dir);
  void update_tasks(const std::vector<Task>& new_tasks) noexcept;
  void update_barriers(const std::vector<uint32_t>& barrier_pos,
                       const std::vector<dna::Direction>& barrier_dir) noexcept;

  bool simulate_one_epoch();
  void setup_burnin_phase();
  void sort_lefs();
};

namespace CTCF {
enum State : uint_fast8_t { NOT_OCCUPIED = 0, OCCUPIED = 1 };
}

struct BlockState {  // NOLINT
  __host__ __device__ BlockState() = default;
  __device__ BlockState(uint32_t block_size, uint32_t nbarriers, uint32_t nlefs,
                        uint32_t contact_buff_size = 0);
  __host__ __device__ BlockState(const BlockState& other) = delete;
  __host__ __device__ BlockState(BlockState&& other) = delete;
  __host__ __device__ ~BlockState() = default;

  __host__ __device__ BlockState& operator=(const BlockState& other) = delete;
  __host__ __device__ BlockState& operator=(BlockState&& other) = delete;

  bool* barrier_mask{nullptr};

  bp_t* rev_unit_pos{nullptr};
  bp_t* fwd_unit_pos{nullptr};
  uint32_t* lef_rev_unit_idx{nullptr};
  uint32_t* lef_fwd_unit_idx{nullptr};

  uint32_t* rev_moves_buff{nullptr};
  uint32_t* fwd_moves_buff{nullptr};
  collision_t* rev_collision_mask{nullptr};
  collision_t* fwd_collision_mask{nullptr};

  float* lef_unloader_affinities{nullptr};
  uint2* contact_local_buff{nullptr};
  uint32_t contact_local_buff_size{};
  uint32_t contact_local_buff_capacity{};
  uint32_t* epoch_buff{nullptr};

  curandStatePhilox4_32_10_t* rng_state{nullptr};

  uint32_t num_active_lefs{0};
  bool burnin_completed{false};
};

class GlobalStateHost;

struct GlobalStateDev {  // NOLINT
  __host__ __device__ GlobalStateDev() = default;
  __host__ __device__ GlobalStateDev(const GlobalStateDev& other) = default;
  __host__ __device__ GlobalStateDev(GlobalStateDev&& other) = default;

  Config* config{nullptr};  // This field points to a region in __constant__ memory

  size_t grid_size{};
  size_t block_size{};

  BlockState* block_states{nullptr};
  Simulation::Task* tasks{nullptr};

  uint32_t nblock_states{};
  uint32_t ntasks{};

  bp_t* barrier_pos{nullptr};
  dna::Direction* barrier_directions{nullptr};
  // float* barrier

  uint32_t* large_uint_buff1{nullptr};
  uint32_t* large_uint_buff2{nullptr};
  uint32_t* large_uint_buff3{nullptr};
  uint32_t* large_uint_buff4{nullptr};

  uint32_t large_uint_buff_chunk_alignment{};

  uint32_t* sorting_offset1_buff{nullptr};
  uint32_t* sorting_offset2_buff{nullptr};

  void* tmp_sorting_storage{nullptr};
  size_t tmp_sorting_storage_bytes{0};
};

class GlobalStateHost {  // NOLINT
  friend Simulation;

 public:
  GlobalStateHost() = default;
  GlobalStateHost(size_t grid_size_, size_t block_size_, size_t max_grid_size_, const Config& c,
                  size_t max_nlefs, size_t max_nbarriers,
                  cuda::device_t dev = cuda::device::current::get());
  GlobalStateHost(const GlobalStateHost& other) = delete;
  GlobalStateHost(GlobalStateHost&& other) = default;
  ~GlobalStateHost();

  [[nodiscard]] size_t get_grid_size() const;
  [[nodiscard]] size_t get_block_size() const;

  [[nodiscard]] GlobalStateDev* get_ptr_to_dev_instance();
  [[nodiscard]] GlobalStateDev get_copy_of_device_instance();
  void write_tasks_to_device(const std::vector<Simulation::Task>& new_tasks);
  void write_barriers_to_device(const std::vector<uint32_t>& new_barrier_positions,
                                const std::vector<dna::Direction>& new_barrier_directions);

  void sync_state_with_device();
  void sync_state_with_device(const GlobalStateDev& state);
  void init(size_t grid_size_, size_t block_size_,
            const std::vector<Simulation::Task>& tasks_host = {},
            const std::vector<uint32_t>& barrier_pos_host = {},
            const std::vector<dna::Direction>& barrier_dir_host = {});

  Config* config;  // This field points to a region in __constant__ memory

  size_t grid_size;
  size_t block_size;
  size_t max_grid_size;

 private:
  cuda::memory::region_t _block_states_buff;
  cuda::memory::region_t _tasks_buff;

 public:
  BlockState* block_states{nullptr};
  Simulation::Task* tasks{nullptr};

  cuda::memory::device::unique_ptr<bp_t[]> barrier_pos{nullptr};                   // NOLINT
  cuda::memory::device::unique_ptr<dna::Direction[]> barrier_directions{nullptr};  // NOLINT
  // float* barrier

  uint32_t nblock_states{};
  uint32_t ntasks{};

  cuda::memory::device::unique_ptr<uint32_t[]> large_uint_buff1{nullptr};      // NOLINT
  cuda::memory::device::unique_ptr<uint32_t[]> large_uint_buff2{nullptr};      // NOLINT
  cuda::memory::device::unique_ptr<uint32_t[]> large_uint_buff3{nullptr};      // NOLINT
  cuda::memory::device::unique_ptr<uint32_t[]> large_uint_buff4{nullptr};      // NOLINT
  cuda::memory::device::unique_ptr<uint32_t[]> sorting_offset1_buff{nullptr};  // NOLINT
  cuda::memory::device::unique_ptr<uint32_t[]> sorting_offset2_buff{nullptr};  // NOLINT
  cuda::memory::device::unique_ptr<uint32_t[]> tmp_sorting_storage{nullptr};   // NOLINT
  size_t tmp_sorting_storage_bytes = 0;

  uint32_t large_uint_buff_chunk_alignment{};

 private:
  cuda::memory::region_t _global_state_dev{};
  cuda::device_t _device{cuda::device::current::get()};

  std::vector<uint32_t> _buff1;
};

[[nodiscard]] Config* write_config_to_device(const Config& c);

}  // namespace modle::cu
