#pragma once

#include <cuda.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

#include <cuda/std/atomic>       // for atomic
#include <cuda/std/cstddef>      // IWYU pragma: keep for size_t
#include <cuda/std/cstdint>      // for uint64_t
#include <cuda/std/type_traits>  // for is_integral
#include <cuda/std/utility>      // for pair

namespace modle::cu {

template <typename I>
class ContactMatrix {
  static_assert(std::is_integral_v<I>,
                "ContactMatrix requires an integral type as template argument.");

 public:
  // Constructors
  __device__ inline ContactMatrix() = default;
  __device__ inline ContactMatrix(size_t nrows, size_t ncols);
  __device__ inline ~ContactMatrix();
  __device__ inline ContactMatrix(const ContactMatrix& other) = delete;
  __device__ inline ContactMatrix(ContactMatrix&& other) = delete;
  __device__ inline ContactMatrix& operator=(const ContactMatrix& other) = delete;
  __device__ inline ContactMatrix& operator=(ContactMatrix&& other) = delete;

  __device__ inline void print() const;

  // Counts getters and setters
  [[nodiscard]] __device__ inline I get(size_t row, size_t col);
  __device__ inline void set(size_t row, size_t col, I n);
  __device__ inline void add(size_t row, size_t col, I n);
  __device__ inline void increment(size_t row, size_t col);

  // Shape/statistics getters

  [[nodiscard]] __device__ inline constexpr size_t ncols() const;
  [[nodiscard]] __device__ inline constexpr size_t nrows() const;
  [[nodiscard]] __device__ inline constexpr size_t npixels() const;
  [[nodiscard]] __device__ inline constexpr size_t get_n_of_missed_updates() const;
  [[nodiscard]] __device__ inline constexpr size_t get_tot_contacts() const;

  // Misc
  __device__ inline void reset(cuda::std::atomic<I>* buff = nullptr, size_t nrows = 0,
                               size_t ncols = 0);
  [[nodiscard]] __device__ inline bool empty() const;

 private:
  size_t _nrows{0};
  size_t _ncols{0};
  cuda::std::atomic<I>* _contacts{nullptr};
  cuda::std::atomic<size_t> _tot_contacts{0};
  cuda::std::atomic<size_t> _updates_missed{0};

  [[nodiscard]] __device__ static inline cuda::std::pair<size_t, size_t> transpose_coords(
      size_t row, size_t col);
  [[nodiscard]] __device__ inline size_t encode_coords(size_t row, size_t col) const;
  [[nodiscard]] __device__ inline cuda::std::pair<size_t, size_t> decode_coords(size_t n) const;
};

template <typename I>
__device__ ContactMatrix<I>::ContactMatrix(size_t nrows, size_t ncols)
    : _nrows(nrows), _ncols(ncols), _contacts(new I[_nrows * _ncols]) {}

template <typename I>
__device__ ContactMatrix<I>::~ContactMatrix() {
  delete[] this->_contacts;
}

template <typename I>
__device__ void ContactMatrix<I>::print() const {
  const auto id = threadIdx.x + blockIdx.x * blockDim.x;
  if (id == 0) {
    for (auto i = 0UL; i < this->nrows() * this->ncols(); ++i) {
      printf("i=%lu -> %u\n", i, this->_contacts[i]);
    }
  }
}

template <typename I>
__device__ constexpr size_t ContactMatrix<I>::nrows() const {
  return this->_nrows;
}

template <typename I>
__device__ constexpr size_t ContactMatrix<I>::ncols() const {
  return this->_ncols;
}

template <typename I>
__device__ constexpr size_t ContactMatrix<I>::npixels() const {
  return this->_nrows * this->_ncols;
}

template <typename I>
__device__ constexpr size_t ContactMatrix<I>::get_n_of_missed_updates() const {
  return this->_updates_missed;
}

template <typename I>
__device__ constexpr size_t ContactMatrix<I>::get_tot_contacts() const {
  return this->_tot_contacts;
}

template <typename I>
__device__ I ContactMatrix<I>::get(size_t row, size_t col) {
  const auto [i, j] = transpose_coords(row, col);
  assert(j < this->_ncols);
  if (i > this->nrows()) {
    return 0;
  }
  const auto idx = this->encode_coords(i, j);
  return this->_contacts[idx];
}

template <typename I>
__device__ void ContactMatrix<I>::set(size_t row, size_t col, I n) {
  const auto [i, j] = transpose_coords(row, col);
  assert(j < this->_ncols);
  if (i > this->nrows()) {
    ++this->_updates_missed;
    return;
  }
  const auto idx = this->encode_coords(i, j);
  this->_tot_contacts -= this->_contacts[idx];
  this->_contacts[idx] = n;
  this->_tot_contacts += n;
}

template <typename I>
__device__ void ContactMatrix<I>::add(size_t row, size_t col, I n) {
  const auto [i, j] = transpose_coords(row, col);
  assert(j < this->_ncols);
  if (i > this->nrows()) {
    ++this->_updates_missed;
    return;
  }
  const auto idx = this->encode_coords(i, j);
  this->_contacts[idx] += n;
  this->_tot_contacts += n;
}

template <typename I>
__device__ void ContactMatrix<I>::increment(size_t row, size_t col) {
  this->add(row, col, 1);
}

template <typename I>
[[nodiscard]] __device__ cuda::std::pair<size_t, size_t> ContactMatrix<I>::transpose_coords(
    size_t row, size_t col) {
  if (row > col) {
    return cuda::std::make_pair(row - col, row);
  }
  return cuda::std::make_pair(col - row, col);
}

template <typename I>
[[nodiscard]] __device__ size_t ContactMatrix<I>::encode_coords(size_t row, size_t col) const {
  const auto [i, j] = transpose_coords(row, col);
  return (j * this->_nrows) + i;
}

template <typename I>
[[nodiscard]] __device__ cuda::std::pair<size_t, size_t> ContactMatrix<I>::decode_coords(
    size_t n) const {
  const auto row = n / this->_nrows;
  const auto col = n % this->_ncols;
  assert(col < this->_ncols);
  return cuda::std::make_pair(row, col);
}

template <typename I>
__device__ bool ContactMatrix<I>::empty() const {
  for (auto i = 0UL; i < this->npixels(); ++i) {
    if (this->_contacts[i] != 0) {
      return false;
    }
  }
  return true;
}

template <typename I>
__device__ inline void ContactMatrix<I>::reset(cuda::std::atomic<I>* buff, size_t nrows,
                                               size_t ncols) {
  const auto id = threadIdx.x + blockIdx.x * blockDim.x;
  if (id == 0) {
    if (buff) {
      this->_contacts = buff;
      this->_nrows = nrows;
      this->_ncols = ncols;
      this->_updates_missed = 0UL;
    } else {
      assert(!buff);       // NOLINT
      assert(nrows == 0);  // NOLINT
      assert(ncols == 0);  // NOLINT
      memset(this->_contacts, static_cast<I>(0), this->_nrows * this->_ncols * sizeof(I));
    }
  }
  //__syncthreads();
}

}  // namespace modle::cu
