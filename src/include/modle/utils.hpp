#pragma once

#include <algorithm>
#include <bitset>
#include <boost/dynamic_bitset.hpp>
#include <cstdint>
#include <execution>
#include <iostream>
#include <numeric>
#include <type_traits>
#include <vector>

namespace modle {
template <typename N, typename = typename std::enable_if<std::is_arithmetic<N>::value, N>::type>
void print_vect(const std::vector<N> &v) {
  for (const auto &n : v) {
    std::cout << n << "\t";
  }
  std::cout << std::endl;
}
/*
// Return the indices corresponding to the sorted vector
template <typename N, typename = typename std::enable_if<std::is_arithmetic<N>::value, N>::type>
std::vector<uint32_t> sort_vector_by_idx(const std::vector<N> &v) {
  std::vector<uint32_t> vi(v.size());
  std::iota(vi.begin(), vi.end(), 0);
  std::sort(std::execution::par_unseq, vi.begin(), vi.end(),
            [&](uint32_t i1, uint32_t i2) { return v[i1] < v[i2]; });
  return vi;
}
 */

// Return the indices corresponding to the sorted vector
template <typename Iterator,
          typename = typename std::enable_if<
              std::is_base_of<std::random_access_iterator_tag,
                              typename std::iterator_traits<Iterator>::iterator_category>::value>>
std::vector<uint32_t> sort_vector_by_idx(Iterator begin, Iterator end) {
  std::vector<uint32_t> vi(std::distance(begin, end));
  std::iota(vi.begin(), vi.end(), 0);
  std::sort(vi.begin(), vi.end(),
            [&](uint32_t i1, uint32_t i2) { return *(begin + i1) < *(begin + i2); });
  return vi;
}
// Given to vectors A and B, sort pairs of values based on vector A, and solve ties with vector B
template <typename Iterator,
          typename = typename std::enable_if<
              std::is_base_of<std::random_access_iterator_tag,
                              typename std::iterator_traits<Iterator>::iterator_category>::value>>
void sort_pair_of_vectors(Iterator v1_begin, Iterator v1_end, Iterator v2_begin, Iterator v2_end) {
  assert(std::distance(v1_begin, v1_end) == std::distance(v2_begin, v2_end));
  const auto vi = sort_vector_by_idx(v1_begin, v2_end);
  for (auto i = 0UL; i < vi.size(); ++i) {
    if (auto j = vi[i]; j != i && j > i) {
      std::swap(v1_begin + i, v1_begin + j);
      std::swap(v2_begin + i, v2_end + j);
    }
  }
  auto start = v2_begin;
  auto end = start;
  auto size = std::distance(v1_begin, v1_end);
  for (auto i = 1UL; i < size; ++i) {
    if (*(v1_begin + i) == *(v1_begin + i - 1)) {
      start = v2_begin + i - 1;
      for (++i; i < size; ++i) {
        if (*start != *(v1_begin + i)) {
          // end is 1 place after the last element we want to sort,
          // so we don't subtract 1
          end = v2_begin + i;
          std::sort(start, end);
          break;
        }
      }
    }
  }
}

// This returns a vector of integers corresponding to the rank of the values of
// vector v
template <typename Iterator,
          typename = typename std::enable_if<
              std::is_base_of<std::random_access_iterator_tag,
                              typename std::iterator_traits<Iterator>::iterator_category>::value>>
std::vector<double> compute_element_ranks(Iterator begin, Iterator end) {
  const auto vi = sort_vector_by_idx(begin, end);
  std::vector<double> vr(vi.size());

  vr[vi[0]] = 0;
  for (auto i = 1U; i < vr.size(); ++i) {
    auto prev = *(begin + vi[i - 1]);
    auto current = *(begin + vi[i]);
    if (prev != current) {
      vr[vi[i]] = i;
    } else {
      const auto min = i - 1;
      auto max = i;
      for (auto j = max + 1; current == *(begin + vi[j]) && j < vr.size(); ++j) max = j;
      for (auto j = min; j <= max; ++j) vr[vi[j]] = (max + min) / 2.0;
      i = max;
    }
  }
  //  print_vect(vr);
  return vr;
}

/*
template <typename Iterator,
    typename = typename std::enable_if<
        std::is_base_of<std::random_access_iterator_tag,
            typename std::iterator_traits<Iterator>::iterator_category>::value>>
std::vector<double> compute_element_ranks_w_ties(Iterator begin, Iterator end) {
  const auto vi = sort_vector_by_idx(begin, end);
  std::vector<double> vr(vi.size());
  absl::flat_hash_map<int64_t, std::pair<uint32_t, uint32_t>>
      ties;  // key = tied_value; value = <lowest_rank, highest_rank>
  ties.reserve(vi.size());
  // This bitset is just used to avoid expensive lookups in the hashmap
  // Bits corresponding to the indices of tied values are set to 1;
  boost::dynamic_bitset<> ties_mask(vr.size());

  for (auto i = 0U; i < vi.size(); ++i) {
    auto n = *(begin + vi[i]);                     // Number whose rank is being computed
    if (auto m = ties.find(n); m != ties.end()) {  // n is tied. Update min/max rank
      m->second.first = std::min(m->second.first, i);
      m->second.second = std::max(m->second.second, i);
      ties_mask[m->second.first] = true;
      ties_mask[m->second.second] = true;
    } else  // First time that we see n. Set min and max ramks to i
      ties.emplace(n, std::make_pair(i, i));
  }

  for (auto i = 0U; i < vr.size(); ++i) {
    if (!ties_mask[i])
      vr[vi[i]] = i;  // Value is not tied: assign normal rank
    else {
      const auto &[first, last] = ties.at(*(begin + vi[i]));
      vr[vi[i]] = (first + last) / 2.0;  // Value is tied: assign the corrected rank
    }
  }
  //  print_vect(vr);
  return vr;
}
 */

}  // namespace modle