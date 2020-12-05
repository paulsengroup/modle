#pragma once
#include <fmt/format.h>

#include <limits>
#include <stdexcept>
#include <string_view>
#include <type_traits>

#include "modle/suppress_compiler_warnings.hpp"

namespace modle {

DISABLE_WARNING_PUSH
DISABLE_WARNING_CONVERSION
template <typename I1, typename I2>
DNA::DNA(I1 length, I2 bin_size)
    : _bins(make_bins(length, bin_size)), _length(length), _bin_size(bin_size) {
  DISABLE_WARNING_POP
  static_assert(std::is_integral<I1>::value && std::is_integral<I2>::value,
                "I1 and I2 should be an integral numeric type.");
#ifndef NDEBUG
  if (length >= std::numeric_limits<decltype(this->_length)>::max()) {
    std::runtime_error(fmt::format(
        FMT_STRING("DNA::DNA(length={}, bin_size={}): Overflow detected: unable to represent {} "
                   "in the range {}-{}"),
        length, bin_size, length, std::numeric_limits<decltype(this->_length)>::min(),
        std::numeric_limits<decltype(this->_length)>::max()));
  }
  if (length >= std::numeric_limits<decltype(this->_length)>::max()) {
    std::runtime_error(fmt::format(
        FMT_STRING("DNA::DNA(length={}, bin_size={}): Overflow detected: unable to represent {} "
                   "in the range {}-{}"),
        length, bin_size, bin_size, std::numeric_limits<decltype(this->_bin_size)>::min(),
        std::numeric_limits<decltype(this->_bin_size)>::max()));
  }
#endif
}

template <typename I1, typename I2>
void validate_params(std::string_view func_signature, I1 start, I2 end, const DNA::Bin* const bin) {
  if (start < std::numeric_limits<decltype(bin->get_start())>::min() ||
      start > std::numeric_limits<decltype(bin->get_start())>::max()) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("{}: Overflow detected: start={} cannot be represented in the range {}-{}"),
        func_signature, start, std::numeric_limits<decltype(bin->get_start())>::min(),
        std::numeric_limits<decltype(bin->get_start())>::max()));
  }
  if (end < std::numeric_limits<decltype(bin->get_end())>::min() ||
      end > std::numeric_limits<decltype(bin->get_end())>::max()) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("{}: Overflow detected: end={} cannot be represented in the range {}-{}"),
        func_signature, end, std::numeric_limits<decltype(bin->get_end())>::min(),
        std::numeric_limits<decltype(bin->get_end())>::max()));
  }
  if (start >= end) {
    throw std::logic_error(
        fmt::format(FMT_STRING("{}: start coordinate should be smaller than the end one: {} >= {}"),
                    func_signature, start, end));
  }
}

DISABLE_WARNING_PUSH
DISABLE_WARNING_CONVERSION
template <typename I>
DNA::Bin::Bin(std::size_t idx, I start, I end, const std::vector<ExtrusionBarrier>& barriers)
    : _idx(idx),
      _start(start),
      _end(end),
      _extr_barriers(std::make_unique<absl::InlinedVector<ExtrusionBarrier, 3>>(barriers.begin(),
                                                                                barriers.end())) {
  DISABLE_WARNING_POP
  static_assert(std::is_integral<I>::value, "I should be an integral numeric type.");

#ifndef NDEBUG
  validate_params(fmt::format(FMT_STRING("DNA::Bin::Bin(idx={}, start={}, end={}, const "
                                         "std::vector<ExtrusionBarrier>& barriers)"),
                              idx, start, end),
                  start, end, this);
#endif
}

DISABLE_WARNING_PUSH
DISABLE_WARNING_CONVERSION
template <typename I>
DNA::Bin::Bin(std::size_t idx, I start, I end)
    : _idx(idx), _start(start), _end(end), _extr_barriers(nullptr) {
  DISABLE_WARNING_POP
#ifndef NDEBUG
  validate_params(
      fmt::format(FMT_STRING("DNA::Bin::Bin(idx={}, start={}, end={})"), idx, start, end), start,
      end, this);
#endif
}

}  // namespace modle