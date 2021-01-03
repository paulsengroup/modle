#pragma once

#include <cstdint>
#include <string>
#include <string_view>
#include <thread>

#include "modle/contacts.hpp"

namespace modle::dcic4dn {

struct short_record {
  short_record() = default;
  std::string readID;
  std::string chr1;
  uint64_t position1;
  std::string chr2;
  uint64_t position2;
  std::string strand1;
  std::string strand2;

  [[nodiscard]] inline std::string to_string() const;
};

class Parser {
  // TODO
};

namespace converter {
inline void from_modle(const std::vector<std::string>& path_to_cmatrices,
                       std::string_view path_to_output_file, uint64_t noise_range, bool force_write,
                       double noise_stddev, uint8_t compression_lvl = 9,
                       std::size_t nthreads = std::thread::hardware_concurrency());

}  // namespace converter

}  // namespace modle::dcic4dn

#include "../../4dn_dcic_impl.hpp"