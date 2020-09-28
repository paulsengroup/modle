#ifndef MODLE_CONFIG_HPP
#define MODLE_CONFIG_HPP

#include <string>

namespace modle {
class Config {
 public:
 explicit Config(std::string_view path_to_chr_bed, std::string_view config_file = "");

 uint32_t bin_size{20'000};
 private:
  std::string _config_file;
  std::string _path_to_chr_bed;
};

}  // namespace modle

#endif  // MODLE_CONFIG_HPP