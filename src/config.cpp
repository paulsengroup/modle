#include "modle/config.hpp"

namespace modle {
Config::Config(std::string_view path_to_chr_bed, std::string_view config_file)
    : _path_to_chr_bed(path_to_chr_bed), _config_file(config_file) {}

}  // namespace modle