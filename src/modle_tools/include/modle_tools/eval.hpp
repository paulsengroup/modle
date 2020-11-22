#pragma once

#include "modle/contacts.hpp"
#include "modle_tools/cli.hpp"

namespace modle::tools {
ContactMatrix<uint32_t> parse_hic_matrix(const modle::tools::config &c);
}