#pragma once

namespace modle::tools {

struct config;  // Pre-declaration

inline void convert_subcmd(const modle::tools::config &c);
inline void eval_subcmd(const modle::tools::config &c);

}  // namespace modle::tools

#include "../../tools_impl.hpp"
