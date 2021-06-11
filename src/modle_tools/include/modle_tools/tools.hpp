#pragma once

namespace modle::tools {

struct config;  // Pre-declaration

inline void eval_subcmd(const modle::tools::config &c);
inline void filter_barriers_subcmd(const modle::tools::config &c);
inline void noisify_subcmd(const modle::tools::config &c);
inline void stats_subcmd(const modle::tools::config &c);

}  // namespace modle::tools

#include "../../tools_impl.hpp"  // IWYU pragma: keep
