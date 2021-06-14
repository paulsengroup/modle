#pragma once

namespace modle::tools {

struct config;  // Pre-declaration

void eval_subcmd(const modle::tools::config &c);
void filter_barriers_subcmd(const modle::tools::config &c);
void noisify_subcmd(const modle::tools::config &c);
void stats_subcmd(const modle::tools::config &c);

}  // namespace modle::tools
