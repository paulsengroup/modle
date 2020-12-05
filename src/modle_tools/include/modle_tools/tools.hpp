#pragma once

namespace modle::tools {

struct config;  // Pre-declaration

void convert_subcmd(const modle::tools::config &c);
void eval_subcmd(const modle::tools::config &c);

}  // namespace modle::tools
