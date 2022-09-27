// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

namespace modle::tools {

// Pre-declare config structs
struct eval_config;
struct transform_config;

void eval_subcmd(const modle::tools::eval_config &c);
void transform_subcmd(const modle::tools::transform_config &c);

}  // namespace modle::tools
