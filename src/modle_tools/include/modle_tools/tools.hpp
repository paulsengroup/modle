// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

namespace modle::tools {

// Pre-declare config structs
struct annotate_barriers_config;
struct eval_config;
struct transform_config;

void annotate_barriers_subcmd(const annotate_barriers_config &c);
void eval_subcmd(const eval_config &c);
void transform_subcmd(const transform_config &c);

}  // namespace modle::tools
