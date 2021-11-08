// Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

namespace modle::tools {

// Pre-declare config structs
struct eval_config;
struct find_barrier_clusters_config;
struct noisify_config;
struct stats_config;
struct transform_config;

void eval_subcmd(const modle::tools::eval_config &c);
void find_barrier_clusters_subcmd(const modle::tools::find_barrier_clusters_config &c);
void noisify_subcmd(const modle::tools::noisify_config &c);
void stats_subcmd(const modle::tools::stats_config &c);
void transform_subcmd(const modle::tools::transform_config &c);

}  // namespace modle::tools
