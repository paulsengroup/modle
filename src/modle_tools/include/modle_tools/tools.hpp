#pragma once

namespace modle::tools {

// Pre-declare config structs
struct eval_config;
struct filter_barrier_config;
struct find_barrier_cluster_config;
struct noisify_config;
struct stats_config;

void eval_subcmd(const modle::tools::eval_config &c);
void filter_barriers_subcmd(const modle::tools::filter_barrier_config &c);
void noisify_subcmd(const modle::tools::noisify_config &c);
void stats_subcmd(const modle::tools::stats_config &c);

}  // namespace modle::tools
