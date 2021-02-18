using combine_bwig_stats

append!(ARGS, ["--sort-by", "mean_spearman,mean_pearson", "--column-labels-to-use-for-joining", "chr,start,end,bin_size,avg_lef_processivity,prob_of_block,nlefs,prob_of_bypass,target_contact_density,lef_unloader_strength", "--column-suffixes", "_pearson,_spearman", "/tmp/dummy_data/pearson_cross.tsv.gz", "/tmp/dummy_data/spearman_cross.tsv.gz"])

combine_bwig_stats.julia_main()

