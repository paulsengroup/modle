using compute_statistics_for_bwigs

append!(ARGS, ["--chr-subranges-bed", "/tmp/dummy_data/chr_subranges.bed", "--extrusion-barriers-bed", "/tmp/dummy_data/extr_barriers.bed", "/tmp/dummy_data/sample.txt"])
compute_statistics_for_bwigs.julia_main()

