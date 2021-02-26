#!/usr/bin/env bash

set -e

# ./testing_settings ~/github/modle/test/data/

wd="$1"
container="$wd/modle_v0.0.1-x86-64.sif"
chr_sizes="$wd/test2.chr.sizes"
base_out_dir="$wd/"
p_of_lef_bypass=0.25
n_barr=40
n_lefs=10
sampl_interval=125
lef_unloader_strn=0.5

cpus="$(nproc)"

bin_sizes=(1000 2000 5000 10000)
iters=(15000000 7500000 3000000 1500000)
p_of_block=(0.5 0.6 0.7 0.75 0.8 0.85 0.9 0.95 0.975 0.99 1.0)
avg_proc=(10000 25000 50000 75000 100000 200000)

# https://unix.stackexchange.com/questions/103920/parallelize-a-bash-for-loop

for pr in "${avg_proc[@]}"; do
  for pb in "${p_of_block[@]}"; do
    for i in "${!bin_sizes[@]}"; do

      bin_size="${bin_sizes[$i]}"
      iter="${iters[$i]}"
      printable_bin_size="$((bin_size / 1000))kb"
      printable_proc="$(bc -l <<< "scale=1; $pr/1000")kb"
      printable_iter="$(bc -l <<< "scale=1; $iter/1000000")e6"
      out_dir="$base_out_dir/$(basename "$chr_sizes" '.chr.sizes')/pblock_$pb-proc_$printable_proc-binsize_$printable_bin_size-iters_$printable_iter/"
      # echo \
      mkdir -p "$out_dir"
      # echo \
      rm -f "$out_dir/log"
      command time -v \
      singularity exec "$container" modle                                             \
                                    -c "$chr_sizes"                                   \
                                    --avg-lef-processivity "$pr"                      \
                                    --probability-of-lef-bypass "$p_of_lef_bypass"    \
                                    --number-of-randomly-generated-barriers "$n_barr" \
                                    --number-of-randomly-generated-lefs "$n_lefs"     \
                                    --randomize-contact-sampling-interval             \
                                    --contact-sampling-interval "$sampl_interval"     \
                                    --hard-stall-multiplier "$lef_unloader_strn"      \
                                    --force                                           \
                                    -b "$bin_size"                                    \
                                    --number-of-iterations "$iter"                    \
                                    --probability-of-barrier-block "$pb"              \
                                    -o "$out_dir" 2> "$out_dir/log" &

      if [[ $(jobs -r -p | wc -l) -ge $cpus ]]; then
        # now there are $cpus jobs already running, so wait here for any job
        # to be finished so there is a place to start next one.
        wait -n
    fi
    done;
  done;
done;
