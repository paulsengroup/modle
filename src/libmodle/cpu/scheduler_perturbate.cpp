// clang-format off
#include "modle/simulation.hpp"
// clang-format on

#include <absl/container/fixed_array.h>          // for FixedArray
#include <absl/strings/str_cat.h>                // for StrAppend, StrCat
#include <absl/time/clock.h>                     // for Now
#include <absl/time/time.h>                      // for FormatDuration, operator-, Time
#include <absl/types/span.h>                     // for Span, MakeConstSpan
#include <fmt/compile.h>                         // for format, FMT_COMPILE
#include <fmt/format.h>                          // for vformat_to, FMT_STRING, join
#include <fmt/ostream.h>                         // for formatbuf<>::int_type
#include <moodycamel/blockingconcurrentqueue.h>  // for BlockingConcurrentQueue
#include <moodycamel/concurrentqueue.h>          // for ConsumerToken, ProducerToken
#include <spdlog/spdlog.h>                       // for info

#include <algorithm>                        // for min, max, transform
#include <array>                            // for array, array<>::value_type
#include <atomic>                           // for atomic
#include <boost/filesystem/operations.hpp>  // for exists, remove
#include <boost/filesystem/path.hpp>        // for operator<<, path
#include <cassert>                          // for assert
#include <chrono>                           // for microseconds, milliseconds
#include <cmath>                            // for round, floor, log
#include <cstdint>                          // for uint64_t
#include <cstdio>                           // for size_t, stdout
#include <exception>                        // for exception_ptr, exception, current_exception
#include <iterator>                         // for move_iterator, make_move_iterator
#include <limits>                           // for numeric_limits
#include <mutex>                            // for mutex, scoped_lock
#include <stdexcept>                        // for runtime_error
#include <string>                           // for string
#include <string_view>                      // for string_view
#include <thread_pool/thread_pool.hpp>      // for thread_pool
#include <vector>                           // for vector

#include "modle/bed.hpp"                 // for BED, BED_tree, BED_tree::at, BED_tree::c...
#include "modle/common/common.hpp"       // for contacts_t, bp_t
#include "modle/compressed_io.hpp"       // for Writer
#include "modle/contacts.hpp"            // for ContactMatrix
#include "modle/cooler.hpp"              // for Cooler, Cooler::WRITE_ONLY
#include "modle/extrusion_barriers.hpp"  // for ExtrusionBarrier
#include "modle/genome.hpp"              // for Chromosome, Genome
#include "modle/interval_tree.hpp"       // for IITree, IITree::empty, IITree::equal_range

namespace modle {
void Simulation::run_perturbate() {
  if (!this->skip_output) {  // Write simulation params to file
    assert(boost::filesystem::exists(this->path_to_output_prefix.parent_path()));  // NOLINT
    if (this->force) {
      boost::filesystem::remove(this->path_to_output_file_bedpe);
    }
  }
  // TODO Do proper error handling
  // For now we support only the case where exactly two files are specified
  assert(this->path_to_feature_bed_files.size() == 2);  // NOLINT
  const size_t task_batch_size_enq = 32;                // NOLINTNEXTLINE
  moodycamel::BlockingConcurrentQueue<TaskPW> task_queue(this->nthreads * 2, 1, 0);
  moodycamel::ProducerToken ptok(task_queue);
  std::array<TaskPW, task_batch_size_enq> tasks;

  std::mutex out_stream_mutex;
  std::mutex cooler_mutex;

  auto out_bedpe_stream = compressed_io::Writer(this->path_to_output_file_bedpe);
  auto out_task_stream = compressed_io::Writer(this->path_to_task_file);

  this->_tpool.reset(this->nthreads);
  for (uint64_t tid = 0; tid < this->nthreads; ++tid) {  // Start simulation threads
    this->_tpool.push_task([&, tid]() {
      this->perturbate_worker(tid, task_queue, out_bedpe_stream, out_stream_mutex, cooler_mutex);
    });
  }

  const auto all_deletions =
      this->path_to_deletion_bed.empty() ? this->generate_deletions() : this->import_deletions();

  size_t task_id = 0;
  size_t num_tasks = 0;
  for (auto& chrom : this->_genome) {
    if (!this->ok()) {
      this->handle_exceptions();
    }
    const auto chrom_name = std::string{chrom.name()};
    const auto& features = chrom.get_features();
    const auto& barriers = chrom.barriers();
    if (features.size() != 2 || barriers.empty()) {
      continue;  // Skip chromosomes that have less than two "kinds" of features or no barriers
    }

    if (!all_deletions.contains(chrom_name)) {
      continue;
    }

    // The idea here is that given a diagonal width D and a feature F1, in order to track all
    // possible interactions between F1 and nearby features we need to simulate at least window
    // 2*D wide centered around F1. If we now consider another feature F2 located 100 bp
    // downstream of F1, in order to track all possible interactions we need to simulate a 2*D
    // window centered around F2. Given how close F1 and F2 we basically end up simulating the
    // same region twice. This can be avoided by extending the window we are simulating by 1*D
    // left and right, increasing the window width to 4*D. We then define an active window, which
    // is a 2*D wide region centered inside the 4*D window. Now we simulate loop extrusion over a
    // 4*D region and record all possible, non-zero contacts between features falling in the 2*D
    // region. Next, we advance both windows by 1*D, and repeat the same procedure. With this
    // logic the overlap between subsequent windows is always 3*D. Deletions are performed on the
    // entire 4*D region, this is because deletions outside of the active window can still affect
    // interactions between features in the active window. The last problem we have to solve, is
    // to figure out how to ensure that we don't output the number of contacts for a pair of
    // features and a given barrier configuration twice. This is likely to happen, as the active
    // window is 2*D and windows are advanced by only 1*D at a time. The solution to this problem
    // was to use the center of the two windows as partition point. For features upstream of the
    // partition point we output all possible pairwise contacts, while for features located
    // downstream of the partition point we output contacts only if one of the feature
    // participating in the pair is upstream of the partition point. Consider the following
    // example: Given four features F1, F2, F3 and F4, where F1 and F2 are located upstream of the
    // partition point P, while F3 and F4 are located downstream of said point. All features fall
    // in a D wide region in the active window. For each distinct extrusion barrier configuration
    // we output the number of contacts for the following pairs (assumin nnz contacts):
    //
    // - F1:F2 - both are upstream of P
    // - F1:F3 - F1 is upstream of P
    // - F1:F4 - F1 is upstream of P
    // - F2:F3 - F2 is upstream of P
    // - F2:F4 - F2 is upstream of P
    //
    // - Contacts between F3 and F4 will be produced when simulating the next window, as in that
    // case they will both be located upstream of the new partition point
    //
    // It should be noted that the first and last window on a chromosome need special handling, as
    // the outer window cannot extend past the chromosomal boundaries. In these cases the active
    // window is extended such that it goes from 5'-end->3*D for the first window and from P +
    // 1*D->3'-end, where P is the partition point of the previous window.

    size_t cell_id = 0;

    // Setup the
    TaskPW base_task{};

    base_task.chrom = &chrom;
    // Initialize task with the initial window
    base_task.window_start = 0;
    base_task.window_end = 4 * this->diagonal_width;
    base_task.active_window_start = 0;
    base_task.active_window_end = 3 * this->diagonal_width;

    do {
      if (!this->ok()) {
        this->handle_exceptions();
      }
      // Find all barriers falling within the outer window. Skip over windows with 0 barriers
      if (!Simulation::map_barriers_to_window(base_task, chrom)) {
        continue;
      }

      // Find all features of type 1 falling within the outer window. Skip over windows with 0
      // features
      if (!Simulation::map_features_to_window(base_task, chrom)) {
        continue;
      }

      assert(!base_task.barriers.empty());  // NOLINT
      assert(!base_task.feats1.empty());    // NOLINT
      assert(!base_task.feats2.empty());    // NOLINT

      const auto deletions = [&]() {  // Find deletions mapping to the outer window
        const auto [first_deletion, last_deletion] =
            all_deletions.at(chrom_name).equal_range(base_task.window_start, base_task.window_end);
        return absl::MakeConstSpan(&(*first_deletion), &(*last_deletion));
      }();

      for (const auto& deletion : deletions) {
        // Add task to the current batch
        auto& t = (tasks[num_tasks++] = base_task);  // NOLINT

        // Complete task setup
        t.id = task_id++;
        t.cell_id = cell_id++;
        t.deletion_begin = deletion.chrom_start;
        t.deletion_size = deletion.chrom_end - deletion.chrom_start;
        t.window_end = std::min(t.window_end, chrom.end_pos());
        t.active_window_end = std::min(t.active_window_end, chrom.end_pos());

        // Compute the number of simulation rounds required to reach the target contact density
        if (this->target_contact_density != 0) {
          // Compute the number of pixels mapping to the outer window
          const auto npix1 = (t.window_end - t.window_start + this->bin_size - 1) / this->bin_size;
          const auto npix2 = (this->diagonal_width + this->bin_size - 1) / this->bin_size;

          t.num_target_contacts = static_cast<size_t>(std::max(
              1.0, std::round(this->target_contact_density * static_cast<double>(npix1 * npix2))));
        }

        // Compute the number of LEFs based on the window size
        t.num_lefs = static_cast<size_t>(
            std::round((static_cast<double>(t.window_end - t.window_start) / Mbp) *
                       this->number_of_lefs_per_mbp));

        // Compute the target number of epochs based on the target number of contacts
        t.num_target_epochs = t.num_target_contacts == 0UL ? this->simulation_iterations
                                                           : (std::numeric_limits<size_t>::max)();

        if (num_tasks == tasks.size()) {  // Enqueue a batch of tasks
          if (out_task_stream) {
            out_task_stream.write(fmt::format(FMT_STRING("{}\n"), fmt::join(tasks, "\n")));
          }
          auto sleep_us = 100;  // NOLINT(readability-magic-numbers)
          while (!task_queue.try_enqueue_bulk(ptok, std::make_move_iterator(tasks.begin()),
                                              num_tasks)) {
            if (!this->ok()) {
              this->handle_exceptions();
            }  // NOLINTNEXTLINE(readability-magic-numbers)
            sleep_us = std::min(100000, sleep_us * 2);
            std::this_thread::sleep_for(std::chrono::microseconds(sleep_us));
          }
          num_tasks = 0;
        }
      }
    } while (Simulation::advance_window(base_task, chrom));
  }

  // Submit any remaining task
  if (num_tasks != 0) {
    if (out_task_stream) {
      out_task_stream.write(fmt::format(FMT_STRING("{}\n"),
                                        fmt::join(tasks.begin(), tasks.begin() + num_tasks, "\n")));
    }
    while (!task_queue.try_enqueue_bulk(ptok, std::make_move_iterator(tasks.begin()), num_tasks)) {
      if (!this->ok()) {
        this->handle_exceptions();
      }  // NOLINTNEXTLINE(readability-magic-numbers)
      std::this_thread::sleep_for(std::chrono::microseconds(100));
    }
  }

  this->_end_of_simulation = true;
  out_task_stream.close();
  this->_tpool.wait_for_tasks();     // Wait on simulate_worker threads
  assert(!this->_exception_thrown);  // NOLINT
}

void Simulation::perturbate_worker(
    const uint64_t tid, moodycamel::BlockingConcurrentQueue<Simulation::TaskPW>& task_queue,
    compressed_io::Writer& out_stream, std::mutex& out_stream_mutex, std::mutex& cooler_mutex,
    const size_t task_batch_size) {
  spdlog::info(FMT_STRING("Spawning simulation thread {}..."), tid);
  moodycamel::ConsumerToken ctok(task_queue);

  absl::FixedArray<TaskPW> task_buff(task_batch_size);  // Tasks are dequeue in batch.

  Simulation::StatePW local_state{};

  try {
    while (this->ok()) {  // Try to dequeue a batch of tasks
      const auto avail_tasks = task_queue.wait_dequeue_bulk_timed(
          ctok, task_buff.begin(), task_buff.size(), std::chrono::milliseconds(10));
      // Check whether dequeue operation timed-out before any task became available
      if (avail_tasks == 0) {
        if (!this->ok()) {
          // Reached end of simulation (i.e. all tasks have been processed)
          return;
        }
        // Keep waiting until one or more tasks become available
        continue;
      }

      // Loop over new tasks
      for (const auto& task : absl::MakeConstSpan(task_buff.data(), avail_tasks)) {
        if (!this->ok()) {
          return;
        }
        assert(!task.barriers.empty());  // NOLINT
        assert(!task.feats1.empty());    // NOLINT
        assert(!task.feats2.empty());    // NOLINT

        local_state = task;  // Set simulation local_state based on task data
        local_state.contacts.resize(local_state.window_end - local_state.window_start,
                                    this->diagonal_width, this->bin_size);

        Simulation::simulate_window(local_state, out_stream, out_stream_mutex, cooler_mutex);
      }
    }
  } catch (const std::exception& e) {
    std::scoped_lock<std::mutex> l(this->_exceptions_mutex);
    this->_exceptions.emplace_back(std::make_exception_ptr(
        std::runtime_error(fmt::format(FMT_STRING("Detected an error in worker thread {}:\n{}\n{}"),
                                       tid, local_state.to_string(), e.what()))));
    this->_exception_thrown = true;

  } catch (...) {
    std::scoped_lock<std::mutex> l(this->_exceptions_mutex);
    this->_exceptions.emplace_back(std::current_exception());
    this->_exception_thrown = true;
  }
}

void Simulation::simulate_window(Simulation::StatePW& state, compressed_io::Writer& out_stream,
                                 std::mutex& out_stream_mutex, std::mutex& cooler_mutex,
                                 bool write_contacts_to_cooler) const {
  spdlog::info(FMT_STRING("Processing {}[{}-{}]; outer_window=[{}-{}]; deletion=[{}-{}];"),
               state.chrom->name(), state.active_window_start, state.active_window_end,
               state.window_start, state.window_end, state.deletion_begin,
               state.deletion_begin + state.deletion_size);
  // Simulating a window consists in generating all valid combinations of barriers after deleting
  // a portion of the DNA from the window. The idea here is that given that biologically relevant
  // boundaries are often marked by multiple extrusion barriers, deleting one barrier at a time is
  // not very interesting, while by deleting regions spanning several kbs we are able to delete
  // clusters of barriers, which should produce a more dramatic change in the number of contacts.
  // Deletions are performed according to the following algorithm:
  //
  // Going from the first barrier at the 5-end of the window and moving towards the last barrier
  // at the 3-'end of the window do:
  // - Delete a region going from the BP + 1 - DS to BP, where BP is the position of the extrusion
  //   barrier and DS is the deletion size. This will produce a novel extrusion barrier
  //   configuration
  // - Simulate loop extrusion on the current window given the above barrier configuration
  // - Move to the next barrier and start over from step #1
  //
  // Thew idea here is that if there's a cluster of barriers going from BP1 until BP10, and
  // BP10 - BP1 < DS, when BP points to BP10, the deletion will include every barrier in the
  // cluster.

  const auto write_to_stdout = out_stream.path().empty();

  std::string out_buffer;
  std::string barrier_str_buff;
  // Figure out whether we are processing the first or last window and compute the partition
  // point
  // const auto first_window = state.active_window_start == state.chrom->start_pos();
  const auto last_window = state.active_window_end == state.chrom->end_pos();
  const auto partition_point = state.active_window_start + this->diagonal_width;

  // Generate barrier configuration
  size_t num_active_barriers = 0;
  state.barrier_tmp_buff.resize(state.barriers.size());
  std::transform(state.barriers.begin(), state.barriers.end(), state.barrier_tmp_buff.begin(),
                 [&](const auto& barrier) {
                   if (barrier.pos() >= state.deletion_begin &&
                       barrier.pos() < state.deletion_begin + state.deletion_size) {
                     // Do not change ( -> {. We don't want to construct the object using the
                     // initializer list
                     return ExtrusionBarrier(barrier.pos(), 0.0, 1.0,
                                             barrier.blocking_direction_major());
                   }
                   ++num_active_barriers;
                   return barrier;
                 });

  state.barriers = absl::MakeConstSpan(state.barrier_tmp_buff);

  // Resize and reset state buffers
  state.resize_buffers();
  state.reset_buffers();
  state.contacts.resize(state.window_end - state.window_start, this->diagonal_width,
                        this->bin_size);
  state.contacts.reset();

  Simulation::simulate_one_cell(state);

  if (out_stream) {  // Output contacts for valid pairs of features
    for (auto i = 0UL; i < state.feats1.size(); ++i) {
      const auto& feat1 = state.feats1[i];
      const auto feat1_abs_center_pos = (feat1.chrom_start + feat1.chrom_end + 1) / 2;
      for (auto j = i; j < state.feats2.size(); ++j) {
        const auto& feat2 = state.feats2[j];
        const auto feat2_abs_center_pos = (feat2.chrom_start + feat2.chrom_end + 1) / 2;

        // Don't output contacts when both feat1 and feat2 are located downstream of the partition
        // point. This does not apply when we are processing the last window
        if (!last_window && feat1.chrom_start >= partition_point &&
            feat2.chrom_start >= partition_point) {
          continue;
        }

        // Convert absolute positions to relative positions, as the contact matrix does not refer
        // to an entire chromosome, but only to a 4*D wide region
        const auto feat1_rel_center_pos = feat1_abs_center_pos - state.window_start;
        const auto feat2_rel_center_pos = feat2_abs_center_pos - state.window_start;

        const auto feat1_rel_bin = feat1_rel_center_pos / this->bin_size;
        const auto feat2_rel_bin = feat2_rel_center_pos / this->bin_size;

        const auto contacts = state.contacts.get(feat1_rel_bin, feat2_rel_bin, this->block_size);
        if (contacts == 0) {  // Don't output entries with 0 contacts
          continue;
        }

        // Compute the absolute bin coordinates. This are used to compute the positions shown in
        // the BEDPE file
        const auto feat1_abs_bin = (feat1_abs_center_pos + this->bin_size - 1) / this->bin_size;
        const auto feat2_abs_bin = (feat2_abs_center_pos + this->bin_size - 1) / this->bin_size;

        // Generate the name field. The field will be "none;none" in case both features don't have
        // a name
        const auto name = absl::StrCat(feat1.name.empty() ? "none" : feat1.name, ";",
                                       feat2.name.empty() ? "none" : feat2.name);
        absl::StrAppend(  // Append a BEDPE record to the local buffer
            &out_buffer,
            fmt::format(
                FMT_COMPILE(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n"),
                feat1.chrom, feat1_abs_bin * bin_size, (feat1_abs_bin + 1) * bin_size, feat2.chrom,
                feat2_abs_bin * bin_size, (feat2_abs_bin + 1) * bin_size, name, contacts,
                feat1.strand, feat2.strand, state.deletion_begin,
                state.deletion_begin + state.deletion_size, num_active_barriers,
                state.barriers.size(), state.active_window_start, state.active_window_end,
                state.window_start, state.window_end));
      }
    }

    // Write the buffer to the appropriate stream
    if (!out_buffer.empty()) {
      if (std::scoped_lock l(out_stream_mutex); write_to_stdout) {
        fmt::print(stdout, FMT_STRING("{}"), out_buffer);
      } else {
        out_stream.write(out_buffer);
      }
      out_buffer.clear();
    }
  }

  if (write_contacts_to_cooler) {
    const auto file_name = boost::filesystem::path(fmt::format(
        FMT_STRING("{}_{:06d}_{}_window_{}-{}_deletion_{}-{}.cool"),
        this->path_to_output_prefix.string(), state.id, state.chrom->name(), state.window_start,
        state.window_end, state.deletion_begin, state.deletion_begin + state.deletion_size));

    if (!this->force && boost::filesystem::exists(file_name)) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("Refusing to overwrite output file {}. Pass --force to overwrite."),
          file_name));
    }

    spdlog::info(FMT_STRING("Writing contacts for {} to file {}..."), state.chrom->name(),
                 file_name);
    const auto t0 = absl::Now();
    {
      std::scoped_lock l(cooler_mutex);
      auto c = cooler::Cooler(file_name, cooler::Cooler::WRITE_ONLY, this->bin_size,
                              this->_genome.chromosome_with_longest_name().name().size());
      for (const auto& chrom : this->_genome) {
        if (&chrom == state.chrom) {
          c.write_or_append_cmatrix_to_file(state.contacts, state.chrom->name(), state.window_start,
                                            state.window_end, state.chrom->size());
          continue;
        }
        assert(chrom.contacts_ptr() == nullptr);  // NOLINT
        c.write_or_append_empty_cmatrix_to_file(chrom.name(), chrom.start_pos(), chrom.end_pos(),
                                                chrom.size());
      }
    }
    spdlog::info(FMT_STRING("DONE writing file {} in {}"), file_name,
                 absl::FormatDuration(absl::Now() - t0));
  }
}

}  // namespace modle
