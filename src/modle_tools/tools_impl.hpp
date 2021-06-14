#pragma once

// IWYU pragma: private, include "modle_tools/tools.hpp"

#include <absl/container/flat_hash_map.h>  // for flat_hash_map, BitMask, operator!=
#include <absl/container/flat_hash_set.h>  // for flat_hash_set
#include <absl/meta/type_traits.h>         // for remove_reference_t
#include <absl/strings/match.h>            // for StartsWith
#include <absl/strings/str_cat.h>          // for StrCat
#include <absl/strings/str_join.h>         // for StrJoin
#include <absl/strings/string_view.h>      // for string_view
#include <absl/strings/strip.h>            // for StripPrefix, StripSuffix
#include <absl/time/clock.h>               // for Now
#include <absl/time/time.h>                // for FormatDuration, operator-, Time
#include <absl/types/span.h>               // for Span
#include <fmt/format.h>                    // for print, FMT_STRING, format
#include <fmt/ostream.h>                   // for print, formatbuf<>::int_type

#include <algorithm>                                // for fill, transform, max
#include <array>                                    // for array, array<>::value_type
#include <boost/dynamic_bitset/dynamic_bitset.hpp>  // for dynamic_bitset
#include <cassert>                                  // for assert
#include <cgranges/IITree.hpp>                      // for IITree
#include <cmath>                                    // for sqrt
#include <cstdint>                                  // for uint32_t
#include <cstdio>                                   // for stderr, stdout
#include <filesystem>                               // for path, operator<<, create_directories
#include <fstream>                                  // for size_t, ofstream, basic_ofstream, str...
#include <functional>                               // for ref
#include <initializer_list>                         // for initializer_list
#include <iterator>                                 // for insert_iterator, inserter
#include <memory>                                   // for unique_ptr, make_unique
#include <numeric>                                  // for iota
#include <sstream>                                  // for basic_stringbuf<>::int_type, basic_st...
#include <stdexcept>                                // for runtime_error
#include <string>                                   // for string, basic_string
#include <thread>                                   // for thread
#include <utility>                                  // for pair, make_pair, move
#include <vector>                                   // for vector

#include "modle/bed.hpp"            // for Parser
#include "modle/bigwig.hpp"         // for write_range, init_bigwig_file, close_...
#include "modle/contacts.hpp"       // for ContactMatrix
#include "modle/cooler.hpp"         // for Cooler, Cooler::READ_ONLY, Cooler::WR...
#include "modle/hdf5.hpp"           // for read_attribute_int
#include "modle_tools/config.hpp"   // for config
#include "modle_tools/eval.hpp"     // for Transformation, Cross, Linear, comput...
#include "modle_tools/noisify.hpp"  // for noisify_contacts
#include "modle_tools/stats.hpp"    // for compute_number_of_contacts_after_depl...

namespace modle::tools {

}  // namespace modle::tools
