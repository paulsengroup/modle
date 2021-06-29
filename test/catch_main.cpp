#define CATCH_CONFIG_MAIN

#include <catch2/catch.hpp>

#include "modle/common/smartdir.hpp"     // IWYU pragma: keep
#include "modle/test/contacts.hpp"       // IWYU pragma: keep
#include "modle/test/correlation.hpp"    // IWYU pragma: keep
#include "modle/test/interval_tree.hpp"  // IWYU pragma: keep
#include "modle/test/libio.hpp"          // IWYU pragma: keep
#include "modle/test/libmodle.hpp"       // IWYU pragma: keep

namespace modle::test {

constexpr auto cleanup_on_exit{true};     // Useful for debugging
const SmartDir testdir{cleanup_on_exit};  // NOLINT Using auto here upsets GCC8

}  // namespace modle::test
