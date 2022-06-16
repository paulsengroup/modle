# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

include(CheckCXXSourceCompiles)

set(CMAKE_REQUIRED_FLAGS, "-std=c++17")

#
# Detect stdlib on *NIX systems
#
check_cxx_source_compiles(
  "#include <ciso646>

   #if defined __GLIBCXX__
      int a = 1;
   #else
      int a = x;
   #endif

   int main(int argc, char* argv[]) {
      return 0;
   }"
  USING_STDLIBCXX)

check_cxx_source_compiles(
  "#include <ciso646>

   #if defined _LIBCPP_VERSION
      int a = 1;
   #else
      int a = x;
   #endif

   int main(int argc, char* argv[]) {
      return 0;
   }"
  USING_LIBCXX)

#
# Detect <variant> availablilty
#
check_cxx_source_compiles(
  "#include <variant>

   int main(int argc, char* argv[]) {
      std::variant<int, long> x;
      (void)x;
      return 0;
   }"
  VARIANT_AVAILABLE)

unset(CMAKE_REQUIRED_FLAGS)
