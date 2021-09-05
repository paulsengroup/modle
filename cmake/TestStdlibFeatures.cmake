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
# Detect <charconv> availability
#
check_cxx_source_compiles(
  "#include <charconv>

   int main(int argc, char* argv[]) {
       int n;
       char* p;
       (void)std::from_chars(p, p, n);

       return 0;
   }"
  CHARCONV_INT_AVAILABLE)

check_cxx_source_compiles(
  "#include <charconv>
   #include <string>

   int main(int argc, char* argv[]) {
       double n;
       char* p;
       (void)std::from_chars(p, p, n);

       return 0;
   }"
  CHARCONV_FP_AVAILABLE)

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

#
# Detect whether <filesystem> requires explicit linking
#
check_cxx_source_compiles(
  "#include <ciso646>

     #if defined _LIBCPP_VERSION && _LIBCPP_VERSION < 9000
        int a = 1;
     #else
        int a = x;
     #endif

     int main(int argc, char* argv[]) {
        return 0;
     }"
  LIBCXX_REQUIRES_LCXXFS_LINK_FLAG)

check_cxx_source_compiles(
  "#include <filesystem>

     std::filesystem::path p(\"foo\");

     int main(int argc, char* argv[]) {
        return 0;
     }"
  STD_FILESYSTEM_AVAILABLE)

unset(CMAKE_REQUIRED_FLAGS)

if(STDLIBCXX_REQUIRES_LSTDCXXFS_LINK_FLAG OR LIBCXX_REQUIRES_LCXXFS_LINK_FLAG)
  set(EXPLICIT_LINK_VS_FILESYSTEM_REQUIRED ON)
endif()

if(NOT STD_FILESYSTEM_AVAILABLE)
  set(EXPLICIT_LINK_VS_FILESYSTEM_REQUIRED OFF)
endif()

unset(STDLIBCXX_REQUIRES_LSTDCXXFS_LINK_FLAG)
unset(LIBCXX_REQUIRES_LCXXFS_LINK_FLAG)
