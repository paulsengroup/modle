include(CheckCXXSourceCompiles)

check_cxx_source_compiles(
  "#include <ciso646>

     #if defined __GLIBCXX__ && (__GLIBCXX__ < 20190503 || _GLIBCXX_RELEASE < 9)
        int a = 1;
     #else
        int a = x;
     #endif

     int main(int argc, char* argv[]) {
        return 0;
     }"
  STDLIBCXX_REQUIRES_LSTDCXXFS_LINK_FLAG)

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

if(STDLIBCXX_REQUIRES_LSTDCXXFS_LINK_FLAG OR LIBCXX_REQUIRES_LCXXFS_LINK_FLAG)
  set(EXPLICIT_LINK_VS_FILESYSTEM_REQUIRED TRUE)
endif()

if(NOT STD_FILESYSTEM_AVAILABLE)
  set(EXPLICIT_LINK_VS_FILESYSTEM_REQUIRED FALSE)
endif()

unset(STDLIBCXX_REQUIRES_LSTDCXXFS_LINK_FLAG)
unset(LIBCXX_REQUIRES_LCXXFS_LINK_FLAG)
