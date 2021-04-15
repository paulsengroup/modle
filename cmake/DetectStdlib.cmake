include(CheckCXXSourceCompiles)

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
