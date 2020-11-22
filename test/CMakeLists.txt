# automatically enable catch2 to generate ctest targets
if(CONAN_CATCH2_ROOT_DEBUG)
  include(${CONAN_CATCH2_ROOT_DEBUG}/lib/cmake/Catch2/Catch.cmake)
else()
  include(${CONAN_CATCH2_ROOT}/lib/cmake/Catch2/Catch.cmake)
endif()

add_executable(catch_main catch_main.cpp)
target_link_libraries(
  catch_main
  PUBLIC CONAN_PKG::catch2
  PRIVATE project_warnings
          project_options
          Modle::cmatrix_test
          Modle::correlation_test
          Modle::io_test
          # This allows Catch2 to discover test units by preventing unused symbols from being discarded by the linker
          -Wl,--whole-archive
          $<TARGET_FILE:Modle::cmatrix_test>
          $<TARGET_FILE:Modle::correlation_test>
          $<TARGET_FILE:Modle::io_test>
          -Wl,--no-whole-archive)

# automatically discover tests that are defined in catch based test files you can modify the unittests. TEST_PREFIX to
# whatever you want, or use different for different binaries
catch_discover_tests(
  catch_main
  EXTRA_ARGS
  -s
  --reporter=xml
  --out=tests.xml)
