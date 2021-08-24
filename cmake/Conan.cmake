macro(run_conan)
  # Download automatically, you can also just copy the conan.cmake file
  if(NOT EXISTS "${CMAKE_BINARY_DIR}/conan.cmake")
    message(STATUS "Downloading conan.cmake from https://github.com/conan-io/cmake-conan")
    file(
      DOWNLOAD "https://raw.githubusercontent.com/conan-io/cmake-conan/v0.16.1/conan.cmake"
      "${CMAKE_BINARY_DIR}/conan.cmake"
      EXPECTED_HASH SHA256=396e16d0f5eabdc6a14afddbcfff62a54a7ee75c6da23f32f7a31bc85db23484
      TLS_VERIFY ON)
  endif()

  set(ENV{CONAN_REVISIONS_ENABLED} 1)
  include(${CMAKE_BINARY_DIR}/conan.cmake)

  conan_add_remote(
    NAME
    conancenter
    URL
    https://center.conan.io
    VERIFY_SSL
    TRUE
    INDEX
    1)

  conan_add_remote(
    NAME
    bincrafters
    URL
    https://bincrafters.jfrog.io/artifactory/api/conan/public-conan
    VERIFY_SSL
    TRUE
    INDEX
    2)

  conan_cmake_run(
    CONANFILE
    conanfile.py
    SETTINGS
    compiler.cppstd=${CMAKE_CXX_STANDARD} # For some reason setting this in conanfile doesn't work
    OPTIONS
    enable_testing=${ENABLE_TESTING}
    ENV
    "CC=${CMAKE_C_COMPILER}"
    "CXX=${CMAKE_CXX_COMPILER}"
    BASIC_SETUP
    CMAKE_TARGETS # individual targets to link to
    BUILD
    outdated
    PROFILE_AUTO
    ALL)
endmacro()
