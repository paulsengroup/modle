include(ExternalProject)
include(CheckLibraryExists)
find_program(MAKE_EXE NAMES gmake nmake make)
find_program(SED_EXE NAMES sed)

ExternalProject_Add(
  libBigWig
  GIT_REPOSITORY https://github.com/dpryan79/libBigWig.git
  GIT_TAG 3a756fad0c75b2c748667481f16b200331b0a8bb
  CONFIGURE_COMMAND cp "${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists_libBigWig.txt" "<SOURCE_DIR>/CMakeLists.txt"
  COMMAND ${CMAKE_COMMAND} <SOURCE_DIR> -Wno-dev
  BUILD_BYPRODUCTS ${BINARY_DIR}/libBigWig.so ${BINARY_DIR}/libBigWig.a
  # BUILD_COMMAND ${CMAKE_COMMAND} <SOURCE_DIR> -Wno-dev
  INSTALL_COMMAND ""
  # PATCH_COMMAND ${SED_EXE} -i "s/echo \"YES\"\)$/echo \"NO\")/" "/Makefile" BUILD_IN_SOURCE ON BUILD_COMMAND
  # ${MAKE_EXE} lib-static INSTALL_COMMAND "")
)

ExternalProject_Get_Property(libBigWig BINARY_DIR)
set(LIB_BIGWIG_ROOT ${BINARY_DIR})
unset(BINARY_DIR)

if(BUILD_SHARED_LIBS) # Generator expr do not work very well here
  add_library(libstraw SHARED IMPORTED)
  set_target_properties(libBigWig_lib PROPERTIES IMPORTED_LOCATION ${STRAW_ROOT}/libstraw.so)
else()
  add_library(
    libBigWig_lib
    STATIC
    IMPORTED
    include/modle/cooler.hpp)
  set_target_properties(libBigWig_lib PROPERTIES IMPORTED_LOCATION ${LIB_BIGWIG_ROOT}/libBigWig.a)
endif()

add_dependencies(libBigWig_lib libBigWig)
target_include_directories(libBigWig_lib INTERFACE ${LIB_BIGWIG_ROOT}/../)

add_library(libio)
add_library(Modle::io ALIAS libio)
target_sources(
  libio
  PRIVATE 4dn_dcic.cpp
          bed.cpp
          chr_sizes.cpp
          juicer_contacts.cpp
          bigwig.cpp)
add_dependencies(libio libBigWig_lib) # libstraw)

target_include_directories(libio PUBLIC include/)

target_link_libraries(
  libio
  PRIVATE project_warnings project_options
  PUBLIC
    Modle::cmatrix
    libBigWig_lib
    CONAN_PKG::abseil
    CONAN_PKG::boost
    CONAN_PKG::fmt
    CONAN_PKG::hdf5
    CONAN_PKG::range-v3
    CONAN_PKG::zlib
    # Deal with std::filesystem quirkiness on older compilers
    $<$<AND:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<VERSION_LESS_EQUAL:$<CXX_COMPILER_VERSION>,9.0.0>>:-lc++fs>
    $<$<AND:$<CXX_COMPILER_ID:GNU>,$<VERSION_LESS_EQUAL:$<CXX_COMPILER_VERSION>,9.1.0>>:-lstdc++fs>)

set_target_properties(libio PROPERTIES OUTPUT_NAME io)
target_compile_definitions(libio PRIVATE -DNOCURL)

if(ENABLE_TESTING)
  add_subdirectory(test)
endif()

#[[
install(
  TARGETS libio
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES
  DESTINATION include)
]]
