include(ExternalProject)
include(CheckLibraryExists)
find_program(MAKE_EXE NAMES gmake nmake make)
find_program(SED_EXE NAMES sed)

ExternalProject_Add(
  libBigWig
  GIT_REPOSITORY https://github.com/dpryan79/libBigWig.git
  GIT_TAG 3a756fad0c75b2c748667481f16b200331b0a8bb
  CONFIGURE_COMMAND ""
  PATCH_COMMAND ${SED_EXE} -i "s/echo \"YES\"\)$/echo \"NO\")/" "<SOURCE_DIR>/Makefile"
  BUILD_IN_SOURCE ON
  BUILD_COMMAND ${MAKE_EXE} lib-static
  INSTALL_COMMAND "")

add_library(libBigWig_static STATIC IMPORTED)
add_dependencies(libBigWig_static libBigWig)
ExternalProject_Get_Property(libBigWig SOURCE_DIR)
set(LIB_BIGWIG_ROOT ${SOURCE_DIR})
unset(SOURCE_DIR)
target_include_directories(libBigWig_static INTERFACE ${LIB_BIGWIG_ROOT}/../)
set_target_properties(libBigWig_static PROPERTIES IMPORTED_LOCATION ${LIB_BIGWIG_ROOT}/libBigWig.a)

add_library(libio)
add_library(Modle::io ALIAS libio)
target_sources(
  libio
  PRIVATE bed.cpp
          chr_sizes.cpp
          juicer_contacts.cpp
          bigwig.cpp)
add_dependencies(libio libBigWig_static)

target_include_directories(
  libio
  PUBLIC include/
  PRIVATE libBigWig_static)

target_link_libraries(
  libio
  PRIVATE project_warnings project_options m # -lm
  PUBLIC Modle::cmatrix
         libBigWig_static
         CONAN_PKG::abseil
         CONAN_PKG::fmt
         CONAN_PKG::range-v3
         CONAN_PKG::zlib)

set_target_properties(libio PROPERTIES OUTPUT_NAME parser)

if(ENABLE_TESTING)
  add_subdirectory(test)
endif()
