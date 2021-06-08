add_library(modle_common INTERFACE)

add_library(Modle::common ALIAS modle_common)

target_sources(
  modle_common
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include/modle/common/common.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/include/modle/common/config.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/include/modle/common/genextreme_value_distribution.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/include/modle/common/random_sampling.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/include/modle/common/suppress_compiler_warnings.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/include/modle/common/utils.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/utils_impl.hpp)

target_link_libraries(
  modle_common
  INTERFACE project_options
            project_warnings
            $<$<BOOL:${ENABLE_XOSHIRO}>:XoshiroCpp::XoshiroCpp>
            CONAN_PKG::abseil
            CONAN_PKG::fmt
            CONAN_PKG::boost)

target_include_directories(modle_common INTERFACE include/)
