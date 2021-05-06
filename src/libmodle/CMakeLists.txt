add_library(libmodle INTERFACE)
add_library(Modle::Modle ALIAS libmodle)

target_sources(
  libmodle
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/dna_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/extrusion_barriers_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/extrusion_factors_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/simulation_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/simulation_correct_moves_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/simulation_detect_collisions_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/simulation_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/simulation_scheduler_impl.hpp
            ${CMAKE_CURRENT_SOURCE_DIR}/simulation_setup_impl.hpp)

target_include_directories(libmodle INTERFACE include/)

if(ENABLE_CUDA)
  add_library(libmodle_cuda STATIC)
  target_sources(
    libmodle_cuda
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/simulation_impl.cu
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include/modle/extrusion_barriers.cuh
           ${CMAKE_CURRENT_SOURCE_DIR}/include/modle/simulation.cuh)
  set_source_files_properties(${CMAKE_CURRENT_SOURCE_DIR}/include/modle/extrusion_barriers.cuh
                              ${CMAKE_CURRENT_SOURCE_DIR}/include/modle/simulation.cuh PROPERTIES LANGUAGE CUDA)
  target_include_directories(libmodle_cuda PUBLIC include/)

  target_link_libraries(
    libmodle_cuda
    PRIVATE project_warnings project_options
    PUBLIC Modle::cmatrix $<IF:$<BOOL:${BUILD_SHARED_LIBS}>,CUDA::cudart,CUDA::cudart_static>
           $<IF:$<BOOL:${BUILD_SHARED_LIBS}>,CUDA::curand,CUDA::curand_static>)
endif()

target_link_libraries(
  libmodle
  INTERFACE project_warnings
            project_options
            $<$<BOOL:${ENABLE_CUDA}>:libmodle_cuda>
            Modle::cmatrix
            Modle::io
            $<$<BOOL:${ENABLE_XOSHIRO}>:XoshiroCpp::XoshiroCpp>
            CONAN_PKG::abseil
            CONAN_PKG::boost
            CONAN_PKG::concurrentqueue
            CONAN_PKG::cpp-sort
            # Deal with std::filesystem quirkiness on older compilers
            $<$<AND:$<BOOL:${USING_LIBCXX}>,$<BOOL:${EXPLICIT_LINK_VS_FILESYSTEM_REQUIRED}>>:-lc++fs>
            $<$<AND:$<BOOL:${USING_STDLIBCXX}>,$<BOOL:${EXPLICIT_LINK_VS_FILESYSTEM_REQUIRED}>>:-lstdc++fs>)

if(ENABLE_PCH)
  target_precompile_headers(
    libmodle
    INTERFACE
    <algorithm>
    <cstddef>
    <cstdint>
    <filesystem>
    <limits>
    <memory>
    <string>
    <string_view>
    <utility>
    <vector>
    <absl/types/span.h>
    <fmt/format.h>
    [["modle/common.hpp"]]
    [["modle/config.hpp"]])
endif()

if(ENABLE_TESTING)
  add_subdirectory(test)
endif()
