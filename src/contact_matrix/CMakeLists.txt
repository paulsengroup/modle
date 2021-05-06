add_library(cmatrix INTERFACE)
add_library(Modle::cmatrix ALIAS cmatrix)

target_sources(cmatrix INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/contacts_impl.hpp)

target_include_directories(cmatrix INTERFACE include/)

if(ENABLE_CUDA)
  add_library(cmatrix_cuda INTERFACE)

  target_sources(cmatrix_cuda INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include/modle/contacts.cuh)
  set_source_files_properties(${CMAKE_CURRENT_SOURCE_DIR}/include/modle/contacts.cuh PROPERTIES LANGUAGE CUDA)

  target_include_directories(cmatrix_cuda INTERFACE include/)

  target_link_libraries(cmatrix_cuda INTERFACE project_warnings project_options
                                               $<IF:$<BOOL:${BUILD_SHARED_LIBS}>,CUDA::cudart,CUDA::cudart_static>)
endif()

target_link_libraries(
  cmatrix
  INTERFACE project_warnings
            project_options
            Modle::utils
            $<$<BOOL:${ENABLE_XOSHIRO}>:XoshiroCpp::XoshiroCpp>
            $<$<BOOL:${ENABLE_CUDA}>:cmatrix_cuda>
            CONAN_PKG::boost)

if(ENABLE_PCH)
  target_precompile_headers(
    cmatrix
    INTERFACE
    <cstdint>
    <utility>
    <vector>
    <boost/dynamic_bitset/dynamic_bitset.hpp>
    [["modle/contacts.hpp"]])
endif()

if(ENABLE_TESTING)
  add_subdirectory(test)
endif()
