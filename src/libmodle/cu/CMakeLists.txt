add_library(libmodle_cuda)

target_sources(
  libmodle_cuda
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/simulation.cu
          ${CMAKE_CURRENT_SOURCE_DIR}/simulation_device_functions.cu
          ${CMAKE_CURRENT_SOURCE_DIR}/simulation_kernels.cu
          ${CMAKE_CURRENT_SOURCE_DIR}/simulation_scheduler.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/simulation_states.cu)

target_include_directories(libmodle_cuda PUBLIC include/)

target_link_libraries(
  libmodle_cuda
  PRIVATE project_options
  PUBLIC Modle::common
         Modle::modle_internal
         CudaAPIWrappers::Runtime
         CONAN_PKG::abseil
         CONAN_PKG::fmt
         $<IF:$<BOOL:${BUILD_SHARED_LIBS}>,CUDA::cudart,CUDA::cudart_static>
         $<IF:$<BOOL:${BUILD_SHARED_LIBS}>,CUDA::curand,CUDA::curand_static>
         # Deal with std::filesystem quirkiness on older compilers
         $<$<AND:$<BOOL:${USING_LIBCXX}>,$<BOOL:${EXPLICIT_LINK_VS_FILESYSTEM_REQUIRED}>>:-lc++fs>
         $<$<AND:$<BOOL:${USING_STDLIBCXX}>,$<BOOL:${EXPLICIT_LINK_VS_FILESYSTEM_REQUIRED}>>:-lstdc++fs>)

# TODO Figure out how to disable at least -Wpedantic target_compile_options(libmodle_cuda PRIVATE -Wno-pedantic -Wall
# -Wextra)
set_target_properties(libmodle_cuda PROPERTIES CUDA_SEPARABLE_COMPILATION ON CUDA_RESOLVE_DEVICE_SYMBOLS ON)

# if(ENABLE_TESTING) add_subdirectory(test) endif()
