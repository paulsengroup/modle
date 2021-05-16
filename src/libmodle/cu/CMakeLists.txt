add_library(libmodle_cuda INTERFACE)

target_sources(
  libmodle_cuda
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/simulation.cu
            #${CMAKE_CURRENT_SOURCE_DIR}/simulation_device_functions.cu
            ${CMAKE_CURRENT_SOURCE_DIR}/simulation_kernels.cu
            ${CMAKE_CURRENT_SOURCE_DIR}/simulation_scheduler.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/simulation_states.cu)

target_include_directories(libmodle_cuda INTERFACE include/)

target_link_libraries(
  libmodle_cuda
  INTERFACE project_warnings
            project_options
            Modle::config
            Modle::modle_io_internal
            CudaAPIWrappers::Runtime
            CONAN_PKG::abseil
            CONAN_PKG::fmt

            $<IF:$<BOOL:${BUILD_SHARED_LIBS}>,CUDA::cudart,CUDA::cudart_static>
            $<IF:$<BOOL:${BUILD_SHARED_LIBS}>,CUDA::curand,CUDA::curand_static>
            # Deal with std::filesystem quirkiness on older compilers
            $<$<AND:$<BOOL:${USING_LIBCXX}>,$<BOOL:${EXPLICIT_LINK_VS_FILESYSTEM_REQUIRED}>>:-lc++fs>
            $<$<AND:$<BOOL:${USING_STDLIBCXX}>,$<BOOL:${EXPLICIT_LINK_VS_FILESYSTEM_REQUIRED}>>:-lstdc++fs>)

# TODO Move these options into a separate target
target_compile_options(libmodle_cuda INTERFACE -Wno-pedantic)

# if(ENABLE_TESTING) add_subdirectory(test) endif()
