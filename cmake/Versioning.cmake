set(MODLE_PROJECT_VERSION_MAJOR 1)
set(MODLE_PROJECT_VERSION_MINOR 0)
set(MODLE_PROJECT_VERSION_PATCH 0)
set(MODLE_PROJECT_VERSION_SUFFIX rc.4)

function(ConfigureVersioning input_config_folder output_config_folder)
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.git")
    # cmake-format: off
    FetchContent_Declare(
            _cmake-git-version-tracking
            URL ${CMAKE_CURRENT_SOURCE_DIR}/external/cmake-git-version-tracking.tar.xz
            URL_HASH SHA512=84d990063503c75eaedab5645cc21ca86738d84011ef0d7b630853a1bcab37693dc941333fc1e31ad421f625bb4d8b5f908bb2fa2618eff9bf1370e2ba76a1ba
    )
    # cmake-format: on
    FetchContent_MakeAvailable(_cmake-git-version-tracking)

    set(GIT_IGNORE_UNTRACKED ON)
    set(PRE_CONFIGURE_FILE "${input_config_folder}/git.cpp.in")
    set(POST_CONFIGURE_FILE "${output_config_folder}/git.cpp")
    include(${_cmake-git-version-tracking_SOURCE_DIR}/git_watcher.cmake)
  else()
    # Add dummy target
    add_custom_target(check_git)

    set(GIT_RETRIEVED_STATE false)
    set(GIT_HEAD_SHA1 "unknown")
    set(GIT_IS_DIRTY false)
    set(GIT_AUTHOR_NAME "unknown")
    set(GIT_AUTHOR_EMAIL "unknown")
    set(GIT_COMMIT_DATE_ISO8601 "unknown")
    set(GIT_COMMIT_SUBJECT "unknown")
    set(GIT_COMMIT_BODY "unknown")
    set(GIT_DESCRIBE "unknown")
    set(GIT_BRANCH "unknown")
    set(GIT_TAG "unknown")

    configure_file("${input_config_folder}/git.cpp.in" "${output_config_folder}/git.cpp" @ONLY)
  endif()

  configure_file("${input_config_folder}/version.cpp.in" "${output_config_folder}/version.cpp" @ONLY)
endfunction()

configureversioning("${CMAKE_CURRENT_SOURCE_DIR}/src/config" "${CMAKE_CURRENT_BINARY_DIR}/src/config")
