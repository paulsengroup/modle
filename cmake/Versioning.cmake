set(MODLE_PROJECT_VERSION_MAJOR 1)
set(MODLE_PROJECT_VERSION_MINOR 0)
set(MODLE_PROJECT_VERSION_PATCH 0)
set(MODLE_PROJECT_VERSION_SUFFIX rc.1)

function(ConfigureVersioning config_folder)
  if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.git")
    # Add dummy target
    add_custom_target(_cmake-git-version-tracking)

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

    configure_file("${config_folder}/git.cpp.in" "${config_folder}/git.cpp" @ONLY)
  else()

    # cmake-format: off
    FetchContent_Declare(
            _cmake-git-version-tracking
            URL ${CMAKE_CURRENT_SOURCE_DIR}/external/cmake-git-version-tracking.tar.xz
            URL_HASH SHA512=59223346fd6199dcbbe856c9f96d0efe75fee544305ea50eddf3fc43c0ee9e789ff5b2d6aa34d97b1de6d861be5af9ae8b82059ab0e04e3045e20aad06b0a541
    )
    # cmake-format: on
    FetchContent_MakeAvailable(_cmake-git-version-tracking)

    set(GIT_IGNORE_UNTRACKED ON)
    set(PRE_CONFIGURE_FILE "${config_folder}/git.cpp.in")
    set(POST_CONFIGURE_FILE "${config_folder}/git.cpp")
    include(${_cmake-git-version-tracking_SOURCE_DIR}/git_watcher.cmake)

  endif()

  configure_file("${config_folder}/version.cpp.in" "${config_folder}/version.cpp" @ONLY)
endfunction()

configureversioning("${CMAKE_CURRENT_SOURCE_DIR}/src/config")
