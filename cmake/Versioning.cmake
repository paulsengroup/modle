# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

if(NOT EXISTS "${PROJECT_SOURCE_DIR}/.git" AND MODLE_ENABLE_GIT_VERSION_TRACKING)
  message(
    WARNING
    "-- Unable to find .git/ under \"${PROJECT_SOURCE_DIR}\". Setting -DMODLE_ENABLE_GIT_VERSION_TRACKING=OFF"
  )
  set(MODLE_ENABLE_GIT_VERSION_TRACKING OFF)
endif()

function(ConfigureVersioning input_config_folder output_config_folder)
  set(PRE_CONFIGURE_FILE "${input_config_folder}/git.cpp.in")
  set(POST_CONFIGURE_FILE "${output_config_folder}/git.cpp")

  if(MODLE_ENABLE_GIT_VERSION_TRACKING)
    include(FetchContent)
    FetchContent_Declare(
      _modle_cmake-git-version-tracking
      URL
        "${PROJECT_SOURCE_DIR}/external/cmake-git-version-tracking.20250721.tar.xz"
      URL_HASH SHA256=0c36e9ed8bb372789850a1877e68dec6131744a3d52e0f873bd20a81d43533d4
    )
    FetchContent_MakeAvailable(_modle_cmake-git-version-tracking)

    set(GIT_IGNORE_UNTRACKED ON)
    include("${_modle_cmake-git-version-tracking_SOURCE_DIR}/git_watcher.cmake")
  else()
    # Add dummy target
    add_custom_target(_modle_check_git)

    if(NOT DEFINED MODLE_GIT_RETRIEVED_STATE)
      set(MODLE_GIT_RETRIEVED_STATE false)
    endif()
    if(NOT DEFINED MODLE_GIT_HEAD_SHA1)
      set(MODLE_GIT_HEAD_SHA1 "unknown")
    endif()
    if(NOT DEFINED MODLE_GIT_IS_DIRTY)
      set(MODLE_GIT_IS_DIRTY false)
    endif()
    if(NOT DEFINED MODLE_GIT_AUTHOR_NAME)
      set(MODLE_GIT_AUTHOR_NAME "unknown")
    endif()
    if(NOT DEFINED MODLE_GIT_AUTHOR_EMAIL)
      set(MODLE_GIT_AUTHOR_EMAIL "unknown")
    endif()
    if(NOT DEFINED MODLE_GIT_COMMIT_DATE_ISO8601)
      set(MODLE_GIT_COMMIT_DATE_ISO8601 "unknown")
    endif()
    if(NOT DEFINED MODLE_GIT_COMMIT_SUBJECT)
      set(MODLE_GIT_COMMIT_SUBJECT "unknown")
    endif()
    if(NOT DEFINED MODLE_GIT_COMMIT_BODY)
      set(MODLE_GIT_COMMIT_BODY "unknown")
    endif()
    if(NOT DEFINED MODLE_GIT_DESCRIBE)
      set(MODLE_GIT_DESCRIBE "unknown")
    endif()
    if(NOT DEFINED MODLE_GIT_BRANCH)
      set(MODLE_GIT_BRANCH "unknown")
    endif()
    if(NOT DEFINED MODLE_GIT_TAG)
      set(MODLE_GIT_TAG "unknown")
    endif()

    if(NOT WIN32)
      file(LOCK "${POST_CONFIGURE_FILE}" GUARD FUNCTION)
    endif()
    configure_file("${PRE_CONFIGURE_FILE}" "${POST_CONFIGURE_FILE}" @ONLY)
  endif()

  set(PRE_CONFIGURE_FILE "${input_config_folder}/version.cpp.in")
  set(POST_CONFIGURE_FILE "${output_config_folder}/version.cpp")

  if(NOT WIN32)
    file(LOCK "${POST_CONFIGURE_FILE}" GUARD FUNCTION)
  endif()
  configure_file("${PRE_CONFIGURE_FILE}" "${POST_CONFIGURE_FILE}" @ONLY)
endfunction()
