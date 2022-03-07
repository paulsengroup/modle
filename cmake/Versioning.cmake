# cmake-format: off
FetchContent_Declare(
        _cmake-git-version-tracking
        URL ${CMAKE_CURRENT_SOURCE_DIR}/external/cmake-git-version-tracking.tar.xz
        URL_HASH SHA512=8800d7ffaad4d5a86fdbbdba807bcac48534402c2b437e3aa244700b0be16961ba82ccf5cc1101d9ec44cb35fbc8f8395b072c6bb8b6e3ae8867169cd57db0a4
)
# cmake-format: on
FetchContent_MakeAvailable(_cmake-git-version-tracking)

set(GIT_IGNORE_UNTRACKED ON)
set(PRE_CONFIGURE_FILE "${CMAKE_CURRENT_SOURCE_DIR}/config/git.hpp.in")
set(POST_CONFIGURE_FILE "${CMAKE_BINARY_DIR}/src/version/include/modle/version/git.hpp")

include(${_cmake-git-version-tracking_SOURCE_DIR}/git_watcher.cmake)

# Get whether or not the working tree is dirty.
# cmake-format: off
if(GIT_IGNORE_UNTRACKED)
  set(untracked_flag "-uno")
else()
  set(untracked_flag "-unormal")
endif()
rungitcommand(status --porcelain ${untracked_flag})
if(NOT exit_code EQUAL 0)
  set(GIT_IS_DIRTY "false")
else()
  if(NOT "${output}" STREQUAL "")
    set(GIT_IS_DIRTY "true")
  else()
    set(GIT_IS_DIRTY "false")
  endif()
endif()

# Get output of git describe
rungitcommand(describe --always ${object})
if(NOT exit_code EQUAL 0)
  set(GIT_DESCRIBE "unknown")
else()
  set(GIT_DESCRIBE "${output}")
endif()


rungitcommand(tag '--sort=-v:refname' --merged)
if(NOT exit_code EQUAL 0)
  set(GIT_LAST_SEMVER_TAG "unknown")
else()
  execute_process(
    COMMAND echo "${output}"
    COMMAND grep 'v[[:digit:]]\\+'
    COMMAND head -n1
    RESULT_VARIABLE exit_code
    OUTPUT_VARIABLE GIT_LAST_SEMVER_TAG
    OUTPUT_QUIET ERROR_QUIET)
  if(NOT exit_code EQUAL 0 OR GIT_LAST_SEMVER_TAG STREQUAL "")
    set(GIT_LAST_SEMVER_TAG "unknown")
  endif()
endif()
# cmake-format: on

set(MODLE_PROJECT_VERSION_MAJOR 1)
set(MODLE_PROJECT_VERSION_MINOR 0)
set(MODLE_PROJECT_VERSION_PATCH 0)
set(MODLE_PROJECT_VERSION_SUFFIX rc.1)
set(MODLE_PROJECT_VERSION
    "${MODLE_PROJECT_VERSION_MAJOR}.${MODLE_PROJECT_VERSION_MINOR}.${MODLE_PROJECT_VERSION_PATCH}")

# cmake-format: off
if(NOT "${MODLE_PROJECT_VERSION_SUFFIX}" STREQUAL "")
  string(CONCAT MODLE_PROJECT_VERSION "${MODLE_PROJECT_VERSION}" "-${MODLE_PROJECT_VERSION_SUFFIX}")
endif()
# cmake-format: on

if(GIT_IS_DIRTY)
  string(CONCAT MODLE_PROJECT_VERSION "${MODLE_PROJECT_VERSION}" "-dirty")
endif()

if(NOT
   "v${MODLE_PROJECT_VERSION}"
   STREQUAL
   "${GIT_LAST_SEMVER_TAG}")
  string(
    REGEX
    REPLACE "-dirty$"
            ""
            MODLE_PROJECT_VERSION
            "${MODLE_PROJECT_VERSION}")
  string(CONCAT MODLE_PROJECT_VERSION "${MODLE_PROJECT_VERSION}" "-${GIT_DESCRIBE}")
  if(GIT_IS_DIRTY)
    string(CONCAT MODLE_PROJECT_VERSION "${MODLE_PROJECT_VERSION}" "-dirty")
  endif()
endif()

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/config/version.hpp.in"
               "${CMAKE_BINARY_DIR}/src/version/include/modle/version/version.hpp" @ONLY)
